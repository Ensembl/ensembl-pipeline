

#!/usr/local/bin/perl

#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new
    ('-genomic'    => $genseq,
     '-features'   => $features,
     '-seqfetcher' => $seqfetcher,
     '-trim'       => 0);
    
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::PrimarySeqI;
use Bio::Tools::Blast;
use Bio::Tools::BPlite;
use Bio::SeqIO;
use Bio::DB::RandomAccessI;

use Data::Dumper;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_TBLASTN
					);

require "Bio/EnsEMBL/Pipeline/pipeConf.pl";

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_idlist'} = []; #create key to an array of feature pairs
  
  my( $genomic, $ids, $seqfetcher, 
      $trim, $endbias) = $self->_rearrange([qw(GENOMIC
                                               IDS
                                               SEQFETCHER
                                               TRIM
                                               ENDBIAS)],
                                           @args);
  
  $self->throw("No genomic sequence input")           unless defined($genomic);
  $self->throw("[$genomic] is not a Bio::PrimarySeqI") unless $genomic->isa("Bio::PrimarySeqI");
  $self->genomic_sequence($genomic) if defined($genomic);
  
  $self->endbias($endbias) if defined($endbias);
  
  $self->throw("No seqfetcher provided")           
    unless defined($seqfetcher);

#  $self->throw("[$seqfetcher] is not a Bio::DB::RandomAccessI") 
#    unless $seqfetcher->isa("Bio::DB::RandomAccessI");
  $self->seqfetcher($seqfetcher) if defined($seqfetcher);
  
  if (defined($ids)) {
    if (ref($ids) eq "ARRAY") {
#      print "Ids @$ids\n";
      push(@{$self->{'_idlist'}},@$ids);
    } else {
      $self->throw("[$ids] is not an array ref.");
    }
  }
  
  if (defined($trim)) {
    $self->trim($trim);
  }
  
  return $self; # success - we hope!
}

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 endbias

    Title   :   endbias
    Usage   :   $self->endbias($endbias)
    Function:   Get/set method for genewise endbias
    Returns :   
    Args    :   

=cut

sub endbias {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{'_endbias'} = $arg;
    }

    if (!defined($self->{'_endbias'})) {
      $self->{'_endbias'} = 0;
    }    

    return $self->{'_endbias'};
}

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::EnsEMBL::Pipeline::SeqFetcher object
    Args    :   Bio::EnsEMBL::Pipeline::SeqFetcher object

=cut

sub seqfetcher {
  my( $self, $value ) = @_;    
  if ($value) {
    #need to check if passed sequence is Bio::DB::RandomAccessI object
#    $value->isa("Bio::DB::RandomAccessI") || $self->throw("Input isn't a Bio::DB::RandomAccessI");
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
}


=head2 get_all_FeatureIds

  Title   : get_all_FeatureIds
  Usage   : my @ids = get_all_FeatureIds
  Function: Returns an array of all distinct feature hids 
  Returns : @string
  Args    : none

=cut

sub get_Ids {
    my ($self) = @_;

    if (!defined($self->{'_idlist'})) {
        $self->{'_idlist'} = [];
    }
    return @{$self->{'_idlist'}};
}


=head2 parse_Header

  Title   : parse_Header
  Usage   : my $newid = $self->parse_Header($id);
  Function: Parses different sequence headers
  Returns : string
  Args    : none

=cut

sub parse_Header {
    my ($self,$id) = @_;

    if (!defined($id)) {        $self->throw("No id input to parse_Header");
    }

    my $newid = $id;

    if ($id =~ /^(.*)\|(.*)\|(.*)/) {
        $newid = $2;
        $newid =~ s/(.*)\..*/$1/;
        
    } elsif ($id =~ /^..\:(.*)/) {
        $newid = $1;
    }
    $newid =~ s/ //g;
    return $newid;
}


=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs est2genome on each distinct feature id
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;

    my @ids = $self->get_Ids;

    my @features = $self->blast_ids(@ids);
    my @newfeatures;

    my %scorehash;

    unless (@features) {
        print STDERR "Contig has no associated features\n";
        return;
    }

    foreach my $f (@features) {
      if (!defined $scorehash{$f->hseqname} || $f->score > $scorehash{$f->hseqname})  {
        $scorehash{$f->hseqname} = $f->score;
      }
    }

    my @forder = sort { $scorehash{$b} <=> $scorehash{$a}} keys %scorehash;

    my $mg      = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise('-genomic'    => $self->genomic_sequence,
                                                                     '-features'   => \@features,
                                                                     '-seqfetcher' => $self->seqfetcher,
                                                                     '-forder'     => \@forder,
                                                                     '-endbias'    => $self->endbias);

    $mg->minirun;
    
    my @f = $mg->output;

    push(@{$self->{'_output'}},@f);

}

sub blast_ids {
    my ($self,@ids) = @_;

    my @seq         = $self->get_Sequences(@ids);
    my @valid_seq   = $self->validate_sequence(@seq);
    
    my @blastseqs   = ($self->genomic_sequence);
    
    my $blastdb     = $self->make_blast_db(@blastseqs);
    my @newfeatures;

    foreach my $seq (@valid_seq) {
        my @tmp = $self->run_blast($seq,$blastdb);
        push(@newfeatures,@tmp);
    }

    unlink $blastdb;
    unlink $blastdb.".csq";
    unlink $blastdb.".nhd";
    unlink $blastdb.".ntb";

    return @newfeatures;
}

sub run_blast {

    my ($self,$seq,$db) = @_;
    my $tmpdir = $::pipeConf{'nfstmp.dir'};
    if(!defined $tmpdir || $tmpdir eq ''){
      $tmpdir = '/tmp';
    }

    my $blastout = $self->get_tmp_file($tmpdir,"blast","out");
    my $seqfile  = $self->get_tmp_file($tmpdir,"seq","fa");
    my @pairs;

    my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
                                -file   => ">$seqfile");

    $seqio->write_seq($seq);
    close($seqio->_filehandle);
    my $tblastn = $GB_TBLASTN;
    # default to tblastn
    if(!defined $tblastn || $tblastn eq ''){
      $tblastn = 'wutblastn';
    }
    my $command  = "$tblastn $db $seqfile B=500 -hspmax 1000 -hitdist=40 > $blastout";

   #print (STDERR "Running command $command\n");
    my $status = system($command );

#    print("Exit status of blast is $status\n");

    my $report = new Bio::Tools::BPlite('-file'=>$blastout);

    while(my $sbjct = $report->nextSbjct){
      while(my $hsp = $sbjct->nextHSP){

	
        # strands
        my $strand = 1;
        if($hsp->subject->strand != $hsp->query->strand){
          $strand = -1;
        }
	
        my $genomic = new Bio::EnsEMBL::SeqFeature(
                                                   -start   => $hsp->subject->start,
                                                   -end     => $hsp->subject->end,
                                                   -strand  => $strand,
                                                   -seqname => $hsp->subject->seqname,
                                                   -score   => $hsp->score,
                                                  );
        # munging protein seqname as BPlite is giving it back like O95793 (577 letters)
        my $protname = $hsp->query->seqname;
        $protname =~ s/^(\S+).+/$1/;
        my $protein = new Bio::EnsEMBL::SeqFeature(
                                                   -start   => $hsp->query->start,
                                                   -end     => $hsp->query->end,
                                                   -strand  => 1,
                                                   -seqname => $protname,
                                                   -score   => $hsp->score,
                                                   );
        my $featurepair = new Bio::EnsEMBL::FeaturePair(
                                                        -feature1 => $genomic,
                                                        -feature2 => $protein
                                                       );

        push (@pairs, $featurepair);
      }
    }
  

    unlink $blastout;
    unlink $seqfile;

    return @pairs;
}

sub print_FeaturePair {
    my ($self,$pair) = @_;

    print STDERR $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . $pair->score . "\t" .
        $pair->strand . "\t" . $pair->hseqname . "\t" . $pair->hstart . "\t" . $pair->hend . "\t" . $pair->hstrand . "\n";
}

sub make_blast_db {
    my ($self,@seq) = @_;

    my $tmpdir = $::pipeConf{'nfstmp.dir'};
    if(!defined $tmpdir || $tmpdir eq ''){
      $tmpdir = '/tmp';
    }
    my $blastfile = $self->get_tmp_file($tmpdir,'blast','fa');
    my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
                               -file   => ">$blastfile");

#    print STDERR "Blast db file is $blastfile\n";

    foreach my $seq (@seq) {
        $seqio->write_seq($seq);
    }

    close($seqio->_filehandle);

    my $status = system("pressdb $blastfile");

    return $blastfile;
}

sub get_tmp_file {
    my ($self,$dir,$stub,$ext) = @_;

    
    if ($dir !~ /\/$/) {
        $dir = $dir . "/";
    }

#    $self->check_disk_space($dir);

    my $num = int(rand(10000));
    my $file = $dir . $stub . "." . $num . "." . $ext;

    while (-e $file) {
        $num = int(rand(10000));
        $file = $stub . "." . $num . "." . $ext;
    }                   
    
    return $file;
}
    
sub get_Sequences {
    my ($self,@ids) = @_;

    my @seq;

    foreach my $id (@ids) {
        my $seq = $self->get_Sequence($id);

        if (defined($seq) && $seq->length > 0) {
            push(@seq,$seq);
        } else {
            print STDERR "Invalid sequence for $id - skipping\n";
        }
    }

    return @seq;

}

sub validate_sequence {
    my ($self,@seq) = @_;
    my @validated;
    foreach my $seq (@seq)
    {
  #      print STDERR ("mrna feature $seq is not a Bio::PrimarySeq or Bio::Seq\n") 
  #                                  unless ($seq->isa("Bio::PrimarySeq") ||
  #                                          $seq->isa("Bio::Seq"));

        my $sequence = $seq->seq;
        if ($sequence !~ /[^acgtn]/i)
        {
            push (@validated, $seq);
        }
        else 
        {
            $_ = $sequence;
            my $len = length ($_);
            my $invalidCharCount = tr/bB/xX/;

            if ($invalidCharCount / $len > 0.05)
            {
                $self->warn("Ignoring ".$seq->display_id()
                    ." contains more than 5% ($invalidCharCount) "
                    ."odd nucleotide codes ($sequence)\n Type returns "
                    .$seq->moltype().")\n");
            }
            else
            {
#               print STDERR ("Cleaned up ".$seq->display_id
 #                  ." for blast : $invalidCharCount invalid chars \n");
                $seq->seq($_);
                push (@validated, $seq);
            }
        }
    } 
    return @validated;  
}

=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id)
  Function: Fetches sequences with id $id
  Returns : Bio::PrimarySeq
  Args    : none

=cut
    
sub get_Sequence {
    my ($self,$id) = @_;
    my $seqfetcher = $self->seqfetcher;
    my $seq;

    if (!defined($id)) {
      $self->warn("No id input to get_Sequence");
    }  
    
#    print(STDERR "Sequence id :  is [$id]\n");

    eval {
      #print STDERR "BlastMiniGenewise: getting sequence for $id\n";
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };
    if($@) {
      $self->warn("Problem fetching sequence for id [$id] $@\n");
      return undef;
    }
    
    if(!defined($seq)){
      $self->warn("Could not find sequence for [$id]");
    }

#    print (STDERR "Found sequence for $id [" . $seq->length() . "]\n");

    return $seq;
}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
    my ($self) = @_;
    if (!defined($self->{'_output'})) {
        $self->{'_output'} = [];
    }
    return @{$self->{'_output'}};
}


sub trim {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{'_trim'} = $arg;
  }
  return $self->{'_trim'};
}

1;
