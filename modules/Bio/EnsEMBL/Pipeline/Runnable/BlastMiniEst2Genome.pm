#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome->new('-genomic'    => $genseq,
									 '-seqfetcher' => $seqfetcher
									 '-blastdb'    => $blastdb);

    $obj->run

    my @features = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Analysis::MSPcrunch;
use Bio::PrimarySeqI;
use Bio::Tools::Blast;
use Bio::SeqIO;
use Bio::DB::RandomAccessI;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-GENOMIC       => $genomicseq,
			   -BLASTDB       => $blastdbname,
                           -SEQFETCHER    => $sf);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome object
    Returns :   A Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome object
    Args    :   -genomic:    Bio::PrimarySeqI object (genomic sequence)
                -blastdb:    Path to Blastable DB
                -seqfetcher  Bio::DB::RandomAccessI object
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_idlist'} = []; #create key to an array of feature pairs
    
    my( $genomic, $blastdb, $seqfetcher) = $self->_rearrange([qw(GENOMIC
								 BLASTDB
								 SEQFETCHER)],
							     @args);
       
    $self->throw("No genomic sequence input")           
      unless defined($genomic);
    $self->throw("[$genomic] is not a Bio::PrimarySeqI") 
      unless $genomic->isa("Bio::PrimarySeqI");
    $self->genomic_sequence($genomic) if defined($genomic);

    $self->throw("No blastdb specified") 
      unless defined($blastdb);
    $self->blastdb($blastdb) if defined($blastdb);

    $self->throw("No seqfetcher provided")           
      unless defined($seqfetcher);
    $self->throw("[$seqfetcher] is not a Bio::DB::RandomAccessI") 
      unless $seqfetcher->isa("Bio::DB::RandomAccessI");
    $self->seqfetcher($seqfetcher) if defined($seqfetcher);
    
    return $self; 
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

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::DB:RandomAccessI
    Args    :   Bio::DB:RandomAccessI

=cut

sub seqfetcher {
    my( $self, $value ) = @_;    
    if ($value) {
      #need to check if passed sequence is Bio::EnsEMBL::Pipeline::SeqFetcherI object
      $self->throw("Input isn't a Bio::DB::RandomAccessI") unless $value->isa("Bio::DB::RandomAccessI");
      $self->{'_seqfetcher'} = $value;
    }
    return $self->{'_seqfetcher'};
}

=head2 blastdb

    Title   :   blastdb
    Usage   :   $self->blastdb($seq)
    Function:   Get/set method for blastdb
    Returns :   
    Args    :   path to blastdb

=cut

sub blastdb {
    my( $self, $blastdb ) = @_;    
    if ($blastdb) {

      # check for presence of relevant files .csq, .nhd, .ntb
      $self->throw("Can't find blastfile [$blastdb]\n") unless -e $blastdb;
      my $csq = $blastdb . ".csq";
      $self->throw("Can't find .csq file [$csq]\n") unless -e $csq;
      my $nhd = $blastdb . ".nhd";
      $self->throw("Can't find .nhd file [$nhd]\n") unless -e $nhd;
      my $ntb = $blastdb . ".ntb";
      $self->throw("Can't find .ntb file [$ntb]\n") unless -e $ntb;

      $self->{'_blastdb'} = $blastdb;
    }
    return $self->{'_blastdb'};
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

    if (!defined($id)) {
	$self->throw("No id input to parse_Header");
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
  Function: Runs blast vs dbEST, and runs a MiniEst2Genome runnable for each appropriate set of blast hits
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;

    # filter ESTs using exonerate
    my @exonerate_res = $self->run_exonerate();

    print STDERR "**got " . scalar(@exonerate_res) . " hits\n";
    
    my %exonerate_ests;
    foreach my $res(@exonerate_res ) {
      my $seqname = $res->hseqname;       #gb|AA429061.1|AA429061
      $seqname =~ s/\S+\|(\S+)\|\S+/$1/;
      print "exonerate: " . $seqname . "\n";

      # score cutoff 500 for exonerate ...
      if($res->score > 500) {
	$exonerate_ests{$seqname} = 1;
      }
    }

    # run blast on genomic seq vs dbEST
    my @blastres = $self->run_blast(keys %exonerate_ests);

    # get list of hits - 
    my %esthash;
    foreach my $result(@blastres) {
      my $seqname = $result->hseqname;       #gb|AA429061.1|AA429061
      $seqname =~ s/\S+\|(\S+)\|\S+/$1/;
      $result->hseqname($seqname);
      
      # score cutoff? percentID? any EST with any hit > ... gets used
      if($result->score > 500 || defined ($esthash{$seqname}) ) {
	push(@{$esthash{$seqname}},$result);
      }
    }
    
    foreach my $id(keys %esthash) {
      print STDERR "id: $id "  . scalar(@{$esthash{$id}}) . "hits\n";

      # make a set of features per EST
      # should we do some ordering? deal with strands?
      my @features = @{$esthash{$id}};
 
      # make MiniEst2Genome runnables
      
      my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome('-genomic'  => $self->genomic_sequence,
								     '-features' => \@features,
								     '-seqfetcher' => $self->seqfetcher);

      # run runnable
      $e2g->run;
      
      # sort out output
      my @f = $e2g->output;
      
      foreach my $f (@f) {
	#      print(STDERR "PogAligned output is $f " . $f->seqname . " " . $f->start . "\t" . $f->end . "\t" . $f->score .  "\n");
      }
      
      push(@{$self->{'_output'}},@f);
    }
}

sub run_exonerate {
  my ($self) = @_;
  my @res;

  my $exr = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
							    '-exonerate' => "/work2/gs2/gs2/bin/exonerate-0.3d",
							    '-genomic'   => $self->genomic_sequence,
							    '-est'      => $self->blastdb
							   );
  
  $exr->run;
  my @res = $exr->output;

  return @res;
  
}

sub run_blast {

    my ($self, @ests) = @_;

    my $genomic = $self->genomic_sequence;
    my @seq = $self->get_Sequences(@ests);
    my $blastdb = $self->make_blast_db(@seq);

    # tmp files
    my $blastout = $self->get_tmp_file("/tmp/","blast","tblastn_dbest.msptmp");
    my $seqfile  = $self->get_tmp_file("/tmp/","seq","fa");

    my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
				-file   => ">$seqfile");

    $seqio->write_seq($genomic);
    close($seqio->_filehandle);

    my $command  = "wublastn $blastdb $seqfile B=500 -hspmax 1000  2> /dev/null |MSPcrunch -d - >  $blastout";

    print (STDERR "Running command $command\n");
    my $status = system( $command );

    print("Exit status of blast is $status\n");
    open (BLAST, "<$blastout") 
        or $self->throw ("Unable to open Blast output $blastout: $!");    
    if (<BLAST> =~ /BLAST ERROR: FATAL:  Nothing in the database to search!?/)
    {
        print "nothing found\n";
        return;
    }

    # process the blast output
    my @pairs;

    eval {
	my $msp = new Bio::EnsEMBL::Analysis::MSPcrunch(-file => $blastout,
							-type => 'DNA-DNA',
							-source_tag => 'e2g',
							-contig_id => $self->genomic_sequence->id,
							);


	@pairs = $msp->each_Homol;
	
	foreach my $pair (@pairs) {
	    my $strand1 = $pair->feature1->strand;
	    my $strand2 = $pair->feature2->strand;
	    
	    print STDERR "***" . $pair->seqname . " " . $pair->hseqname . " " . $pair->score . "\n\n";

	    $pair->invert;
	    $pair->feature2->strand($strand2);
	    $pair->feature1->strand($strand1);
	    $pair->hseqname($genomic->id);
	    $pair->invert;
	    $self->print_FeaturePair($pair);
	}
    };
    if ($@) {
	$self->warn("Error processing msp file for " . $genomic->id . " [$@]\n");
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

    my $blastfile = $self->get_tmp_file('/tmp/','blast','fa');
    my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
			       -file   => ">$blastfile");

    print STDERR "Blast db file is $blastfile\n";

    foreach my $seq (@seq) {
	print STDERR "Writing seq " . $seq->id ."\n";
	$seqio->write_seq($seq);
    }

    close($seqio->_filehandle);

    my $status = system("pressdb $blastfile");
    print (STDERR "Status from pressdb $status\n");

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
        print STDERR ("mrna feature $seq is not a Bio::PrimarySeq or Bio::Seq\n") 
                                    unless ($seq->isa("Bio::PrimarySeq") ||
                                            $seq->isa("Bio::Seq"));
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
                $self->warn ("Cleaned up ".$seq->display_id
                   ." for blast : $invalidCharCount invalid chars \n");
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
    
    print(STDERR "Sequence id :  is [$id]\n");

    eval{
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };
    if ($@){
      $self->throw("Problem with seqfetcher [$id]: [$@]\n");
    }
	
    if(!defined($seq)){
      $self->throw("Could not find sequence for [$id]");
    }

    print (STDERR "Found sequence for $id [" . $seq->length() . "]\n");

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


