
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
                                                                         '-seqfetcher' => $seqfetcher,
                                                                         '-queryseq'   => $queryseq,
                                                                         '-threshold'  => $threshold);

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

use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;

use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_INPUTID_REGEX
					);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-GENOMIC       => $genomicseq,
                           -QUERYSEQ      => $queryseq,
                           -SEQFETCHER    => $sf,
                           -THRESHOLD     => $threshold);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome object
    Returns :   A Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome object
    Args    :   -genomic:    Bio::PrimarySeqI object (genomic sequence)
                -queryseq:   Either path to file containing query seqs or reference to aray of Bio::Seq
                -seqfetcher  Bio::DB::RandomAccessI object
                -threshold   minimum e-value for blast hits - defaults to 10 e-60
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_idlist'} = []; #create key to an array of feature pairs
    #print "@args\n";
    my( $genomic, $queryseq, $seqfetcher, $threshold) = $self->_rearrange([qw(GENOMIC
                                                                               QUERYSEQ
                                                                               SEQFETCHER
                                                                               THRESHOLD)],
                                                                           @args);
         
    $self->throw("No genomic sequence input")           
      unless defined($genomic);
    $self->throw("[$genomic] is not a Bio::PrimarySeqI") 
      unless $genomic->isa("Bio::PrimarySeqI");
    $self->genomic_sequence($genomic) if defined($genomic);

    $self->throw("No queryseq specified") 
      unless defined($queryseq);
    $self->queryseq($queryseq) if defined($queryseq);

    $self->throw("No seqfetcher provided")           
      unless defined($seqfetcher);
    $self->throw("[$seqfetcher] is not a Bio::DB::RandomAccessI") 
      unless $seqfetcher->isa("Bio::DB::RandomAccessI");
    $self->seqfetcher($seqfetcher) if defined($seqfetcher);

    if(defined $threshold){
      $self->blast_threshold($threshold);
    }
    else{
      $self->blast_threshold(10e-60);
    } 

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

=head2 blast_threshold

    Title   :   blast_threshold
    Usage   :   $self->blast_threshold($num)
    Function:   Get/set method for blast_threshold
    Returns :   numerical value
    Args    :   numerical value

=cut

sub blast_threshold {
    my( $self, $value ) = @_;    
    if ($value) {
        $self->{'_blast_threshold'} = $value;
    }
    return $self->{'_blast_threshold'};
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

=head2 queryseq

    Title   :   queryseq
    Usage   :   $self->queryseq($seq)
    Function:   Get/set method for queryseq
    Returns :   
    Args    :   name of a file containing query seq(s), OR reference to an array of Bio::Seq

=cut

sub queryseq {
    my( $self, $queryseq ) = @_;   
    if ($queryseq) { 
      if (ref($queryseq) eq 'ARRAY') {

        # I'm not at all sure this is right
        my $time = time; chomp($time);
        my $estfile = "/tmp/estfile_.$$.$time.fn";
        $self->queryfilename($estfile);

        foreach my $est(@$queryseq) {
          $est->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        }

        $self->{'_query_sequences'} = $queryseq;
      }
      else {
        # it's a filename - check the file exists
        $self->throw("[$queryseq] : file does not exist\n") unless -e $queryseq;
        $self->queryfilename($queryseq);
        $self->{'_query_sequences'} = $queryseq;
    }
  }
  
  #NB ref to an array of Bio::Seq
  return $self->{'_query_sequences'};

  }

=head2 queryfilename

    Title   :   queryfilename
    Usage   :   $self->queryfilename($filename)
    Function:   Get/set method for queryfilename
    Returns :   
    Args    :   

=cut

sub queryfilename {
  my ($self, $queryfilename) = @_;
  $self->{'_queryfilename'} = $queryfilename if ($queryfilename);
  return $self->{'_queryfilename'};
}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs blast vs input seqs, and runs a MiniEst2Genome runnable for each appropriate set of blast hits
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;

    # filter ESTs using blast
    my @blast_res = $self->run_blast();

    print STDERR "BlastMiniEst2Genome->run ***Have " . scalar(@blast_res) . " blast hits.\n";
    
    my %blast_ests;
    my @feat;
    foreach my $res(@blast_res ) {
      my $seqname = $res->hseqname;       #gb|AA429061.1|AA429061

      # will this break for cDNAs with sensible names?
      $seqname =~ s/\S+\|(\S+)\|\S+/$1/;

      $res->hseqname($seqname);

      push(@{$blast_ests{$seqname}}, $res);
    }

    foreach my $est(keys %blast_ests) {
      my @features = @{$blast_ests{$est}};
 
      # make MiniEst2Genome runnables
      
      my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome(
                           '-genomic'    => $self->genomic_sequence,
                           '-features'   => \@features,
                           '-seqfetcher' => $self->seqfetcher);

      # run runnable
      $e2g->run;
      
      # sort out output
      my @f = $e2g->output;
      print "Features returned from e2g: " . scalar @f . "\n";;      
      push(@{$self->{'_output'}},@f);
    }

    return 1;
}

sub run_blast {
  
  my ($self) = @_;
  
  my $estseqs = $self->queryseq;
  my @valid_estseqs = $self->validate_sequence(@{$estseqs});

  my $blastdb     = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(
				  -sequences => \@valid_estseqs,
				  -type      => 'DNA');
  
  $blastdb->run;
  
  my $dbname = $blastdb->dbname;
  
  
  my $genomic = $self->genomic_sequence;
  
  my $blast = new Bio::EnsEMBL::Pipeline::Runnable::Blast(
				  -query    => $genomic,
				  -program  => 'wublastn',
				  -database => $blastdb->dbfile,
				  -filter => 1,
				 );
  my $regex;
  
  if(($GB_INPUTID_REGEX)&&
     ($self->genomic_sequence->id =~ /$GB_INPUTID_REGEX/)){
    $regex = $GB_INPUTID_REGEX;
  }elsif ($self->genomic_sequence->id =~ /^(.*)\|(.*)\|(.*)/) {
    $regex = '^.*\|(.*)\|.*';
  } elsif ($self->genomic_sequence->id =~ /^..\:(.*)/) {
    $regex = '^..\:(.*)';
  }else {
    $regex = '^(\w+)\s+';
  }

  $blast->add_regex($dbname, $regex);
  $blast->run;
  
  my @blast_features = $blast->output;
  
  $blastdb->remove_index_files;
  unlink $blastdb->dbfile;
  
  return @blast_features;
}


sub print_FeaturePair {
    my ($self,$pair) = @_;

    print STDERR $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . $pair->score . "\t" .
        $pair->strand . "\t" . $pair->hseqname . "\t" . $pair->hstart . "\t" . $pair->hend . "\t" . $pair->hstrand . "\n";
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

1;
