package Bio::EnsEMBL::Pipeline::Runnable::Blastz_m;

use strict;

use Bio::SeqIO;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::DnaDnaAlignFeature;

use vars qw(@ISA);

# Object preamble - inherits from Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;

  my $self = bless {}, $class;

  my ($query,$program,$database,$analysis) = 
    $self->_rearrange([qw(QUERY PROGRAM DATABASE ANALYSIS)], @args);

  $self->throw("Must pass in query sequence and database file") if (! defined $query || ! defined $database);

  $self->query   ($query);
  $self->database($database);
  $self->analysis($analysis);

  $self->program($self->find_executable($program)) if ($program);

  return $self;
}

=head2 analysis

 Title   : analysis
 Usage   : $obj->analysis($newval)
 Function: 
 Example : 
 Returns : value of analysis
 Args    : newvalue (optional)


=cut

sub analysis{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'analysis'} = $value;
    }
    return $obj->{'analysis'};

}

sub query {
  my ($obj,$value) = @_;
  if( defined $value) {
    if (! ref $value || ! $value->isa('Bio::PrimarySeqI')) {
      $obj->throw("$value is not a PrimarySeqI object. Cannot throw");
    }
    my $selfseq = Bio::PrimarySeq->new(-display_id => $value->id,
				       -seq => $value->seq);
    if ($selfseq->length == 0) {
      $obj->throw("attempting to bl2seq seemingly 0 length sequence!");
    }
    $obj->{'_seq1'} = $selfseq;
  }
  return $obj->{'_seq1'};
}


sub database {
  my ($self,$value) = @_;

  if (defined($value)) {
     if ( -e $value) {
       $self->{_database} = $value;
     } else {
       $self->throw("Database file [$value] doesn't exist. Exiting\n");
     }
  }
  return $self->{_database};
}

=head2 program

    Title   :   program
    Usage   :   $obj->program('/usr/local/ensembl/bin/bl2seq');
    Function:   Get/set method for the location of the bl2seq executable
    Returns :   string
    Args    :   string

=cut

sub program {
  my ($self, $location) = @_;

  if ($location) {
    $self->throw("executable not found at $location: $!\n") unless (-e $location && -x $location);
    $self->{'_program'} = $location ;
  }
  return $self->{'_program'};
}

=head2 run

 Title   : run
 Usage   : $self->run();
 Function: run bl2seq program (NCBI) on two previously defined sequences
 Example :
 Returns : 1
 Args    : None


=cut

sub run {
  my ($self) = @_;

  # Soft mask the sequence N - n

  my $seq = $self->query->seq;

  $seq =~ tr/N/n/;
  $self->query->seq($seq);

  my $dbname = $self->analysis->db;

  $dbname =~ s/.*\/(.*?)/$1/;

  my $query  = "/tmp/" . $self->query->id . ".blastz." . $dbname . ".fa";
  my $output = "/tmp/" . $self->query->id . ".blastz." . $dbname . ".out";

  my $seqio = new Bio::SeqIO(-file => ">$query",-format => 'fasta');
  $seqio->write_seq($self->query);
  $seqio->close;

  $self->run_analysis($query,$output);

  $self->parse_results($output);

  # delete sequence files
  unlink($query);
  unlink($self->results);

  return 1;
}

sub run_analysis {
  my ($self,$query,$outfile) = @_;

  my $param = $self->analysis->parameters;

  my $command = "blastz $query " . $self->database . " T=1 > " . $outfile;

  print STDERR ("Running blastz [$command]\n");

  my $status = system("$command");

}
  
sub parse_results { 
  my ($self,$file) = @_;
  
  open BLASTZ, $file || $self->throw("Couldn't open file ".$file.", $!\n");

  my $score;
  my $start;
  my $qid;
  my $qstart;
  my $qend;
  my $qstrand;
  my $hid;
  my $hstart;
  my $hend;
  my $hstrand;
  my $len;
  my $pid;

  while (<BLASTZ>) {
    print $_;
    chomp;
    
    my @f = split(' ',$_);

    if ($f[0] eq 'h') {
      my $idline  = <BLASTZ>;
      my $hidline = <BLASTZ>;

      my @qid = split(' ',$idline);
      $qid = $qid[0];

      my @hid = split(' ',$hidline);
      $hid = $hid[0];

      if ($hid eq "\"REVCOMP") {
	$hstrand = -1;
	$hid = $hid[2];
      }

      $qid =~ s/\"//g;
      $hid =~ s/\"//g;

    }

    if ($f[0] eq 'l') {
	$qstart   = $f[1];
	$qend     = $f[3];
	$hstart   = $f[2];
	$hend     = $f[4];
	$pid      = $f[5];
        
	my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start      =>   $qstart,
						   -end         =>   $qend,
						   -seqname     =>   $qid,
						   -strand      =>   1,
						   -score       =>   $score,
                                                   -percent_id  =>   $pid,
						   -analysis => $self->analysis,
						  );
	
	my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $hstart,
						   -end         =>   $hend,
						   -seqname     =>   $hid,
						   -strand      =>   $hstrand,
						   -score       =>   $score,
                                                   -percent_id  =>   $pid,
						  -analysis => $self->analysis);
	
	#create featurepair

        eval {
	my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
						 -feature2 => $feat2,
						);
	
	my $df = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => [$fp]);


	$self->output($df);
        };
        if ($@) {
          print STDERR "Error creating feature from :\n" . $feat1->gffstring . "\n" . $feat2->gffstring . "\n";
        }
      } elsif ($f[0] eq 's') {
        $score = $f[1];
      }
    }
  close(BLASTZ);
}
  


=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   None

=cut

sub output {
  my ($self,$val) = @_;

    if (!defined($self->{_output})) {
      $self->{_output} = [];
    }

  if (defined($val)) {
    push(@{$self->{_output}},$val);

  }

  return @{$self->{'_output'}};
}

=head2 workdir

 Title   : workdir
 Usage   : $obj->workdir($newval)
 Function: 
 Example : 
 Returns : value of workdir
 Args    : newvalue (optional)


=cut

sub workdir{
   my ($self,$value) = @_;
   if( defined $value) {
       $self->{'_workdir'} = $value;
   }
   return $self->{'_workdir'};
}

1;
