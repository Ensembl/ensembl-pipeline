package Bio::EnsEMBL::Pipeline::RunnableDB::Blastz_m;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blastz_m;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::SeqFeature;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);  

  if(!defined $self->seqfetcher) {
    # will look for pfetch in $PATH - change this once PipeConf up to date
    my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch; 
    $self->seqfetcher($seqfetcher);
  }  
  
  return $self; 
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for exonerate from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my( $self) = @_; 
  print STDERR "Fetching input \n";
  
  $self->throw("No input id") unless defined($self->input_id);

  my $input_id  = $self->input_id;
  

  my $chr;
  my $chrstart;
  my $chrend;
 
  my $contig;
 
  if ($input_id =~ /(\S+)\.(\S+)-(\S+)/) {
    $chr      = $1;
    $chrstart = $2;
    $chrend   = $3;

    $contig  = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$chrstart,$chrend);
  } else {
    $contig = $self->db->get_RawContigAdaptor->fetch_by_name($input_id);
  }

  $self->slice($contig);

  my $genseq = $contig->get_repeatmasked_seq();

  $self->{_genseq} = $genseq;
  
  my @db;

  if (!defined $contig || !defined($self->analysis->db)) {
    $self->throw("Can't run blat if no sequence [$contig] or database [" . $self->analysis->db ."] exists");
  }

  my $executable =  $self->analysis->program_file();


  if ( -d $self->analysis->db) {
    print "Found database dir " . $self->analysis->db . "\n";
    my $db = $self->analysis->db;
    my $files = `ls -1 $db/*`;

    @db = split(/\n/,$files);
  } else {
    push(@db,$self->analysis->db);
  }

  foreach my $db (@db) {
    my $blastz = new Bio::EnsEMBL::Pipeline::Runnable::Blastz_m(-query         => $genseq,
							    -program       => $executable,
							    -database      => $db,
							    -analysis      => $self->analysis);
    
    $self->runnable($blastz);
  }
}


=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the exonerate analysis, producing Bio::EnsEMBL::Gene predictions
    Returns :   Nothing, but $self{_output} contains the predicted genes.
    Args    :   None

=cut

sub run {
  my ($self) = @_;
  
  $self->throw("Can't run - no runnable objects") unless defined($self->runnable);
  
  foreach my $run ($self->runnable) {
    $run->run;
  }


}
sub slice {
  my ($self,$slice) = @_;

  if (defined($slice)) {
    $self->{_slice} = $slice;
  }
  return $self->{_slice};
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->{_output} into $self->db
    Returns :   1
    Args    :   None

=cut

sub write_output {

  my($self) = @_;
  
  my @features = $self->output();
  my $db       = $self->db;

  my $fa = $db->get_DnaAlignFeatureAdaptor;

  foreach my $output (@features) {
    $output->contig($self->slice);
    $output->attach_seq($self->slice);

    $output->analysis($self->analysis);

    if ($self->slice->isa("Bio::EnsEMBL::Slice")) {
    my @mapped = $output->transform;

    if (@mapped == 0) {
      $self->warn("Couldn't map $output - skipping");
      next;
    }
    if (@mapped == 1 && $mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")) {
      $self->warn("$output seems to be on a gap - something bad has happened ...");
      next;
    }

    $fa->store(@mapped);
    } else {
      $fa->store($output);
    }
  }
}


1;
