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

Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -analysis   => $analysis
                                             );
    $obj->fetch_input();
    $obj->run();

    my @newfeatures = $obj->output();
    
    $obj->write_output();

=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::ExonerateArray;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for exonerate from the query_file
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my( $self) = @_;
  print STDERR "Fetching input \n";
  
  my $input_id = $self->input_id;
  my $analysis = $self->analysis;
  my $exonerate = $analysis->program;
  my $version = $analysis->program_version;
  my $options;
  my $query_file = $MICROARRAY_INPUT_DIR.$input_id;
  my $target_dir =$MICROARRAY_TARGET_DIR;
  
  ###Runnable::ExonerateArray take a array of query_seq_obj, so it's need to be generated here###

  my @query_seqs;

  my $in = Bio::SeqIO->newFh(
			     -FILE => $query_file,
			     -FORMAT => 'Fasta',
			    );

  while (my $seq = <$in>) {
    push (@query_seqs, $seq);
  }

  $self->total_query_seq(scalar @query_seqs);

  # prepare runnable
  
  $self->throw("Can't run Exonerate without both query and target sequences") 
    unless (defined($query_file) && defined($target_dir));
  
  print "exonerate is '$exonerate', target_dir is $target_dir, query_file is $ query_file\n";
  my $executable =  $BIN_DIR."/".$exonerate."-".$version;
  
  $self->throw("No exonerate executable") unless defined($executable);

  #my $target_file = $target_dir . "*";###exonerate-0.8.2 can use both file and dir
  my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::ExonerateArray(
								      '-db'           => $self->db,
								      '-database'  => $target_dir,
								      '-query_seqs'   => \@query_seqs,
								      '-query_type'   => 'dna',
								      '-target_type'  => 'dna',
								      '-exonerate'    => $executable,
								      '-options'      => $options,
								      '-analysis'     => $analysis,
								     );
  $self->runnables($runnable);
}


=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the exonerate analysis, producing Bio::EnsEMBL::MiscFeature
    Returns :   Nothing, but $self{_output} contains the features
    Args    :   None

=cut

sub run {
  my ($self) = @_;
  
  $self->throw("Can't run - no runnable objects") unless defined($self->runnable);
  
  foreach my $runnable ($self->runnables){
    $runnable->run();
    my @out = $runnable->output();
    my $match = $runnable->output_match_count;
    $self->output(\@out);
    $self->output_match_count($match);
  }
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns the contents of $self->{_output}, which holds predicted features.
    Returns :   Array of Bio::EnsEMBL::DnaDnaAlignFeature
    Args    :   None

=cut

sub output {

  my ($self, $output) = @_;
  if ($output){
    push(@{$self->{_output}}, @$output);
  }
  return @{$self->{_output}};
}

sub output_match_count {

  my ($self, $match) = @_;
  if ($match){
    $self->{_match} = $match;
  }
  return $self->{_match};
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->{_output} into $self->dbobj
    Returns :   1
    Args    :   None

=cut

sub write_output {

  my($self) = @_;
  
  my @misc_features = $self->output(); 
  my $match = $self->output_match_count;
  
  my $mfa = $self->db->get_MiscFeatureAdaptor();
  $mfa->store( @misc_features );
  
  my ($tot_pass_ids,$only_25,$only_24,$both,$total_query_seq,%done);

  foreach my $q_id (keys %{$match}) {
    $tot_pass_ids++;
    if ($match->{$q_id}->{'full_match_count'} and !$match->{$q_id}->{'mis_match_count'}) {
      $only_25 += $match->{$q_id}->{'full_match_count'};
    }
    elsif ($match->{$q_id}->{'full_match_count'} and $match->{$q_id}->{'mis_match_count'}) {
      $both +=$match->{$q_id}->{'full_match_count'};
      $both +=$match->{$q_id}->{'mis_match_count'};
    }
    elsif ($match->{$q_id}->{'mis_match_count'} and !$match->{$q_id}->{'fullmatch_count'}) {
      $only_24 += $match->{$q_id}->{'mis_match_count'};
    }
  }
  
  
  my $ratio_25 = $only_25/$tot_pass_ids;
  my $ratio_24 = $only_24/$tot_pass_ids;
  my $ratio_both = $both/$tot_pass_ids;
  $total_query_seq += $self->total_query_seq;
  printf "total_query_seq is $total_query_seq, total pass ids is $tot_pass_ids, with 25 bases exact match is $only_25 (%.2f), with 24 bases exact match is $only_24 (%.2f) and both is $both (%.2f)\n", $ratio_25,$ratio_24,$ratio_both;
  
  return 1;
}


############################################################
#
# get/set methods
#
############################################################

sub query_file {
  
  my ($self,$file) = @_;
  if( defined $file) {
  $self->{'_query_file'} = $file;
}
  return $self->{'_query_file'};
}

############################################################

sub target_file {
  
  my ($self,$file) = @_;
  if( defined $file) {
  $self->{'_target_file'} = $file;
}
  return $self->{'_target_file'};
}


############################################################

sub exonerate {
  my ($self, $location) = @_;
  if ($location) {
    $self->throw("Exonerate not found at $location: $!\n") unless (-e $location);
    $self->{_exonerate} = $location ;
  }
  return $self->{_exonerate};
}

############################################################

sub options {
  my ($self, $options) = @_;
  if ($options) {
    $self->{_options} = $options ;
  }
  return $self->{_options};
}

############################################################

sub db {
  my ($self, $db) = @_;
  if ($db) {
    $self->{_db} = $db ;
  }
  return $self->{_db};
}

############################################################

sub input_id {
  my ($self, $input_id) = @_;
  if ($input_id) {
    $self->{_input_id} = $input_id ;
  }
  return $self->{_input_id};
}

############################################################

sub total_query_seq {
  my ($self, $total_query_seq) = @_;
  if ($total_query_seq) {
    $self->{_total_query_seq} = $total_query_seq ;
  }
  return $self->{_total_query_seq};
}


############################################################
    
sub runnables{
  my ($self, $runnable) = @_;
  if($runnable) {
    unless($runnable->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
      $self->throw("$runnable is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
    push @{$self->{_runnable}}, $runnable;
  }
  return @{$self->{_runnable}};
}


1;
