#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Exofish

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Exofish;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Exofish;
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);



@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
  my($self) = @_;
  
  $self->throw("No input id") unless defined($self->input_id);
  
  my $contig    = $self->db->get_RawContigAdaptor->fetch_by_name($self->input_id);
  
  if (@$PIPELINE_REPEAT_MASKING) {
    $contig = $contig->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING);
  }

  $self->query($contig);

  my %parameters;
  $parameters{-query} = $self->query;
  $parameters{-database} = $self->analysis->db_file;
  $parameters{-program} = $self->analysis->program_file;
  $parameters{-options} = $self->analysis->parameters;

  my $run = Bio::EnsEMBL::Pipeline::Runnable::Exofish->new(%parameters);
  $self->runnable($run);

  return 1;
}


sub write_output {
  my($self) = @_;
  
  my @features = $self->output();
  my $db       = $self->db;
  
  my $fa = $db->get_DnaAlignFeatureAdaptor;
  
  foreach my $output (@features) {
    
    $output->contig($self->query);    
    $output->analysis($self->analysis);
    
    if ($self->query->isa("Bio::EnsEMBL::Slice")) {
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
    }
    else {
      $fa->store($output);      
    }
  }
}


1;
