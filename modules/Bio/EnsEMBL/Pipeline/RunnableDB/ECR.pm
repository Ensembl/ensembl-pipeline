# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::ECR

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $genscan = Bio::EnsEMBL::Pipeline::RunnableDB::ECR->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->write_output(); #writes to DB


=head1 DESCRIPTION

This is a RunnableDB wrapper for Bio::EnsEMBL::Pipeline::Runnable::ECR,
allowing query sequences to fetched from, and results to be written to,
the given database handle. 


=cut
package Bio::EnsEMBL::Pipeline::RunnableDB::ECR;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::ECR;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::ECR qw(ECR_FILTERS);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for exonerate from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my( $self) = @_; 
  
  $self->throw("No input id") unless defined($self->input_id);

  my $input_id  = $self->input_id;
  
  my ($contig, $chr_name, $chr_start, $chr_end);
 
  if ($input_id =~ /$SLICE_INPUT_ID_REGEX/) {
    ($chr_name, $chr_start, $chr_end) = ($1, $2, $3);
    $contig  = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr_name,$chr_start,$chr_end);
  } else {
    $contig = $self->db->get_RawContigAdaptor->fetch_by_name($input_id);
  }

  $self->query($contig);
  
  my (@feat_lists, $merge, $union);

  if (exists($ECR_FILTERS->{$self->analysis->logic_name}) and
      exists($ECR_FILTERS->{$self->analysis->logic_name}->{sources})) {

    if (exists($ECR_FILTERS->{$self->analysis->logic_name}->{union})) {
      $union = $ECR_FILTERS->{$self->analysis->logic_name}->{union};
    } else {
      $union = 0;
    }
    if (exists($ECR_FILTERS->{$self->analysis->logic_name}->{merge})) {
      $merge = $ECR_FILTERS->{$self->analysis->logic_name}->{merge};
    } else {
      $merge = 0;
    }

    foreach my $source (@{$ECR_FILTERS->{$self->analysis->logic_name}->{sources}}) {      
       push @feat_lists, $contig->get_all_DnaAlignFeatures($source);
    }

  } else {
    $self->throw("Dont know how to retrieve features for " . 
                 $self->analysis->logic_name .
                 " ; no entry in config with 'sources' filled in\n");
  }

  my $ecr = Bio::EnsEMBL::Pipeline::Runnable::ECR->new(-features => \@feat_lists,
                                                       -union    => $union,
                                                       -merge    => $merge);

  $self->runnable($ecr);
}


=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the Blastz analysis, producing Bio::EnsEMBL::DnaDnaAlignFeature
    Returns :   
    Args    :   

=cut

sub run {
  my ($self) = @_;
  
  $self->throw("Can't run - no runnable objects") unless defined($self->runnable);
  
  foreach my $run ($self->runnable) {
    $run->run;
  }

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
  
  my $fa = $db->get_SimpleFeatureAdaptor;
  
  foreach my $output (@features) {
    $output->attach_seq($self->query);    
    $output->seqname($self->query->id);
    $output->analysis($self->analysis);
    $output->display_label($self->analysis->logic_name);

    if ($self->query->isa("Bio::EnsEMBL::Slice")) {
      $output->is_splittable(1);
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
