# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::AxtFilter

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

This is a RunnableDB wrapper for Bio::EnsEMBL::Pipeline::Runnable::AxtFilter


=cut
package Bio::EnsEMBL::Pipeline::RunnableDB::AxtFilter;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::AxtFilter;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::ECR qw(AXT_FILTERS SUBSETAXT_MATRIX SUBSETAXT_CUTOFF);

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
  my $logic_name = $self->analysis->logic_name;

  my ($contig, $chr_name, $chr_start, $chr_end);
 
  if ($input_id =~ /$SLICE_INPUT_ID_REGEX/) {
    ($chr_name, $chr_start, $chr_end) = ($1, $2, $3);
    $contig  = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr_name,$chr_start,$chr_end);
  } else {
    $contig = $self->db->get_RawContigAdaptor->fetch_by_name($input_id);
  }

  $self->query($contig);
  
  my %parameters = (-query => $contig);


  if (exists($AXT_FILTERS->{$logic_name}) and
      exists($AXT_FILTERS->{$logic_name}->{source})) {

    if ($AXT_FILTERS->{$logic_name}->{best}) {
      $parameters{-best} =1;
    }
    if ($AXT_FILTERS->{$logic_name}->{subset}) {
      $parameters{-subsetmatrix} = $SUBSETAXT_MATRIX;
      $parameters{-subsetcutoff} = $SUBSETAXT_CUTOFF;
    }
    if ($AXT_FILTERS->{$logic_name}->{target_lengths}) {
      $parameters{-targetlengths} = $AXT_FILTERS->{$logic_name}->{targetlengths};
    }

    $parameters{-targetnib} = $AXT_FILTERS->{$logic_name}->{targetnib};

    my @features = @{$contig->get_all_DnaAlignFeatures($AXT_FILTERS->{$logic_name}->{source})};
    $parameters{-features} = \@features;

    if (not @features) {
      $self->input_is_void(1);
      return;
    }

  } else {
    $self->throw("Dont know how to retrieve features for " . 
                 $logic_name .
                 " ; no entry in config with 'source' filled in\n");
  }

  my $run = Bio::EnsEMBL::Pipeline::Runnable::AxtFilter->new(%parameters);

  $self->runnable($run);
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
  
  my $fa = $db->get_DnaAlignFeatureAdaptor;
  
  foreach my $output (@features) {
    $output->attach_seq($self->query);    
    $output->analysis($self->analysis);

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
