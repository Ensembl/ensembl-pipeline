# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Blastz

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $genscan = Bio::EnsEMBL::Pipeline::RunnableDB::Blastz->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->write_output(); #writes to DB


=head1 DESCRIPTION

This is a RunnableDB wrapper for Bio::EnsEMBL::Pipeline::Runnable::Blastz,
allowing query sequences to fetched from, and results to be written to,
the given database handle. 


=cut
package Bio::EnsEMBL::Pipeline::RunnableDB::Blastz;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Blastz;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


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
  
  &throw("No input id") unless defined($self->input_id);

  my $sa = $self->db->get_SliceAdaptor;
  my $slice = $sa->fetch_by_name($self->input_id);
  if(@$PIPELINE_REPEAT_MASKING){
    # blastz works best with soft masking
    my $sequence = $slice->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING, 1);
    $self->query($sequence);
  }else{
    $self->query($slice);
  }

  my $executable =  $self->analysis->program_file;
  $executable = "$BIN_DIR/blastz" if not $executable;

  my $database = $self->analysis->db;
  &throw("You must define a database to search against in the analysis") if not $database;
  
  my @db;
  if ( -d $database) {
    @db = glob("$database/*");
    
    # if the files have standard names, try to sort them 
    # for consistency
    @db = sort { my ($o) = ($a =~ /\/([^\.\/]+)[^\/]*\.fa$/); 
                 my ($t) = ($b =~ /\/([^\.\/]+)[^\/]*\.fa$/); 
                 $o <=> $t } @db;
    
  }
  else {
    push(@db,$database);
  }

  foreach my $db (@db) {
    my $blastz = new Bio::EnsEMBL::Pipeline::Runnable::Blastz(-query         => $self->query,
                                                              -program       => $executable,
                                                              -database      => $db,
                                                              -options       => $self->analysis->parameters);    
    $self->runnable($blastz);
  }
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
  
  &throw("Can't run - no runnable objects") unless defined($self->runnable);
  
  foreach my $run ($self->runnable) {
    $run->run;
  }

}


1;
