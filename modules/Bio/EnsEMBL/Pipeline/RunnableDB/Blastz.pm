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

  # if repeatmasking required, always soft mask for blastz
  my $genseq = (@$PIPELINE_REPEAT_MASKING) 
      ? $contig->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING, 1)
      : $contig;

  my @db;

  my $executable =  $self->analysis->program_file;
  $executable = "$BIN_DIR/blastz" if not $executable;

  my $database = $self->analysis->db;
  $self->throw("RunnableDB/Blastz error: you must define a database in your analysis") if not $database;
  
  if ( -d $database) {
    @db = glob("$database/*");
    
    # do largest files first, so that memory problems, if they occur, occur early
    @db = sort { -s $b <=> -s $a } @db;
    
  } else {
    push(@db,$database);
  }

  foreach my $db (@db) {
    my $blastz = new Bio::EnsEMBL::Pipeline::Runnable::Blastz(-query         => $genseq,
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
