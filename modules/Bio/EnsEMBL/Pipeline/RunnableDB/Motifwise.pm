#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::Motifwise

=head1 SYNOPSIS

  my $mw = Bio::EnsEMBL::Pipeline::RunnableDB::Motifwise->new(
				    -db              => $db,
				    -input_id        => $input_id,
                                    -analysis_tss    => $analysis_tss, 
                                    -analysis_motif  => $analysis_motif,
                                    );

  $mw->fetch_input;
  $mw->run;
  my @output = $mw->output;
  $mw->write_output;

  The rows needed in the analysis table on the core database are:

INSERT INTO analysis (logic_name, db, db_file, program, program_version, program_file, parameters, module, gff_source, gff_feature) VALUES ('motifwise_tss','motif','/ecs2/work1/dta/motifwise/motif_dropout.lr','motifwise','beta','motifwise','-tfm_cutoff 11.0','Motifwise','motifwise','TSS');
INSERT INTO analysis (logic_name, db, db_file, program, program_version, program_file, parameters, module, gff_source, gff_feature) VALUES ('motifwise_motif','motif','/ecs2/work1/dta/motifwise/motif_dropout.lr','motifwise','beta','motifwise','-tfm_cutoff 11.0','Motifwise','motifwise','Motif');

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Motifwise to add
functionality for reading and writing to databases. The two
Bio::EnsEMBL::Analysis objects must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Cared for by Dan Andrews (dta@sanger.ac.uk)

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Motifwise;

use strict;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Motifwise;
use Bio::EnsEMBL::Pipeline::Config::General qw(BIN_DIR
					       SLICE_INPUT_ID_REGEX
					       PIPELINE_REPEAT_MASKING
					      );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new
  
Title   :   new
Usage   :   $self->new(_things_);
Function:   Creates a new RunnableDB::Motifwise object
Returns :   Bio::EnsEMBL::Pipeline::RunnableDB::Motifwise;
Args    :   
  
=cut
  
sub new {
  my ($class, @args) = @_;
    
  # This is a non-standard constructor because we use two analysis
  # objects to store the motifwise output (one for 'regions' or 
  # transcription start sites and another for the actual motif
  # locations).  The two analysis objects have logic names like
  # motifwise_tss and motifwise_motif.
    
  # Get our analysis objects.

  my ($analysis_tss,
      $analysis_motif,
      $workdir) = Bio::EnsEMBL::Root->_rearrange([qw(ANALYSIS_TSS
						      ANALYSIS_MOTIF
						      WORK_DIR)],
						  @args);
  
  Bio::EnsEMBL::Root->throw("Missing either one or more analysis objects.  Need an analysis\n" . 
			    "object for TSS features and another analysis object for Motif\n" .
			    "features")
    unless (($analysis_tss && $analysis_tss->isa("Bio::EnsEMBL::Analysis"))
	    &&($analysis_motif && $analysis_motif->isa("Bio::EnsEMBL::Analysis")));
  
  # Add one of the analysis objects to our array of arguments
  # that are passed to the parent class constructor.

  push @args, ('-analysis', $analysis_tss);


  # Make object using parent constructor.
  
  my $self = $class->SUPER::new(@args);
  
  # Manually add the analysis that wasn't passed to the parental
  # constructor.

  $self->_analysis_motif($analysis_motif);
  
  return $self; 
}


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;    
    
    # Check analysis exists
    $self->throw("Analysis object not specified") 
      unless $self->analysis->isa("Bio::EnsEMBL::Analysis");

    # Check input id
    $self->throw("No input id specified (eg. 1.500000-1000000)") 
      unless $self->input_id;
    
    # Fetch slice specified by input id
    unless ($self->input_id =~ /$SLICE_INPUT_ID_REGEX/ ){
      $self->throw("Input id [".$self->input_id."] not compatible with\n".
		   "Bio::EnsEMBL::Pipeline::Config::General::SLICE_INPUT_ID_REGEX\n".
		   "[$SLICE_INPUT_ID_REGEX]");
    }
    my $chr   = $1;
    my $start = $2;
    my $end   = $3;

    my $slice = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);
    $self->query($slice);

    my $seq;

    if (scalar @$PIPELINE_REPEAT_MASKING) {
      $seq = $slice->get_repeatmasked_seq($PIPELINE_REPEAT_MASKING);
    } else {
      $seq = $slice;
    }


    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Motifwise->new(
		     -query_seq  => $seq,
		     -executable => $BIN_DIR . "/motifwise",
                     -parameters => $self->parameters,
		     -motif_file => $self->_motif_file,
		     );

    $self->runnable($runnable);
    
    return 1;
}


sub write_output {
  my($self) = @_;
  
  my $sfa = $self->db->get_SimpleFeatureAdaptor;
  
  my @mapped_features;
  
  my $slice = $self->query;
  
  foreach my $simple_feat ($self->output) {

    $simple_feat->analysis($self->_analysis_tss) 
      if $simple_feat->seqname eq 'motifwise_tss';

   $simple_feat->analysis($self->_analysis_motif)
      if $simple_feat->seqname eq 'motifwise_motif';

    
    $simple_feat->contig($slice);
    my @mapped = $simple_feat->transform;
    
    if (@mapped == 0) {
      $self->warn("Couldn't map $simple_feat - skipping");
      next;
    }
    if (@mapped == 1 && $mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")) {
      $self->warn("$simple_feat seems to be on a gap - something bad has happened ...");
      next;
    }
    
    push @mapped_features, $mapped[0];    
  }

  $sfa->store(@mapped_features) if @mapped_features;
  
  return 1;
}

### Storage and retrieval

sub _motif_file {
  my $self = shift;

  return $self->analysis->db_file;
}

sub analysis {
  my ($self, @args) = @_;

  $self->_analysis_tss(@args);
}


sub _analysis_tss {
  my ($self, $analysis) = @_;

  if ($analysis) {
    $self->throw("Specified object is not a Bio::EnsEMBL::Analysis object")
      unless ($analysis && 
	      $analysis->isa("Bio::EnsEMBL::Analysis"));
    $self->{_analysis_tss} = $analysis;
    $self->parameters($analysis->parameters);
  }

  return $self->{_analysis_tss}
}


sub _analysis_motif {
  my $self = shift;

  $self->{_analysis_motif} = shift if @_;

  $self->throw("Specified object is not a Bio::EnsEMBL::Analysis object")
    unless ($self->{_analysis_motif} && 
	    $self->{_analysis_motif}->isa("Bio::EnsEMBL::Analysis"));

  return $self->{_analysis_motif}
}

sub _workdir {
  my $self = shift;

  $self->{_workdir} = shift if @_;

  return $self->{_workdir}
}

return 1;

