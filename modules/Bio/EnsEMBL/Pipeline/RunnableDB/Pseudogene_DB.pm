# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Pseudogene_DB.pm

=head1 SYNOPSIS



=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Pseudogene.pm to add
functionality to read from databases (so far).

=head1 CONTACT

Describe contact details here

=head1 APPENDIX



=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Pseudogene_DB;

use strict;
use Carp qw(cluck);
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Pseudogene;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pseudogene_config;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Pseudogene.pm from the database
    Returns :   none
    Args    :   none

=cut



sub fetch_input {
  my( $self) = @_;
    
  $self->throw("No input id") unless defined($self->input_id);

  $self->fetch_sequence;
  my $results = [];		# array ref to store the output
  my %parameters = $self->parameter_hash;
  my %repeat_blocks;
  $parameters{'-query'} = $self->query;



  my $runname = "Bio::EnsEMBL::Pipeline::Runnable::Pseudogene";

  my $rep_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_DBHOST,
     '-user'   => $GB_DBUSER,
     '-dbname' => $GB_DBNAME,
     '-pass'   => $GB_DBPASS,
     '-port'   => $GB_DBPORT,
     '-dnadb'  => $self->db,
    );
  my $genes = $self->query->get_all_Genes;
  foreach my $gene (@{$genes}) {
    my $rsa = $rep_db->get_SliceAdaptor;
    my $rep_gene_slice = $rsa->fetch_by_region(
					       'toplevel',
					       $self->query->chr_name,
					       $gene->start,
					       $gene->end,
					      );

    $repeat_blocks{$gene} = $self->get_all_repeat_blocks($rep_gene_slice->get_all_RepeatFeatures);

  }

  if ($self->validate_genes($genes)) {
    my $runnable = $runname->new
      ( 
       '-max_intron_length' => $PS_MAX_INTRON_LENGTH,
       '-max_intron_coverage' => $PS_MAX_INTRON_COVERAGE,
       '-max_exon_coverage' => $PS_MAX_EXON_COVERAGE,
       '-genes' => $genes,
       '-repeat_features' => \%repeat_blocks
      );
    $self->runnable($runnable);
    $self->results($runnable->output);
  }
  return 1;
}


=head2 get_all_repeat_blocks

  Args       : none
  Description: merges repeats into blocks for each gene
  Returntype : array of Seq_Feature blocks;

=cut 

sub get_all_repeat_blocks {
  my ($self,$repeat_ref) = @_;
  my @repeat_blocks;
  my @repeats = @{$repeat_ref};
  @repeats = sort {$a->start <=> $b->start} @repeats;
  my $curblock = undef;

  REPLOOP: foreach my $repeat (@repeats) {
      my $rc = $repeat->repeat_consensus;
      if ($rc->repeat_class !~ /LINE/ && $rc->repeat_class !~ /LTR/ && $rc->repeat_class !~ /SINE/) { 
	next REPLOOP;
      }
      if ($repeat->start <= 0) { 
	$repeat->start(1); 
      }
      if (defined($curblock) && $curblock->end >= $repeat->start) {
	if ($repeat->end > $curblock->end) { 
	  $curblock->end($repeat->end); 
	}
      }
      else {
	$curblock = Bio::EnsEMBL::SeqFeature->new(
						  -START => $repeat->start,
						  -END => $repeat->end, 
						  -STRAND => $repeat->strand
						 );
	push (@repeat_blocks,$curblock);
      }
    }
  @repeat_blocks = sort {$a->start <=> $b->start} @repeat_blocks;
  return\@repeat_blocks;
}



sub write_output {
  my($self) = @_;
  my @genes = @{$self->output};
  # write genes out to a different database from the one we read genes from.
  my $dbname = $GB_FINALDBNAME;
  my $dbhost = $GB_FINALDBHOST;
  my $dbuser = $GB_FINALDBUSER;
  my $dbpass = $GB_FINALDBPASS;
  my $dbport = $GB_FINALDBPORT;
  
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					      '-host'   => $dbhost,
					      '-user'   => $dbuser,
					      '-dbname' => $dbname,
					      '-pass'   => $dbpass,
					      '-port'   => $dbport,
					      '-dnadb'  => $self->db,
					     );
  # sort out analysis
  my $analysis = $self->analysis;
  unless ($analysis){
    $self->throw("an analysis logic name must be defined in the command line");
  }
  
  my $gene_adaptor = $db->get_GeneAdaptor;
  foreach my $gene (@genes) { 
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->translation;
    }
    $gene->analysis($analysis);
    # store
    eval {
     $gene_adaptor->store($gene);
      print STDERR  "wrote gene " . $gene->dbID . " to database ".
	$gene_adaptor->db->dbname."\n";
    };
    if ( $@ ) {
      $self->warn("UNABLE TO WRITE GENE:\n$@");
    }
  }
}

sub validate_genes {
  my ($self,$genes) =@_;
   foreach my $gene (@{$genes}) {
     unless ($gene->isa('Bio::EnsEMBL::Gene')){
       $self->throw('Object is not a valid gene, it is a $gene');
     }
   }
  return 1;
}

sub run  {

    my ($self) = @_;

    foreach my $runnable ($self->runnable) {

      $self->throw("Runnable module not set") unless ($runnable);

      # Not sure about this
      $self->throw("Input not fetched")       unless ($self->query);

      $runnable->run();
      if($self->validate_genes){
	$self->output($runnable->output);
      }
    }
    return 1;
}

sub output{
    my ($self,$genes ) = @_;
    unless ( $self->{_genes} ){
	$self->{_genes} = [];
    }
    if ($genes){
	$self->{_genes} = $genes;
    }
    return $self->{_genes};
}

1;
