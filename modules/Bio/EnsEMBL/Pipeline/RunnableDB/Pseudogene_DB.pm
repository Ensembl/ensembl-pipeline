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
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pseudogene_config;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace);

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

  my $results = [];		# array ref to store the output
  my %parameters = $self->parameter_hash;
  my %repeat_blocks;
  $parameters{'-query'} = $self->query;
  $self->fetch_sequence;
  my $runname = "Bio::EnsEMBL::Pipeline::Runnable::Pseudogene";
  
  my $rep_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_DBHOST,
     '-user'   => $GB_DBUSER,
     '-dbname' => $GB_DBNAME,
     '-pass'   => $GB_DBPASS,
     '-port'   => $GB_DBPORT,
    );

  my $genes_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_FINALDBHOST,
     '-user'   => $GB_FINALDBUSER,
     '-dbname' => $GB_FINALDBNAME,
     '-pass'   => $GB_FINALDBPASS,
     '-port'   => $GB_FINALDBPORT,
    );

#genes come from final genebuild database

  my $genes_slice = $genes_db->get_SliceAdaptor->fetch_by_region(
								 'chromosome',
								 $self->query->chr_name,
								 $self->query->start,
								 $self->query->end
								);
  my $genes = $genes_slice->get_all_Genes;
  my @transferred_genes;

  foreach my $gene (@{$genes}) {
    my $rsa = $rep_db->get_SliceAdaptor;

#repeats come from core database

    my $rep_gene_slice = $rsa->fetch_by_region(
					       'chromosome',
					       $self->query->chr_name,
					       $gene->start,
					       $gene->end,
					      );

##########################################################################
#transfer gene coordinates to entire chromosome to prevent problems arising 
# due to offset with repeat features 

    my $chromosome_slice = $self->query->adaptor->fetch_by_region(
									'chromosome',
									$self->query->chr_name,
								       );

    my $transferred_gene = $gene->transfer($chromosome_slice);
    push @transferred_genes,$transferred_gene;

    $repeat_blocks{$transferred_gene} = $self->get_all_repeat_blocks($rep_gene_slice->get_all_RepeatFeatures);
  }

  if ($self->validate_genes(\@transferred_genes)) {
    my $runnable = $runname->new
      ( 
       '-max_intron_length' => $PS_MAX_INTRON_LENGTH,
       '-max_intron_coverage' => $PS_MAX_INTRON_COVERAGE,
       '-max_exon_coverage' => $PS_MAX_EXON_COVERAGE,
       '-genes' => \@transferred_genes,
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

=head2 write_output

  Args       : none
  Description: writes gene array into db specified in Bio::EnsEMBL::Config::GeneBuild::Databases.pm
  exception  : warns if it cannot write gene
  Returntype : array of Seq_Feature blocks;

=cut 


sub write_output {
  my($self) = @_;
  my @genes = @{$self->output};

  empty_Analysis_cache();

  # write genes out to a different database from the one we read genes from.
  my $dbname = $PSEUDO_DBNAME;
  my $dbhost = $PSEUDO_DBHOST;
  my $dbuser = $PSEUDO_DBUSER;
  my $dbpass = $PSEUDO_DBPASS;
  my $dbport = $PSEUDO_DBPORT;
  
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					      '-host'   => $dbhost,
					      '-user'   => $dbuser,
					      '-dbname' => $dbname,
					      '-pass'   => $dbpass,
					      '-port'   => $dbport,
					     );
  # sort out analysis
  my $analysis = $self->analysis;
  print STDERR $analysis->logic_name."\n";
  unless ($analysis){
    $self->throw("an analysis logic name must be defined in the command line");
  }

  #make sure all the anaysis analysis objects for objects associated with the gene
  #have adaptors that point to the new database 


  my $gene_adaptor = $db->get_GeneAdaptor;
  foreach my $gene (@genes) { 
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->translation;
      foreach my $exon (@{$trans->get_all_Exons}) {
        foreach my $sf (@{$exon->get_all_supporting_features}) {
          $sf->analysis(get_Analysis_by_type($gene_adaptor->db,
                                             $sf->analysis->logic_name));
        }
      }
    }
    my $anal = get_Analysis_by_type($gene_adaptor->db, 
                                    $analysis->logic_name);

    $gene->analysis($anal);

    eval {
      #make sure display_xrefs point in the right direction
      if (defined ($gene->display_xref->adaptor)){
	$gene->display_xref->adaptor($gene_adaptor);
      }      
      #make sure xrefs point in the right direction
      foreach my $xrefs (@{$gene->get_all_DBLinks}){
	if (defined ($xrefs->analysis)){
	  $xrefs->analysis(get_Analysis_by_type($gene_adaptor->db,
						$xrefs->analysis->logic_name));;
	}
      }
    };
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

{
  my (%ana);

  sub empty_Analysis_cache {
    %ana = ();
  }

  sub get_Analysis_by_type {
    my ($dba, $type) = @_;
    unless (exists($ana{$type})) {
      my $logic    = $type;
      my $ana_aptr = $dba->get_AnalysisAdaptor;
      $ana{$type}  = $ana_aptr->fetch_by_logic_name($logic);
      if (!defined($ana{$type})) {
        $ana{$type} = new Bio::EnsEMBL::Analysis( -logic_name => $type);
        $ana_aptr->store($ana{$type});
      }
    }
    return $ana{$type};
  }
}


=head2 validate_gene
  Args       : array ref to Bio::EnsEMBL::Gene objects
  Description: checks Gene objects
  Exceptions : throws if not Bio::EnsEMBL::Gene
  Returntype : scalar

=cut 


sub validate_genes {
  my ($self,$genes) =@_;
   foreach my $gene (@{$genes}) {
     unless ($gene->isa('Bio::EnsEMBL::Gene')){
       $self->throw('Object is not a valid gene, it is a $gene');
     }
   }
  return 1;
}


=head2 run
  Args       : none
  Description: overrides runnableDb run method to allow gene objects to be validated 
before runnning the runnable
  Returntype : scalar

=cut 

sub run  {

    my ($self) = @_;

    foreach my $runnable ($self->runnable) {

      $self->throw("Runnable module not set") unless ($runnable);

      # Not sure about this
 #     $self->throw("Input not fetched")       unless ($self->query);

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
