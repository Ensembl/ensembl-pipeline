# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Pseudogene_DB.pm


=head1 SYNOPSIS

my $runnabledb = Bio::EnsEMBL::Pipeline::RunnableDB::Pseudogene_DB->new(
									-db => $db_adaptor,
									-input_id => $slice_id,		
									-analysis => $analysis,
								       );

$runnabledb->fetch_input();
$runnabledb->run();
my @array = @{$runnabledb->output};
$runnabledb->write_output();

array ref returned by output contain all the genes found on the slice, modified to relflect pseudogene status

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Pseudogene.pm 

Opens connections to 3 dbs:
1 for repeat sequences (GB_DB)
1 for fetching genes from (GB_FINAL)
1 db connection for querying align_feature tables for spliced elsewhere tests

fetches all the genes on the slice, all the repeats associtaed with each gene and 
collects alignment feature evidence for single exon genes and passes them to the 
runnable.


=head1 CONTACT

Describe contact details here

=head1 APPENDIX



=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Pseudogene_DB;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Pseudogene;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBConnection;
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
  my %homolog_hash;
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
  #store repeat db internally
  $self->rep_db($rep_db);

  my $genes_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_FINALDBHOST,
     '-user'   => $GB_FINALDBUSER,
     '-dbname' => $GB_FINALDBNAME,
     '-pass'   => $GB_FINALDBPASS,
     '-port'   => $GB_FINALDBPORT,
    );

  #genes come from final genebuild database
  my $gene_db_connection = new Bio::EnsEMBL::DBSQL::DBConnection
    (
     '-host'   => $GB_FINALDBHOST,
     '-user'   => $GB_FINALDBUSER,
     '-dbname' => $GB_FINALDBNAME,
     '-pass'   => $GB_FINALDBPASS,
     '-port'   => $GB_FINALDBPORT,
    );


  #store gene db internally
  $self->gene_db($genes_db);
  $self->gene_db_connection($gene_db_connection); 

  my $genedb_sa = $genes_db->get_SliceAdaptor;
  my $genes_slice = $genedb_sa->fetch_by_region(
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

    my $chromosome_slice = $rsa->fetch_by_region(
						 'chromosome',
						 $self->query->chr_name,
						);
    
    my $transferred_gene = $gene->transfer($chromosome_slice);

    # for testing purposes
###########################################################
##########################################################
    $transferred_gene->type('unset');                   #
##########################################################
###########################################################
    push @transferred_genes,$transferred_gene;

    if (scalar(@{$transferred_gene->get_all_Exons()})==1){
      # single exon - check for retrotransposition
      $homolog_hash{$transferred_gene} = $self->spliced_elsewhere($transferred_gene,$genedb_sa);
    }
    else{
      # multiexon - check for repeats in introns
      $repeat_blocks{$transferred_gene} = $self->get_all_repeat_blocks($rep_gene_slice->get_all_RepeatFeatures);
    }
  }
  if ($self->validate_genes(\@transferred_genes)) {

    my $runnable = $runname->new
      ( 
       '-genes' => \@transferred_genes,
       '-repeat_features' => \%repeat_blocks,
       '-homologs'        => \%homolog_hash,
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
    } else {
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

  #  empty_Analysis_cache();

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

  my $gene_adaptor = $db->get_GeneAdaptor;
  foreach my $gene (@genes) { 
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->translation;
    }

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
    if ($self->validate_genes) {
      $self->output($runnable->output);
    }
  }
  return 1;
}

sub output{
  my ($self,$genes ) = @_;
  unless ( $self->{_genes} ) {
    $self->{_genes} = [];
  }
  if ($genes) {
    $self->{_genes} = $genes;
  }
  return $self->{_genes};
}

sub spliced_elsewhere{
  my ($self,$gene,$sa) = @_;
  my $table;
  my %identified_genes;
  my $connection = $self->gene_db_connection;
  my $exon = @{$gene->get_all_Exons()}[0];
  my %homolog;
  

#####################################################
# Fetch supporting features for the single exon

  foreach my $sf (@{$exon->get_all_supporting_features}){

    if ($sf->isa('Bio::EnsEMBL::DnaDnaAlignFeature')){
	$table = "dna_align_feature";
      }
    if ($sf->isa('Bio::EnsEMBL::DnaPepAlignFeature')){
      $table = "protein_align_feature";
    }
    my $sql = "select seq_region.name,$table.seq_region_start,$table.seq_region_end,$table.hit_name,$table.score
from $table,seq_region
where $table.hit_name ='".$sf->hseqname."' 
and $table.score > $PS_FEATURE_SCORE
and $table.seq_region_id = seq_region.seq_region_id 
limit $PS_SQL_LIMIT;";
    
    my $sth = $connection->prepare($sql);
    $sth->execute;
    while (my @row = $sth->fetchrow_array()) {
      my $slice = $sa->fetch_by_region('seqlevel',
				       $row[0],
				       $row[1],
				       $row[2],
				      );

      ############################################################
      # For each hit in the protein / dna align feature table
      # get a slice for each locus that the feature aligns to 
      # return all the genes in that slice and store them in a hash

      my @homologs = @{$slice->get_all_Genes};
      foreach my $homologous_gene(@homologs){ 
	# only look at each homologous gene once
	# print $homologous_gene->display_id."\n";
	next unless ($identified_genes{$homologous_gene->dbID} == 0);
	$identified_genes{$homologous_gene->dbID} = 10;
	if ($gene->stable_id ne $homologous_gene->stable_id){
	  # dont look at self alignments
	  foreach my $homologous_transcript(@{$homologous_gene->get_all_Transcripts}){
	    # homologs trascript must have more than 1 exon
	    next unless scalar(@{$homologous_transcript->get_all_Exons()}) > 1;
	    $homolog{$homologous_transcript} = $homologous_transcript;
	  }
	}
      }
    }
    $sth->finish();
  }

return \%homolog;
}


############################################################################
# container methods

=head2 gene_db

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set gene db adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub gene_db {
  my ($self, $gene_db) = @_;
  
  unless ($gene_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
  $self->throw("gene db is not a Bio::EnsEMBL::DBSQL::DBAdaptor, it is a $gene_db");
}
    $self->{'_gene_db'} = $gene_db;
  return $self->{'_gene_db'};
}

=head2 gene_db_connection

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection
  Description: get/set gene db conection used for sql query in spliced_elsewhere
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : general

=cut

sub gene_db_connection {
  my ($self, $gene_db_connection) = @_;
  if ($gene_db_connection){
  unless ($gene_db_connection->isa("Bio::EnsEMBL::DBSQL::DBConnection")){
    $self->throw("gene db is not a Bio::EnsEMBL::DBSQL::DBConnection, it is a $gene_db_connection");
    }
    $self->{'_gene_db_connection'} = $gene_db_connection;
  }
  return $self->{'_gene_db_connection'};
}

=head2 rep_db

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
  Description: get/set gene db adaptor
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub rep_db {
  my ($self, $rep_db) = @_;
  
  unless ($rep_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
  $self->throw("gene db is not a Bio::EnsEMBL::DBSQL::DBAdaptor, it is a $rep_db");
}
  $self->{'_rep_db'} = $rep_db;
  return $self->{'_rep_db'};
}




1;
