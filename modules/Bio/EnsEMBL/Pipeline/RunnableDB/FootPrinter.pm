bg#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB:FootPrinter

=head1 SYNOPSIS

my $fp = Bio::EnsEMBL::Pipeline::RunnableDB::FootPrinter->new(
                            -db              => $db,
                            -analysis        => $analysis,
                            -multispecies_db => $msdb,
                            -gene_stable_id  => $gene_sid);

$fp->fetch_input;
$fp->run;
my @output = $fp->output;
$fp->write_output;

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::FootPrinter to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::FootPrinter;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::FootPrinter;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

#    Title   :   new
#    Usage   :   $self->new(-db              => $db,
#                           -multispecies_db => $msdb,
#                           -gene_stable_id  => $gene_sid);
#    Function:   
#    Returns :   
#    Args    :   

=cut

sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;

  # Explanation of input variables:
  #
  # write_db        - ensembl core db where footprinter features will be stored.
  # multispecies_db - database of homologous upstream sequences.
  # compara_db      - the compara database used to build the multispecies db.
  # species_list    - (optional) list of species that footprinter should look 
  #                   for conserved promoters et al.
  # analysis        - a footprinter analysis object.
  # gene_stable_id  - id of gene for which one wants to find putative promoters.
  # species         - the organism from which the above gene originates (binomen)
  # score_cutoff    - conserved upstream bases must have a score greater than this.
  #                   This value is the score calculated during the generation of
  #                   the conserved regions stored in the multispecies db.
  # perc_id_cutoff  - The identity cutoff that blastz matches must exceed in order
  #                   to be counted towards scoring conserved upstream regions.
  # upstream_length - Length of region upstream of gene in which to look for 
  #                   promoters.  Note that if this is longer than that for which
  #                   the multispecies database has been calculated then conserved
  #                   regions will not be found in this 'gap'.
  # tree            - Footprinter inherantly needs a tree.  These should be passed
  #                   in as a newick/new hampshire string.  Only use binomens to 
  #                   refer to species. Default tree is: (((('Rattus norvegicus', 
  #                   'Mus musculus'),'Homo sapiens'),('Danio rerio', 'Fugu rubripes')),
  #                   ('Anopheles gambiae', ('Drosophila melanogaster', 
  #                   'Drosophila pseudoobscura')));


    my ($write_db, 
	$multispecies_db, 
	$compara_db,
	$species_list,
	$tree,
	$analysis, 
	$gene_stable_id,
	$species,
	$score_cutoff
	$perc_id_cutoff
	$upstream_length) = 
	  $self->_rearrange ([WRITE_DB
			      MULTISPECIES_DB
			      COMPARA_DB
			      SPECIES_LIST
			      TREE
			      ANALYSIS
			      GENE_STABLE_ID
			      SPECIES
			      SCORE_CUTOFF
			      PERC_ID_CUTOFF
			      UPSTREAM_LENGTH], @args);

  $self->_write_db($write_db)               if $write_db;
  $self->_multispecies_db($multispecies_db) if $multispecies_db;
  $self->_compara_db($compara_db)           if $compara_db;
  $self->_species_list($species_list)       if $species_list;
  $self->_tree($tree)                       if $tree;
  $self->analysis($analysis)                if $analysis;
  $self->_gene_stable_id($gene_stable_id)   if $gene_stable_id;
  $self->_species($species)                 if $species;
  $self->_score_cutoff($score_cutoff)       if $score_cutoff;
  $self->_perc_id_cutoff($perc_id_cutoff)   if $perc_id_cutoff;
  $self->_upstream_length($upstream_length) if $upstream_length;

  return $self
}



=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;    

    # Fetching the seqs we want is a little tricky.  We use a swag of
    # adaptors to connect to the compara database, the multispecies
    # match database and the core database for each species of interest.
    # Using these adaptors the 
    # _multispecies_db->retrieve_homologous_upstream_sequences method 
    # retrieves a set of orthologs for our gene (defined by 
    # primary_gene_stable_id) and fetches the upstream conserved region 
    # from each of these orthologous sequences.  These sequences are
    # then passed to footprinter.
    
    my $conserved_upstream_seqs = 
      $self->_multispecies_db->retrieve_homologous_upstream_sequences(
                                                $self->_species_list,
                                                $self->_core_db_list,
                                                $self->_compara_db,
                                                $self->_gene_stable_id,
                                                $self->_species,
                                                $self->_score_cutoff, 
                                                $self->_perc_id_cutoff, 
                                                $self->_upstream_length);

    # There is a small glitch here.  We have to pass a tree to footprinter
    # so that it will work properly.  Often, more than one ortholog will 
    # turn up for a given species.  This is not a problem if they are the
    # result of a recent duplication, but sometimes these orthologs are
    # 'lost'.  To cope with this and gratuitously make computation
    # easier, we will assume that the sequence with the most conserved
    # bases is the most similar ortholog and discard the rest.

      # Partition by species

    my %seqs_by_species;

    foreach my $ortholog_id (keys %$conserved_upstream_seq){

      $conserved_upstream_seq{$ortholog_id}->{_gene_stable_id} = $ortholog_id;

      push (@{$seqs_by_species{$conserved_upstream_seq{$ortholog_id}->{species}}}, 
	$conserved_upstream_seq{$ortholog_id})
    }

      # Filter for longest conserved region per species.

    my @final_sequence_set;

    foreach my $species (keys %$seqs_by_species) {

      my $longest = 0;
      my $long_seq;

      foreach my $conserved_seq (@{$seqs_by_species{$species}}){

	my $seq = $conserved_seq{_seq};
	$seq =~ s/N//g;

	$long_seq = $conserved_seq 
	  if (length($seq) > $longest);
      }

      # Make a seq that can be passed to footprinter

      my $bio_seq = Bio::Seq->new(-display_id => $long_seq{_species},
				  -seq        => $long_seq{_seq});
      
      push (@final_sequence_set, $bio_seq)
    }


    my $primary_coords = $seqs_by_species{$self->_species}->{_coords};

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::FootPrinter->new(
			    conserved_seqs => \@final_sequence_set,
			    tree           => $self->_tree,
			    parameters     => $self->_parameters,
			    primary_id     => $self->_species,
			    primary_coords => $primary_coords,
			    primary_strand => $self->_primary_strand,
			    db             => $self->_write_db,    ### Purely for getting the analysis.
			    executable     => $self->_executable   ### !!!  Where ?!?!
		    );

							   
    $self->runnable($runnable);
    
    return 1;
}

### Storage ###

sub _write_db {
  my $self =  shift;

  if (@_) {
    $self->{_write_db} = shift;
  }

  $self->throw("Core write db not specified!") 
    unless $self->{_write_db};

  return $self->{_write_db}
}


sub _multispecies_db {
  my $self =  shift;

  if (@_) {
    $self->{_multispecies_db} = shift;
  }

  $self->throw("Need a multispecies db object!") 
    unless $self->{_multispecies_db};

  return $self->{_multispecies_db}
}


sub _compara_db {
  my $self =  shift;

  if (@_) {
    $self->{_compara_db} = shift;
  }

  $self->throw("Need a compara db object!") 
    unless $self->{_compara_db};

  return $self->{_compara_db}
}


sub _gene_stable_id {
  my $self =  shift;

  if (@_) {
    $self->{_gene_stable_id} = shift;
  }

  $self->throw("Need a gene stable id!") 
    unless $self->{_gene_stable_id};

  return $self->{_gene_stable_id}
}


sub _species {
  my $self =  shift;

  if (@_) {
    $self->{_species} = shift;

    # Always lc

    $self->{_species} =~ s/[A-Z]/[a-z]/g;

  }

  $self->throw("Need an input species to go with gene!") 
    unless $self->{_species};

  return $self->{_species}
}

sub _tree {
  my $self =  shift;

  if (@_) {
    $self->{_tree} = shift;

    # Always lc
    $self->{_tree} =~ tr/[A-Z]/[a-z]/;
  }
    
  # Set default tree

  unless ($self->{_tree}) {
    $self->{_tree} = "((((\'rattus norvegicus\', \'mus musculus\'),\'homo sapiens\'),".
      "(\'danio rerio\', \'fugu rubripes\')),(\'anopheles gambiae\', ".
	"(\'drosophila melanogaster\', \'drosophila pseudoobscura\')));";
  }

  # Sanity check - are all species in our
  # species list represented in our tree
  # and vice versa.
  
  my $tree = $self->{_tree};
  $tree =~ s/[\(\)\;\:\d\'\"]; # Should only be text and commas left.
  my @species = split /\,/, $tree;
  
  foreach my $list_species (@{$self->_species_list}) {
    
    my $match = 0;
    
    foreach my $tree_species (@species){
      
      if ($tree_species eq $list_species) {

	$match = 1;
	last;
      }
    }
    $self->throw("Species  [$list_species] is represented in tree [".$self->_tree."].") 
      unless $match;
  }
  
  return $self->{_tree}
}


sub _species_list {
  my $self =  shift;

  if (@_) {
    $self->{_species_list) = shift;

    # Always lc

    my @new_array;

    foreach my $species (@{$self->{_species_list}}) {

      $species =~ s/[A-Z]/[a-z]/g;
      push (@new_array, $species)
    }

    $self->{_species_list} = \@new_array;
  }

  # If the species list is undefined, return all

  unless ($self->{_species_list}) {
    my $genome_adaptor = $self->_compara_db->get_GenomeDBAdaptor;
    my $genomes = $genome_adaptor->fetch_all;

    foreach my $genome (@$genomes) {

      # Always lc
      my $genome_name = $genome->name;
      $genome_name =~ s/[A-Z]/[a-z]/g;
      push (@{$self->{_species_list}}, $genome_name)
    }
  }

  return $self->{_species_list}
}


sub _core_db_list {
  my $self = shift;

  if (@_) {
    $self->{_core_db_list) = shift;
  }

  # The core db list can be set manually, as above, but it is unlikely
  # to be, especially since this information can be derived from
  # the multispecies db.

  unless ($self->{_core_db_list}) {

    foreach my $species (@{$self->{_species_list}) {

      push (@{$self->{_core_db_list}}, 
	    $self->_multispecies_db->_get_core_adaptor($species))
    }
  }

  return $self->{_core_db_list}
}

1;
