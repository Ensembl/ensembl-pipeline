=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor

=head1 SYNOPSIS

$denorm_genes = Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor->new($db) ;

$denorm_genes->store($gene);
my $genes = $denorm_genes->fetch_all_by_Slice_and_type($slice);

=head1 DESCRIPTION

All features written to the exon and dna_align_feature tables are meant to be
represented by RawContig coordinates.  Hence, when running an analysis like BLAT
which is happiest with large slices or whole chromosomes a lot of coordinate
transformation needs to be carried out.  This is database intensive and hence
represents a bottleneck when BLAT jobs are run across many cpus.  When building
EST genes from BLAT results, all building is conducted on minigenomic sequence
via Slices - hence the RawContig coordinates stored with BLAT output are not used
except to map features to Slices.  All of the slow assembly mapping is a fairly 
BIG waste of time.

Using this module and a set of extra tables in the core database allows features
to be written in chromosomal coordinates.  When used sensibly this should allow
for a streamlined EST gene build.  This module basically allows genes to be stored
and retrieved in the coordinate system of the chromosomes or supercontigs fed 
into BLAT.


=head1 CONTACT

  Dan Andrews : dta@sanger.ac.uk

=head1 APPENDIX

=cut


package Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor;

use strict;

use Data::Dumper;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_INPUTID_REGEX
					EST_BLAT_ANALYSIS
				       );


use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 new

  Arg [1]    : list of arguments @args
  Example    : $denorm_gene_adaptor = Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor->new($db);
  Description: Creates a new DenormGeneAdaptor object
  Returntype : Bio::EnsEMBL::Pipeline::DBSQL::DenormGeneAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::Pipeline::RunnableDB::BlatToGenes, Bio::EnsEMBL::Pipeline::RunnableDB::EST_genebuilder

=cut

sub new {
  my($class, @args) = @_;

  $class = ref $class || $class;

  #call superclass constructor
  my $self = $class->SUPER::new(@args);

  return $self;
}

sub store {
  my ($self, $gene) = @_;

  # It is necessary to detach the current adaptors from our
  # gene objects before we store them.  These will need to be
  # re-attached/rejuvenated when they are revived.

  for (my $k = 0; 
       (($k <= scalar @{$gene->{'_transcript_array'}}) 
       && (defined $gene->{'_transcript_array'}->[$k])); 
       $k++){

    my $transcript = $gene->{'_transcript_array'}->[$k];

    for (my $j = 0; 
	 (($j <= scalar $transcript->{'_trans_exon_array'}) 
	 && (defined $transcript->{'_trans_exon_array'}->[$j])); 
	 $j++) {

      my $exon = $transcript->{'_trans_exon_array'}->[$j];

      for (my $i = 0; 
	   (($i <= scalar @{$exon->{'_supporting_evidence'}}) 
	    && (defined $exon->{'_supporting_evidence'}->[$i])); 
	   $i++){
	my $dna_align_feature = $exon->{'_supporting_evidence'}->[$i];

	my $gsf_seq = $gene->{'_transcript_array'}[0]{'_trans_exon_array'}[0]{'_gsf_seq'};

	# All DnaAlignFeatures for the same Gene/Transcript have the 
	# same Slice attached to them.  However, in the aftermath of
	# Data::Dumper all except the first attached Slice are mangled 
	# because they are cross-referenced.  This can be avoided by
	# attaching distinct duplicates before dumping.

	my $new_slice = Bio::EnsEMBL::Slice->new(-name          => $gsf_seq->name,
						 -chr_name      => $gsf_seq->chr_name, 
						 -chr_start     => $gsf_seq->chr_start, 
						 -chr_end       => $gsf_seq->chr_end,
						 -strand        => $gsf_seq->strand,
						 -assembly_type => $gsf_seq->assembly_type);

	$dna_align_feature->{'_gsf_seq'} = $new_slice;
	$dna_align_feature->{'_gsf_seq'}->{'adaptor'} = undef;
	$dna_align_feature->{'_analysis'}->{'_adaptor'} = undef;
      }
      my $old_exon_slice = $gene->{'_transcript_array'}[0]{'_trans_exon_array'}[0]{'_gsf_seq'};

      my $new_exon_slice = Bio::EnsEMBL::Slice->new(-name          => $old_exon_slice->name,
						    -chr_name      => $old_exon_slice->chr_name, 
						    -chr_start     => $old_exon_slice->chr_start, 
						    -chr_end       => $old_exon_slice->chr_end,
						    -strand        => $old_exon_slice->strand,
						    -assembly_type => $old_exon_slice->assembly_type);


      $exon->{'_gsf_seq'} = $new_exon_slice;
      $exon->{'_gsf_seq'}->{'adaptor'} = undef;
      $exon->{'adaptor'} = undef;

    }

  }
  my $gene_dump = Data::Dumper->Dump([$gene],['gene']);

$gene_dump =~ s/[\s\t\n]//g;

  my $storable_gene = $gene_dump;

  my $denorm_gene_dump_sql = q{INSERT INTO denormalised_gene_dump (denorm_gene_dump_id, dumped_gene)
			       VALUES ('/N', ?)};

  my $denorm_gene_sql = q{INSERT INTO denormalised_gene(denorm_gene_id, denorm_gene_dump_id, 
			chr_name, gene_start, gene_end, type)
			VALUES ('/N', ?, ?, ?, ?, ?)};

  my $sth1 = $self->prepare($denorm_gene_dump_sql);
  my $sth2 = $self->prepare($denorm_gene_sql);

  $sth1->execute($storable_gene);
  $sth2->execute($sth1->{'mysql_insertid'}, 
		 $gene->chr_name, 
		 $gene->start, 
		 $gene->end,  
		 $gene->type);

  return $sth2->{'mysql_insertid'};
}


sub get_genes_by_Slice_and_type {
  my ($self, $slice, $type) = @_;

  my $chr_name = $slice->chr_name;
  my $chr_start = $slice->chr_start;
  my $chr_end = $slice->chr_end;

  my $sql = q{SELECT dg.denorm_gene_dump_id, dgd.dumped_gene, gdg.gene_id
	      FROM denormalised_gene dg, denormalised_gene_dump dgd 
	      LEFT JOIN  gene_dumpedgene gdg 
	      ON dgd.denorm_gene_dump_id = gdg.denorm_gene_dump_id
	      WHERE dg.chr_name = ?
	      AND dg.gene_end >= ?
	      AND dg.gene_start <= ?
	      AND dg.type = ?
              AND dg.denorm_gene_dump_id = dgd.denorm_gene_dump_id};
  
  my $sth = $self->prepare($sql);

  $sth->execute($chr_name, $chr_start, $chr_end, $type);

  my @genes;

  # Fresh adaptors to append to objects:

  my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
  my $slice_adaptor = $self->db->get_SliceAdaptor;
  my $exon_adaptor = $self->db->get_ExonAdaptor;
  my $gene_adaptor = $self->db->get_GeneAdaptor;
  
  # Make an appropriate analysis object

  my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($EST_BLAT_ANALYSIS);

 DUMPED_GENE:
  while (my ($dumped_gene_id, $retrieved_gene, $already_mapped) = $sth->fetchrow_array) {

    next DUMPED_GENE if ($already_mapped);

    my $gene_dump = $retrieved_gene;


    my $gene;

    # Thaw out our stored gene objects:

    eval $gene_dump;

    # Tack on a fresh analysis object.

    $gene->analysis($analysis);

    # Replace the adaptor objects that we removed before 
    # dumping the gene objects.

    for (my $k = 0; 
	 (($k <= scalar @{$gene->{'_transcript_array'}}) 
	  && (defined $gene->{'_transcript_array'}->[$k]));
	 $k++){

      my $transcript = $gene->{'_transcript_array'}->[$k];

      for (my $j = 0; 
	   (($j <= scalar $transcript->{'_trans_exon_array'}) 
	    && (defined $transcript->{'_trans_exon_array'}->[$j]));
	   $j++) {
	my $exon = $transcript->{'_trans_exon_array'}->[$j];
	for (my $i = 0; 
	     (($i <= scalar @{$exon->{'_supporting_evidence'}}) 
	      && (defined $exon->{'_supporting_evidence'}->[$i])); 
	     $i++){
	  my $dna_align_feature = $exon->{'_supporting_evidence'}->[$i];
	  $dna_align_feature->{'_gsf_seq'}->{'adaptor'} = $slice_adaptor;
	  $dna_align_feature->{'_analysis'} = $analysis;
	}
	$exon->{'_gsf_seq'}->{'adaptor'} = $slice_adaptor;
	$exon->{'adaptor'} = $exon_adaptor;
      }
      
    }
    
    # We have to transform these genes to raw contig
    # coordinates sometime - now is a good time because
    # these will be retrieved by Slice, which means only
    # a small part of the assembly needs to be loaded.

    my $transformed_gene;

    eval{
      $transformed_gene = $gene->transform();
    };

    if ($@ || !defined $transformed_gene){
      $self->warn("could not transform coordinates of gene\n$@");
      foreach my $tran (@{$gene->get_all_Transcripts}){
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $tran );
      }
    }

    # Store our newly revived and remapped gene:

    my $gene_dbID = $gene_adaptor->store($transformed_gene);

    if ($gene_dbID) {
      my $lookup_sql = qq{INSERT INTO gene_dumpedgene (denorm_gene_dump_id, gene_id) 
			  VALUES ($dumped_gene_id, $gene_dbID)};

      my $update_sth = $self->db->prepare($lookup_sql);
      $update_sth->execute;

    }

  }

  $gene_adaptor->fetch_all_by_Slice($slice);

}


1;
__END__
