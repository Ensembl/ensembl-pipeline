=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::BlatGeneAdaptor

=head1 SYNOPSIS

$blat_genes = Bio::EnsEMBL::Pipeline::DBSQL::BlatGeneAdaptor->new($db) ;

$blat_genes->store($gene);
my $genes = $blat_genes->fetch_all_by_Slice($slice);

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


package Bio::EnsEMBL::Pipeline::DBSQL::BlatGeneAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_INPUTID_REGEX
				       );


use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 new

  Arg [1]    : list of arguments @args
  Example    : $blat_gene_adaptor = Bio::EnsEMBL::Pipeline::DBSQL::BlatGeneAdaptor->new($db);
  Description: Creates a new BlatGeneAdaptor object
  Returntype : Bio::EnsEMBL::Pipeline::DBSQL::BlatGeneAdaptor
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


=head2 store

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $blat_gene_adaptor->store($gene);
  Description: stores a gene in the special BLAT gene tables.
  Returntype : the database identifier of the newly stored gene
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene or if 
               $gene does not have an analysis object
  Caller     : general

=cut

sub store {
   my ($self,$gene) = @_;

   if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Must store a gene object, not a $gene");
   }
   if( !defined $gene->analysis ) {
       $self->throw("Genes must have an analysis object!");
   }

   my $aA = $self->db->get_AnalysisAdaptor();
   my $analysisId = $aA->exists( $gene->analysis() );

   if( defined $analysisId ) {
     $gene->analysis()->dbID( $analysisId );
   } else {
     $analysisId = $self->db->get_AnalysisAdaptor()->store( $gene->analysis );
   }

   if ( !defined $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not writing!");
   }
 
   my $trans_count = scalar(@{$gene->get_all_Transcripts});
   my $type = "";

   if (defined($gene->type)) {
       $type = $gene->type;
   }

   # assuming that the store is used during the EST Genebuild process, set
   # the display_xref_id to 0.  This ought to get re-set during the protein
   # pipeline run.  This probably update to the gene table has yet to be
   # implemented.
   my $xref_id = 0;

   my $sth2 = $self->prepare("INSERT INTO blat_gene(type, analysis_id, 
                                               transcript_count, display_xref_id) 
                              VALUES('$type', $analysisId, $trans_count, $xref_id)" );

   $sth2->execute();

   my $gene_dbID = $sth2->{'mysql_insertid'};


   # write exons transcripts and exon_transcript table

   my $trans = $gene->get_all_Transcripts;

   #force lazy loading of translations before new exon dbIDs are set
   map {$_->translation} @$trans;

   foreach my $t ( @$trans ) {
     $self->_store_Transcript($t,$gene_dbID );
   }

   $gene->adaptor( $self );
   $gene->dbID( $gene_dbID );
   
   return $gene->dbID;
}


=head2 _store_Transcript

 Title   : _store_Transcript
 Usage   : $self->_store_Transcript( $transcript )
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    : needs a gene ...


=cut

sub _store_Transcript {
   my ($self,$transcript,$gene_dbID) = @_;

   if( ! ref $transcript || !$transcript->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$transcript is not a EnsEMBL transcript - not dumping!");
   }

   # store the transcript
   # then store the exon_transcript table

   my ( $exon_count, $exons );
   $exons = $transcript->get_all_Exons();
   $exon_count = scalar( @{$exons} );

   # assuming that the store is used during the EST Genebuild process, set
   # the display_xref_id to 0.  This ought to get re-set during the protein
   # pipeline run.  This probably update to the gene table has yet to be
   # implemented.
   my $xref_id = 0;

   # ok - now load this line in
   my $tst = $self->prepare("
        insert into blat_transcript ( gene_id, translation_id, exon_count, display_xref_id )
        values ( ?, ?, ?, 0)
        ");

   $tst->execute( $gene_dbID, 0, $exon_count );


   my $transc_dbID = $tst->{'mysql_insertid'};

   my $etst = $self->prepare("insert into blat_exon_transcript (exon_id,transcript_id,rank) values (?,?,?)");
   my $rank = 1;
   foreach my $exon (@$exons) {
     my $exon_dbID = $self->_store_Exon( $exon );
     $etst->execute($exon_dbID,$transc_dbID,$rank);
     $rank++;
   }

   $transcript->dbID( $transc_dbID );
   $transcript->adaptor( $self );
   return $transc_dbID;
}


=head2 _store_Exon

  Arg [1]    : Bio::EnsEMBL::Exon $exon
               the exon to store in this database
  Example    : $self->_store_Exon($exon);
  Description: Stores an exon in the database
  Returntype : none
  Exceptions : thrown if exon (or component exons) do not have a contig_id
               or if $exon->start, $exon->end, $exon->strand, or $exon->phase 
               are not defined or if $exon is not a Bio::EnsEMBL::Exon 
  Caller     : general

=cut

sub _store_Exon {
  my ( $self, $exon ) = @_;

  if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
    $self->throw("$exon is not a EnsEMBL exon");
  }

  if( ! $exon->start || ! $exon->end ||
      ! $exon->strand || ! defined $exon->phase ) {
    $self->throw("Exon does not have all attributes to store");
  }

  # trap contig_id separately as it is likely to be a common mistake

  my $exon_sql = q{
    INSERT into blat_exon ( exon_id, chr_id, chr_start, 
		       chr_end, chr_strand, phase, 
		       end_phase, sticky_rank )
    VALUES ( ?, ?, ?, ?, ?, ?, ?,? ) 
  };
  my $exonst = $self->prepare($exon_sql);

  my $exonId = undef;

  if( $exon->isa( 'Bio::EnsEMBL::StickyExon' )) {
    # sticky storing. Sticky exons contain normal exons ...

    my $componentExons = $exon->get_all_component_Exons;
    for my $componentExon ( @$componentExons ) {
      my $slice = $componentExon->contig();

      unless($slice && ref $slice && $slice->name) {
	$self->throw("Component Exon does not have an attached Slice " .
		     "with a valid name. Needs to have one set");
      }

      my $chr_name;
      my $chr_offset;
      if ($slice->name =~ /$EST_INPUTID_REGEX/){
	$chr_name = $1;
	$chr_offset = $2;
      } else {
	$self->throw("Cant parse chromosome input id with regular expression " . 
		     "provided in EST_INPUTID_REGEX in Bio::EnsEMBL::Pipeline::ESTConf.pm");
      }
      
    
      my $chr_start = $chr_offset + $componentExon->start();
      my $chr_end = $chr_offset + $componentExon->end();


      $exonst->execute( $exonId, $chr_name,
			$chr_start,
			$chr_end,
			$componentExon->strand(),
			$componentExon->phase(),
			$componentExon->end_phase(),
			$componentExon->sticky_rank() );
      if( ! $exonId ) {
	$exonId = $exonst->{'mysql_insertid'};
      }
    }
  } else {
    # normal storing

    my $slice = $exon->contig();

    unless( $slice && ref $slice && $slice->name ) {
      $self->throw("Exon does not have an attached Slice with a valid " . 
		   "name.  Needs to have one set");
    }

    my $chr_name;
    my $chr_offset;
    if ($slice->name =~ /$EST_INPUTID_REGEX/){
      $chr_name = $1;
      $chr_offset = $2;
    } else {
      $self->throw("Cant parse chromosome input id with regular expression " . 
		   "provided in EST_INPUTID_REGEX in Bio::EnsEMBL::Pipeline::ESTConf.pm");
    }

    
    my $chr_start = $chr_offset + $exon->start();
    my $chr_end = $chr_offset + $exon->end();

    $exonst->execute( undef,$chr_name,
		      $chr_start,
		      $chr_end,
		      $exon->strand(),
		      $exon->phase(),
		      $exon->end_phase(),
		      $exon->sticky_rank() );
    $exonId = $exonst->{'mysql_insertid'};
  }

  # Now the supporting evidence
  # should be stored from featureAdaptor
  my $sql = "insert into blat_supporting_feature (exon_id, feature_id, feature_type)
             values(?, ?, ?)";  
  
  my $sf_sth = $self->db->prepare($sql);

  my $anaAdaptor = $self->db->get_AnalysisAdaptor();
  my $type;

  my @exons = ();
  if($exon->isa('Bio::EnsEMBL::StickyExon')) {
    @exons = @{$exon->get_all_component_Exons};
  } else {
    @exons = ($exon);
  }

  foreach my $e (@exons) {
    foreach my $sf (@{$e->get_all_supporting_features}) {
      unless($sf->isa("Bio::EnsEMBL::BaseAlignFeature")){
	$self->throw("$sf must be an align feature otherwise" .
		     "it can't be stored");
      }

      #sanity check
      eval { $sf->validate(); };
      if ($@) {
        $self->warn("Supporting feature invalid. Skipping feature\n");
	next;
      }

#      $sf->contig($e->contig);

      my $sf_dbID;
      if($sf->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
	$sf_dbID = $self->_store_DnaDnaAlignFeature($sf);
	$type = 'dna_align_feature';
      } else {
	$self->warn("Supporting feature is not a DnaDnaAlignFeature. Skipping : [$sf]\n");
	next;
      }

      $sf_sth->execute($exonId, $sf_dbID, $type);
    }
  }

  #
  # Finally, update the dbID and adaptor of the exon (and any component exons)
  # to point to the new database
  #
  foreach my $e (@exons) {
    $e->dbID($exonId);
    $e->adaptor($self);
  }

  $exon->adaptor($self);
  $exon->dbID($exonId);
}

=head2 _store_DnaDnaAlignFeature
 
  Arg [1]    : list of Bio::EnsEMBL::DnaAlignFeatures @sf
               the features to store in the database
  Example    : $self->_store_DnaDnaAlignFeature(@features);
  Description: Stores a list of DnaDnaAlignFeatures in the database
  Returntype : none
  Exceptions : 
  Caller     : ?

=cut

sub _store_DnaDnaAlignFeature {
  my ($self, $sf) = @_;

  my $sth = $self->prepare("
     INSERT INTO blat_dna_align_feature (chr_id, chr_start, chr_end,
                             chr_strand, hit_start, hit_end,
                             hit_strand, hit_name, cigar_line,
                             analysis_id, score, evalue, perc_ident) 
     VALUES (?,?,?,?,?,?,?,?,?,?,?, ?, ?)");


  if( !ref $sf || !$sf->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
    $self->throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature," 
		 . " not a [$sf]");
  }
  
  my $slice = $sf->entire_seq();
  unless(defined $slice && $slice->isa("Bio::EnsEMBL::Slice")) {
    $self->throw("A Slice must be attached to the features to be " .
		 "stored via the attach seq method\n");
  }

  my $chr_name;
  my $chr_offset;
  if ($slice->name =~ /$EST_INPUTID_REGEX/){
    $chr_name = $1;
    $chr_offset = $2;
  } else {
    $self->throw("Cant parse chromosome input id with regular expression " . 
		 "provided in EST_INPUTID_REGEX in Bio::EnsEMBL::Pipeline::ESTConf.pm");
  }
  
  
  my $chr_start = $chr_offset + $sf->start();
  my $chr_end = $chr_offset + $sf->end();
  
  if( !defined $sf->analysis ) {
    $self->throw("Cannot store sequence features without analysis");
  }
  
  $sth->execute( $chr_name, $chr_start, $chr_end, $sf->strand,
		 $sf->hstart, $sf->hend, $sf->hstrand, $sf->hseqname,
		 $sf->cigar_string, $sf->analysis->dbID, $sf->score, 
		 $sf->p_value, $sf->percent_id);
  $sf->dbID($sth->{'mysql_insertid'});
}

#writing
###########################################################
# reading


=head2 fetch_all_Genes_by_Slice

  Arg [1]    :                
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub fetch_all_Genes_by_Slice {



}


=head2 fetch_by_dbID

  Arg [1]    : int $geneId 
               the unique internal database id of the Gene to be retrieved
  Arg [2]    : int $chromosomal_coordinates (optional)
               if defined, try to return chromosomal coordinates.
  Example    : $gene = $gene_adaptor->fetch_by_dbID
  Description: Retrieves a gene object from the database using its unique
               internal identifier.
  Returntype : Bio::EnsEMBL::Gene in contig coordinates
  Exceptions : thrown if no exons exist for the gene with dbID $geneId
  Caller     : general

=cut


sub fetch_by_dbID {
  my ( $self, $geneId, $chr_coordinates ) = @_;

  my $exonAdaptor = $self->db->get_ExonAdaptor();
  my @exons = @{$exonAdaptor->fetch_all_by_gene_id( $geneId )};
  my %exonIds;

  #
  # We assumme that if we get no exons this is a throwing situation
  #
  if( scalar( @exons ) == 0 ) {
    $self->throw("No exons for gene $geneId, assumming no gene");
  }

  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  
  # fetching all exons by gene
  # fetching all transcripts
  # adding the exons
  # adding the transcripts
  my %transcriptExons;
  my %transcripts;

  my $query = qq{
    SELECT tscript.gene_id
      , tscript.transcript_id
      , e_t.exon_id, e_t.rank
      , gene.analysis_id
      , gene.type
      , tscript.translation_id
    FROM gene
      , transcript tscript
      , exon_transcript e_t
    WHERE gene.gene_id = tscript.gene_id
      AND tscript.transcript_id = e_t.transcript_id
      AND gene.gene_id = $geneId

    ORDER BY tscript.gene_id
      , tscript.transcript_id
      , e_t.rank
    };

  my $sth = $self->prepare( $query );
  $sth->execute();
  my $ana;
  my $first = 1;
  my $gene;

  while( my @arr = $sth->fetchrow_array() ) {
    # building a gene
    if( $first ) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->adaptor($self);
      $gene->dbID( $geneId );
      $ana = $self->db->get_AnalysisAdaptor->fetch_by_dbID($arr[4]);
      $gene->analysis($ana);
      $gene->type($arr[5]);
      $first = 0;
    }

    # store an array of exon ids for each transcript
    if( !exists $transcriptExons{$arr[1]} ) {
      $transcriptExons{$arr[1]} = [];
    }

    push( @{$transcriptExons{$arr[1]}}, $arr[2] );
    $transcripts{$arr[1]} = $arr[6];
  }

  if( $first ) {
    return undef;
  }
  
  #discard duplicate exons, add analysis object to exons
  foreach my $exon ( @exons ) {
    $exon->analysis($ana);
    $exonIds{$exon->dbID} = $exon;
  }
  
  foreach my $transcriptId ( keys %transcripts ) {
    # should be fetch_by_geneId ..
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->dbID( $transcriptId );
    $transcript->adaptor( $self->db->get_TranscriptAdaptor() );
    $transcript->_translation_id($transcripts{$transcriptId} );

    foreach my $exonId ( @{$transcriptExons{$transcriptId}} ) {
      
      $transcript->add_Exon( $exonIds{$exonId} );
    }
    # store also a type for the transcript
    $transcript->type($gene->type);
    $gene->add_Transcript( $transcript );
  }
  
  # if chromosomal coordinates are needed, transform with empty slice

  if( defined $chr_coordinates ) {
    my $sa = $self->db->get_SliceAdaptor();
    my $empty_slice = Bio::EnsEMBL::Slice->new
      ( 
       -empty => 1,
       -adaptor => $sa 
      );
    $gene->transform( $empty_slice );
  }

  return $gene;
}


=head2 fetch_all_Exons_by_gene_id

  Arg [1]    : int $id
               The identifier of the gene whose exons will be retrieved 
  Example    : @exons = $exon_adaptor->fetch_all_by_gene_id(1234); 
  Description: Retrieves all exons from the gene specified by $geneId
  Returntype : listref of Bio::EnsEMBL::Exon in contig coordinates
  Exceptions : thrown if $geneId is not defined  
  Caller     : general

=cut

sub fetch_all_Exons_by_gene_id {
  my ( $self, $gene_id ) = @_;
  my %exons;
  my $hashRef;
  my ( $currentId, $currentTranscript );

  if( !$gene_id ) {
      $self->throw("Gene dbID not defined");
  }
  $self->{rchash} = {};
  my $query = qq {
    SELECT  e.exon_id
      , e.contig_id
      , e.contig_start
      , e.contig_end
      , e.contig_strand
      , e.phase
      , e.end_phase
      , e.sticky_rank
    FROM blat_exon e
      , blat_exon_transcript et
      , blat_transcript t
    WHERE t.gene_id = ?
      AND et.transcript_id = t.transcript_id
      AND e.exon_id = et.exon_id
    ORDER BY t.transcript_id,e.exon_id
      , e.sticky_rank DESC
  };

  my $sth = $self->prepare( $query );
  $sth->execute($gene_id);

  while( $hashRef = $sth->fetchrow_hashref() ) {
    if( ! exists $exons{ $hashRef->{exon_id} } ) {

      my $exon = $self->_exon_from_sth( $sth, $hashRef );

      $exons{$exon->dbID} = $exon;
    }
  }
  delete $self->{rchash};
  
  my @out = ();

  push @out, values %exons;

  return \@out;
}










1;
__END__
