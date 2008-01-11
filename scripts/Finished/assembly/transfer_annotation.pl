#!/software/bin/perl

=head1 NAME

transfer_annotation.pl - transfer gene annotations across assemblies

=head1 SYNOPSIS

transfer_annotation.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  --assembly=ASSEMBLY                 assembly version ASSEMBLY
  --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional arguments:

  --chromosomes, --chr=LIST           only process LIST chromosomes
  --altchromosomes, --altchr=LIST     supply alternative chromosome names (the two lists must agree)
  --ref_start                         start coordinate on reference chromosomes
  --ref_end                           end coordinate on reference chromosomes
  --alt_start                         start coordinate on alternative chromosomes
  --alt_end                           end coordinate on alternative chromosomes
  --reverse                           use the reverse mapping order (compound to assembly)
  --haplotype                         for haplotype transfer, make a copy of genes (with new stable_id)

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION


=head1 AUTHOR

ml6@sanger.ac.uk

=head1 CONTACT

ml6@sanger.ac.uk

=cut

use strict;
use Bio::Vega::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport("$Bin/../../..");

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
	'assembly=s',         'altassembly=s',
	'chromosomes|chr=s@', 'altchromosomes|altchr=s@',
	'reverse!', 'haplotype!', 'ref_start=i', 'ref_end=i', 'alt_start=i', 'alt_end=i'
);
$support->allowed_params( $support->get_common_params, 'assembly',
	'altassembly', 'chromosomes', 'altchromosomes', );

if ( $support->param('help') or $support->error ) {
	print STDOUT $support->error if $support->error;
	pod2usage(1);
}

$support->comma_to_list( 'chromosomes', 'altchromosomes' );

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params( 'assembly', 'altassembly' );

# database connection
my $vega_db = $support->get_database('ensembl');
bless $vega_db, "Bio::Vega::DBSQL::DBAdaptor";
my $dbh = $vega_db->dbc->db_handle;

my $assembly    = $support->param('assembly');
my $altassembly = $support->param('altassembly');
my $haplo 		= $support->param('haplotype');

my $R_start = $support->param('ref_start') || undef;
my $R_end = $support->param('ref_end') || undef;
my $A_start = $support->param('alt_start') || undef;
my $A_end = $support->param('alt_end') || undef;

my $geneAd   = $vega_db->get_GeneAdaptor;
my $sliceAd  = $vega_db->get_SliceAdaptor;
my $sfeat_Ad = $vega_db->get_SimpleFeatureAdaptor();

# make sure that the coordinate system versions are different
# if they are the same (i.e. both Otter) create a temporary cs version MAPPING
my $sql_meta_insert = qq{
	insert into meta (meta_key, meta_value) values
	('assembly.mapping', 'chromosome:MAPPING|contig'),
	('assembly.mapping', 'chromosome:MAPPING|contig|clone'),};
if ( $support->param('reverse') ) {
	$sql_meta_insert .=
	  qq{('assembly.mapping', 'chromosome:MAPPING#chromosome:Otter')};
}
else {
	$sql_meta_insert .=
	  qq{('assembly.mapping', 'chromosome:Otter#chromosome:MAPPING')};
}
my $sql_meta_delete = qq{
	delete from meta where meta_value like 'chromosome%MAPPING%'};

my $sql_cs_insert = qq{
	insert into coord_system (coord_system_id, name, version, rank, attrib) values (100, 'chromosome', 'MAPPING', 100, '')};
my $sql_cs_delete = qq{
	delete from coord_system where version = 'MAPPING'};

my $sql_mc_insert = qq{
	insert into meta_coord (table_name, coord_system_id, max_length)
	values ('exon', 100, 1), ('gene', 100, 1), ('simple_feature', 100, 1), ('transcript', 100, 1)};
my $sql_mc_delete = qq{
	delete from meta_coord where coord_system_id = 100};

my $sql_sr_update = qq{
	update seq_region s, coord_system cs
	set s.coord_system_id = cs.coord_system_id
	where s.name = ?
	and cs.name = 'chromosome'
	and cs.version = ?};

my $sql_meta_bl_insert = qq{
	insert into meta (meta_key, meta_value) values
	('simple_featurebuild.level', '1'),
	('genebuild.level', '1')};
my $sql_meta_bl_delete =
qq{delete from meta where meta_key in ('simple_featurebuild.level','genebuild.level')};

my $sql_attrib_type_insert =
qq{insert into seq_region_attrib (seq_region_id, attrib_type_id, value) values (?, 6, 1)};
my $sql_attrib_type_delete =
qq{delete from seq_region_attrib where seq_region_id = ? and attrib_type_id = 6};

my $sth_attrib_type_insert = $dbh->prepare($sql_attrib_type_insert);
my $sth_attrib_type_delete = $dbh->prepare($sql_attrib_type_delete);

my $cs_change = 0;

if ( $assembly eq $altassembly ) {
	$cs_change   = 1;
	$altassembly = 'MAPPING';

	$dbh->do($sql_meta_insert);
	$dbh->do($sql_cs_insert);
	$dbh->do($sql_mc_insert);
	$dbh->do($sql_meta_bl_insert);
}

my @R_chr_list = $support->param('chromosomes');
die "Reference chromosome list must be defined!" unless scalar(@R_chr_list);

my @A_chr_list = $support->param('altchromosomes');
if ( scalar(@R_chr_list) != scalar(@A_chr_list) ) {
	die "Chromosome lists do not match by length";
}

SET: for my $i ( 0 .. scalar(@R_chr_list) - 1 ) {
	my $R_chr = $R_chr_list[$i];
	my $A_chr = $A_chr_list[$i];

	# for stats
	my $total_genes      = 0;
	my $transfered_genes = 0;
	my $total_sf         = 0;
	my $transfered_sf    = 0;
	my $skipped_sf       = 0;

	$dbh->begin_work;

	eval {
		$support->log_verbose("Annotation transfer $R_chr => $A_chr...\n");
		my $sth_cs = $dbh->prepare($sql_sr_update);
		$sth_cs->execute( $A_chr, $altassembly ) unless ( !$cs_change );

		my $ref_sl =
	  		$sliceAd->fetch_by_region( 'chromosome', $R_chr, $R_start, $R_end, undef,$assembly );
		my $alt_sl =
	  		$sliceAd->fetch_by_region( 'chromosome', $A_chr, $A_start, $A_end, undef,$altassembly );

		my $alt_seq_region_id = $alt_sl->get_seq_region_id;
		my $ref_seq_region_id = $ref_sl->get_seq_region_id;
		map ( $sth_attrib_type_insert->execute($_),
			( $alt_seq_region_id, $ref_seq_region_id ) );

		my @genes;
		@genes       = @{ $ref_sl->get_all_Genes };
		$total_genes = scalar @genes;
		@genes = sort { $a->start <=> $b->start } @genes;

		# transfer simple features (PolyA_site/_signal)
		my @proj_feat;
		my @simple_features = @{ $sfeat_Ad->fetch_all_by_Slice($ref_sl) };
		$total_sf = scalar @simple_features;
		foreach my $f (@simple_features) {
			my $tf = $f->transfer($alt_sl);
			if ( defined $tf ) {
				$tf->dbID(undef);
				$tf->adaptor(undef);
				if ( can_save( $sfeat_Ad, $tf ) ) {
					push @proj_feat, $tf;
					$transfered_sf++;
				}
				else {
					$support->log_verbose(
						sprintf(
							"SKIP %s %d %d %d, already saved on $A_chr\n",
							$tf->display_label, $tf->start,
							$tf->end,           $tf->strand
						)
					);
					$skipped_sf++;
				}

			}
			else {
				warning(
					sprintf(
						" %s %d %d %d cannot be transfered on $A_chr\n",
						$f->display_label, $f->start, $f->end, $f->strand
					)
				);

			}

		}

	  GENE: foreach my $g (@genes) {
			my $tg = Bio::Vega::Gene->new;
			$tg->analysis( $g->analysis );
			$tg->biotype( $g->biotype );
			$tg->status( $g->status );
			$tg->source( $g->source );
			$tg->stable_id( $g->stable_id ) unless $haplo;
			$tg->gene_author( $g->gene_author );
			$tg->description( $g->description );

			$support->log_verbose(
				sprintf(
					"Transferring gene %s %s %d %d\n",
					$g->stable_id, $g->biotype, $g->start, $g->end
				)
			);

			my @proj_trans;

		  TRANSCRIPT: foreach my $t ( @{ $g->get_all_Transcripts } ) {
				my $tt = $t->transfer($alt_sl);    # ori
				if ( defined $tt ) {
					push @proj_trans, $tt;
					$support->log_verbose(
						sprintf(
							"\t%s %s %d %d transferred successfully:\n",
							$t->stable_id, $t->biotype, $t->start, $t->end
						)
					);
					$tt->status( $t->status );
					&remove_all_db_ids($tt);

					&log_compare_transcripts( $t, $tt );
				}
				else {
					warning(
						sprintf(
							"\t%s %s %d %d cannot be transfered on $A_chr\n",
							$t->stable_id, $t->biotype, $t->start, $t->end
						)
					);
				}
			}

			$support->log_verbose(
				sprintf(
					"\tSummary for %s : %d out of %d transcripts transferred\n",
					$g->stable_id,
					scalar(@proj_trans),
					scalar( @{ $g->get_all_Transcripts } )
				)
			);

			if (@proj_trans) {

			   # essential for loading attribute list so that get_all_Attributes
			   # won't return empty list
				$tg->add_Attributes( @{ $g->get_all_Attributes() } );

				map( $tg->add_Transcript($_), @proj_trans );

				&write_gene($tg);

				if ( $geneAd->store($tg) ) {
					$transfered_genes++;
					$support->log_verbose(
						sprintf(
"GENE %s %s successfully TRANSFERED (%s:%d-%d => %s:%d-%d)\n",
							$tg->stable_id, $tg->biotype, $R_chr,
							$g->start,      $g->end,      $A_chr,
							$tg->start,     $tg->end
						)
					);
				}
				else {
					throw(
						sprintf(
							"GENE %s %s cannot be saved (%s:%d-%d => %s:%d-%d)",
							$tg->stable_id, $tg->biotype, $R_chr,
							$g->start,      $g->end,      $A_chr,
							$tg->start,     $tg->end
						)
					);
				}
			}
		}

		# save polyA features if gene transfered
		$sfeat_Ad->store(@proj_feat) if ( @proj_feat && $transfered_genes );

		$sth_cs->execute( $A_chr, $support->param('altassembly') )
		  unless ( !$cs_change );
		map ( $sth_attrib_type_delete->execute($_),
			( $alt_seq_region_id, $ref_seq_region_id ) );

		$dbh->commit;
	};

	if ($@) {
		$dbh->rollback;
		warning(
			    "UNABLE TO TRANSFER ANNOTATION FROM $R_chr:$assembly to $A_chr:"
			  . $support->param('altassembly')
			  . "\n[$@]\n" );
		next SET;
	}

	# print annotation transfer stats
	$support->log_verbose(
		sprintf(
"Annotation transfer %s:%s => %s:%s\ntransfered Gene: %d/%d\ntransfered PolyA features: %d/%d\nskipped PolyA features: %d/%d\n",
			$R_chr,            $assembly,
			$A_chr,            $support->param('altassembly'),
			$transfered_genes, $total_genes,
			$transfered_sf,    $total_sf,
			$skipped_sf,       $total_sf
		)
	);
}

# remove temporary 'MAPPING' cs from meta and cs table
if ($cs_change) {
	$dbh->do($sql_meta_delete);
	$dbh->do($sql_cs_delete);
	$dbh->do($sql_mc_delete);
	$dbh->do($sql_meta_bl_delete);
}

# check if feature is already on alt slice
sub can_save {
	my ( $sfeat_Ad, $f ) = @_;
	my $a =
	  $sfeat_Ad->fetch_all_by_Slice( $f->feature_Slice,
		$f->analysis->logic_name );

	return @$a ? 0 : 1;
}

sub remove_all_db_ids {
	my $t = shift;

	$t->dbID(undef);
	$t->stable_id(undef) if $haplo;
	$t->adaptor(undef);
	if ( $t->translation ) {
		$t->translation->dbID(undef);
		$t->translation->adaptor(undef);
	}
	foreach my $e ( @{ $t->get_all_Exons } ) {
		$e->dbID(undef);
		$e->adaptor(undef);
	}

	return $t;
}

sub transcript_is_missed {
	my ( $t, $alt_sl ) = @_;

	foreach my $e ( @{ $t->get_all_Exons } ) {
		my $alt_e = $e->transfer($alt_sl);
		if ( $alt_e && scalar @{$alt_e} ) {
			return 0;
		}
	}

	return 1;
}

sub get_coding_part_of_transcript {
	my ($t) = @_;

	my @cds = @{ $t->get_all_translateable_Exons };
	print STDOUT $t;
	print STDOUT $t->stable_id;
	foreach my $e ( @{ $t->get_all_Exons } ) {

		#    print STDOUT "EXON: $e", " ", $e->stable_id;
	}
	$t->{_trans_exon_array} = [];
	foreach my $e (@cds) {

		#   print STDOUT "CODING $e", " ", $e->stable_id;
		$t->add_Exon($e);
	}

	my $tr = $t->translation;

	#map { print STDOUT "$_ ", $tr->{$_} } keys %$tr; die;

	$tr->start_Exon( $cds[0] );
	$tr->start(1);
	$tr->end_Exon( $cds[-1] );
	$tr->end( $tr->end_Exon->length );
	$t->translation($tr);

	$t->stable_id( $t->stable_id . "_CDS" );
	return $t;

	#  #my $cds_t = Bio::EnsEMBL::Transcript->new();
##  my $cds_t = Bio::Vega::Transcript->new();

	#  for(my $i=0; $i < @cds; $i++) {
	#    $cds_t->add_Exon($cds[$i]);
	#  }

	#  $cds_t->analysis($t->analysis);
	#  $cds_t->biotype($t->biotype);
	#  $cds_t->slice($t->slice);
	#  $cds_t->status($t->status);
	# # my $tr = Bio::EnsEMBL::Translation->new;
	#  my $tr = Bio::Vega::Translation->new;

	#  $tr->start_Exon($cds[0]);
	#  $tr->start(1);
	#  $tr->end_Exon($cds[-1]);
	#  $tr->end($tr->end_Exon->length);
	#  $cds_t->translation($tr);

	#  $cds_t->stable_id($t->stable_id . "_CDS");

	#  return $cds_t;
}

sub project_transcript_the_hard_way {
	my ( $t, $cs_name, $cs_version, $alt_sl ) = @_;

	#map { print STDOUT "$_: ", $t->{$_} } keys %$t;die;
	#  my $new_t = Bio::EnsEMBL::Transcript->new();
	#  my $new_t = Bio::Vega::Transcript->new();

	#  $new_t->analysis($t->analysis);
	#  $new_t->stable_id($t->stable_id);
	#  $new_t->biotype($t->biotype);
	#  $new_t->status($t->status);

	my @new_e;
	foreach my $e ( @{ $t->get_all_Exons } ) {
		my $te = $e->transform($alt_sl);

		if ( defined $te ) {
			push @new_e, $te;
		}
		else {
			my @bits = @{ $e->project( $cs_name, $cs_version ) };

			if (@bits) {

				# need to make a new exon from the bits
				my ( $ex_st, $ex_en, $strand );
				foreach my $bit (@bits) {
					if ( not defined $ex_st
						or $bit->to_Slice->start < $ex_st )
					{
						$ex_st = $bit->to_Slice->start;
					}
					if ( not defined $ex_en
						or $bit->to_Slice->end > $ex_en )
					{
						$ex_en = $bit->to_Slice->end;
					}
					if ( not defined $strand ) {
						$strand = $bit->to_Slice->strand;
					}
				}
				if ( defined $ex_st and defined $ex_en ) {

					#my $new_e = Bio::EnsEMBL::Exon->new;
					my $new_e = Bio::Vega::Exon->new;

					$new_e->start($ex_st);
					$new_e->end($ex_en);
					$new_e->strand($strand);
					$new_e->phase(-1);
					$new_e->end_phase(-1);
					$new_e->slice($alt_sl);
					$new_e->stable_id( $e->stable_id );
					push @new_e, $new_e;
				}
			}
		}
	}

	$t->{_trans_exon_array} = [];
	foreach my $e (@new_e) {
		$t->add_Exon($e);
	}
	return $t;

	#  map { $new_t->add_Exon($_) } @new_e;
	#  return $new_t;
}

sub log_summarise_projection_by_exon {
	my ( $t, $cs_name, $cs_version ) = @_;

	my $exons_total       = 0;
	my $exons_transferred = 0;
	my $exons_missed      = 0;
	my $exons_part_missed = 0;
	my $exons_over_indel  = 0;

	foreach my $e ( @{ $t->get_all_Exons } ) {
		$exons_total++;

		my $te = $e->transform( $cs_name, $cs_version );

		if ( defined $te ) {
			$exons_transferred++;
		}
		else {
			my @bits = @{ $e->project( $cs_name, $cs_version ) };

			if ( not @bits ) {
				$exons_missed++;
			}
			elsif ($bits[0]->from_start != 1
				or $bits[-1]->from_end != $e->length )
			{
				$exons_part_missed++;
			}
			else {
				$exons_over_indel++;
			}
		}
	}

	$support->log_verbose(
		sprintf(
" %d exons; %d transferred, %d missed, %d part missed, %d over indel\n",
			$exons_total,       $exons_transferred, $exons_missed,
			$exons_part_missed, $exons_over_indel
		)
	);
}

sub log_compare_transcripts {
	my ( $old, $new ) = @_;

	my $cdnaseq_old = $old->spliced_seq;
	my $cdnaseq_new = $new->spliced_seq;

	my @exons_old = @{ $old->get_all_Exons };
	my @exons_new = @{ $new->get_all_Exons };

	my ( $pepseq_old, $pepseq_new );
	if ( $old->translation ) {
		$pepseq_old = $old->translate;
	}
	eval {
		if ( $new->translation )
		{
			$pepseq_new = $new->translate;
		}
	};

	unless ($@) {
		if ( $cdnaseq_old eq $cdnaseq_new ) {
			$support->log_verbose("\t\ttransfer to identical DNA\n");
		}
		else {
			my $cdna_diffs = &compare_seqs( $cdnaseq_old, $cdnaseq_new );
			$support->log_verbose("\t\ttransfer to $cdna_diffs cDNA diffs");
			if ( defined $pepseq_old and defined $pepseq_new ) {
				my $prot_diffs =
				  &compare_seqs( $pepseq_old->seq, $pepseq_new->seq );
				$support->log_verbose(" and $prot_diffs protein diffs");
				my $pepseq_new_seq = $pepseq_new->seq
				  || die "something's wrong with the sequence";
				my $stop_count = $pepseq_new_seq =~ tr/\*/\*/;
				if ($stop_count) {
					$support->log_verbose(" ($stop_count STOPS)");
				}
			}
			$support->log_verbose("\n");
		}
	}
}

sub compare_seqs {
	my ( $seq1, $seq2 ) = @_;

	my @seq1 = split( //, $seq1 );
	my @seq2 = split( //, $seq2 );

	my $total = 0;
	my $diffs = 0;

	for ( my $i = 0 ; $i < @seq1 ; $i++ ) {
		if ( $seq1[$i] ne $seq2[$i] ) {
			$diffs++;
		}
		$total++;
	}

	return $diffs;
}

sub write_gene {
	my $g = shift;

	# seq source type start end score strand frame group
	my $seq_id = $g->slice->seq_region_name;

	printf( "\t$seq_id\tgene\t%s\t%s\t%s\n",
		$g->stable_id, $g->biotype, $g->status );

	my %exon_hash;

	foreach my $t ( @{ $g->get_all_Transcripts } ) {

		foreach my $e ( @{ $t->get_all_Exons } ) {
			my $hash_key = sprintf( "\t%d:%d:%d:%d",
				$e->start, $e->end, $e->strand, $e->phase, $e->end_phase );
			push @{ $exon_hash{ $e->stable_id }->{$hash_key} }, $e;
		}
	}
	foreach my $sid ( keys %exon_hash ) {
		if ( scalar( keys %{ $exon_hash{$sid} } ) > 1 ) {

			# different projections for same exon; we
			# need to give them all distinct stable ids
			my $count = 1;
			foreach my $key ( keys %{ $exon_hash{$sid} } ) {
				my $elist = $exon_hash{$sid}->{$key};

				foreach my $mem (@$elist) {
					$mem->stable_id( $mem->stable_id . "_" . $count );
				}
				$count++;
			}
		}
		foreach my $key ( keys %{ $exon_hash{$sid} } ) {
			my ($e) = @{ $exon_hash{$sid}->{$key} };

			printf( "\t$seq_id\texon\t%s\t%d\t%d\t%d\t%d\t%d\n",
				$e->stable_id, $e->start, $e->end, $e->strand, $e->phase,
				$e->end_phase );
		}
	}

	foreach my $t ( @{ $g->get_all_Transcripts } ) {
		printf( "\t$seq_id\ttranscript\t%s\t%s\t%s\n",
			$t->stable_id, $t->biotype, $t->status );

		foreach my $e ( @{ $t->get_all_Exons } ) {
			printf( "\t$seq_id\texon_transcript\t%s\t%s\n",
				$e->stable_id, $t->stable_id );
		}
		eval {
			if ( $t->translation )
			{
				printf(
					"\t$seq_id\ttranslation\t%s\t%d\t%s\t%d\n",
					$t->translation->start_Exon->stable_id,
					$t->translation->start,
					$t->translation->end_Exon->stable_id,
					$t->translation->end
				);

			}
		};
	}
}

