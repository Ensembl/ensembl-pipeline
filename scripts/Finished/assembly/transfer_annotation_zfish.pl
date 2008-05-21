#!/software/bin/perl

=head1 NAME

transfer_annotation_zfish.pl - transfer gene annotations from one set of zfish assemblies to another, step 4

=head1 SYNOPSIS

transfer_annotation_zfish.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  --from_assembly                     old assembly date
  --to_assembly                       new assembly date
  --author                            author login used to lock the assemblies

Optional arguments:


  --dry=0|1                           don't write results to database (default: false)

  --assembly=ASSEMBLY                 assembly version ASSEMBLY (default: Otter)
  --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY (default: assembly)
  --email                             email used to lock the assemblies (default: author)
  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: true)
  -i, --interactive=0|1               run script interactively (default: false)
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION


=head1 RELATED FILES

The whole process of creating a whole genome alignment between two sets of zfish assemblies
is done by a series of scripts. Please see scripts in

  ensembl-pipeline/scripts/Finished/assembly/

for a high-level description of this process, and POD in the individual scripts
for the details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

modified by Mustapha Larbaoui <ml6@sanger.ac.uk>

=head1 CONTACT

Please post comments/questions to Anacode
<anacode-people@sanger.ac.uk>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);
BEGIN {
    $SERVERROOT = "$Bin/../../../..";
    unshift(@INC, "$Bin");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "/software/anacode/otter/otter_production_main/ensembl-otter/modules/");
    unshift(@INC, "$SERVERROOT/bioperl-0.7.2");
    unshift(@INC, "$SERVERROOT/bioperl-1.2.3-patched");
}

use Getopt::Long;
use Pod::Usage;
use Sys::Hostname;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Attribute;
use Bio::Vega::DBSQL::DBAdaptor;
use Bio::Vega::ContigLockBroker;
use Bio::Vega::Author;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
	'from_assembly=s', 'to_assembly=s',	'assembly=s', 'altassembly=s',
	'author=s', 'email=s'
);
$support->allowed_params( $support->get_common_params, 'from_assembly=s', 'to_assembly=s', 'assembly=s','altassembly=s' );

if ( $support->param('help') or $support->error ) {
	print STDOUT $support->error if $support->error;
	pod2usage(1);
}

$support->param('verbose', 1);
$support->param('interactive', 0);

# get log filehandle and print heading and parameters to logfile
$support->init_log;


# database connection
my $vega_db = $support->get_database('ensembl');
bless $vega_db, "Bio::Vega::DBSQL::DBAdaptor";
my $dbh = $vega_db->dbc->db_handle;

my $from_assembly = $support->param('from_assembly');
my $to_assembly = $support->param('to_assembly');

my $assembly    = $support->param('assembly') || 'Otter';
my $altassembly = $support->param('altassembly') || $assembly;

my $write_db      = not $support->param('dry_run');
my $author 		= $support->param('author');
my $email 		= $support->param('email') || $author;


throw("must set author name to lock the assemblies") if (!$author);


my $geneAd   = $vega_db->get_GeneAdaptor;
my $sliceAd  = $vega_db->get_SliceAdaptor;
my $sfeat_Ad = $vega_db->get_SimpleFeatureAdaptor();

# make sure that the coordinate system versions are different
# if they are the same (i.e. both Otter) create a temporary cs version MAPPING
my $sql_meta_insert = qq{
	insert ignore into meta (meta_key, meta_value) values
	('assembly.mapping', 'chromosome:MAPPING|contig'),
	('assembly.mapping', 'chromosome:MAPPING|contig|clone'),
	('assembly.mapping', 'chromosome:Otter#chromosome:MAPPING')};

my $sql_meta_delete = qq{
	delete from meta where meta_value like 'chromosome%MAPPING%'};

my $sql_cs_insert = qq{
	insert ignore into coord_system (coord_system_id, name, version, rank, attrib) values (100, 'chromosome', 'MAPPING', 100, '')};
my $sql_cs_delete = qq{
	delete from coord_system where version = 'MAPPING'};

my $sql_mc_insert = qq{
	insert ignore into meta_coord (table_name, coord_system_id, max_length)
	values ('exon', 100, 1), ('gene', 100, 1), ('simple_feature', 100, 1), ('transcript', 100, 1)};
my $sql_mc_delete = qq{
	delete from meta_coord where coord_system_id = 100};

my $sql_sr_update = qq{
	update seq_region s, coord_system cs
	set s.coord_system_id = cs.coord_system_id
	where s.name = ?
	and cs.name = 'chromosome'
	and cs.version = ?};

# 1 and 2 below must be specifically set in order to
# fetch the features on a particular seq_region and disable
# the use of the mapping mechanism (used to test duplicates)

# 1. sql queries to enable gene/simple_feature build
my $sql_meta_bl_insert = qq{
	insert ignore into meta (meta_key, meta_value) values
	('simple_featurebuild.level', '1'),
	('genebuild.level', '1')};
my $sql_meta_bl_delete =
qq{delete from meta where meta_key in ('simple_featurebuild.level','genebuild.level')};

# 2. attrib_type_id 6 is a Top Level Non-Redundant Sequence Region
my $sql_attrib_type_select =
qq{SELECT COUNT(*) FROM seq_region_attrib WHERE seq_region_id = ? AND attrib_type_id = 6;};
my $sql_attrib_type_insert =
qq{insert into seq_region_attrib (seq_region_id, attrib_type_id, value) values (?, 6, 1)};
my $sql_attrib_type_delete =
qq{delete from seq_region_attrib where seq_region_id = ? and attrib_type_id = 6};

my $sth_attrib_type_select = $dbh->prepare($sql_attrib_type_select);
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

my $chr_list = $sliceAd->fetch_all('chromosome');

my @from_chrs = sort {lc($a->seq_region_name) cmp lc($b->seq_region_name)}  grep ($_->seq_region_name =~ /$from_assembly/, @$chr_list);
my @to_chrs = sort {lc($a->seq_region_name) cmp lc($b->seq_region_name)} grep ($_->seq_region_name =~ /$to_assembly/, @$chr_list);


if( scalar(@from_chrs) != scalar(@to_chrs) ) {
	throw("Chromosome lists do not match by length:\n[".
	join(" ",map($_->seq_region_name,@from_chrs))."]\n[".
	join(" ",map($_->seq_region_name,@to_chrs))."]\n");
}

# Check that the chromosome names match
for my $i ( 0 .. scalar(@from_chrs) - 1 ) {
	my $R_sr_name = $from_chrs[$i]->seq_region_name;
	my $A_sr_name = $to_chrs[$i]->seq_region_name;
	my ($R_chr) = $R_sr_name =~ /(.*)_$from_assembly/;
	my ($A_chr) = $A_sr_name =~ /(.*)_$to_assembly/;
	throw("chromosome names don't match $R_chr != $A_chr\n[".
	join(" , ",map($_->seq_region_name,@from_chrs))."]\n[".
	join(" , ",map($_->seq_region_name,@to_chrs))."]\n") unless $R_chr eq $A_chr;

	$support->log("$R_sr_name	=>	$A_sr_name\n");
}

my @R_chr_list = map $_->seq_region_name , @from_chrs;
my @A_chr_list = map $_->seq_region_name , @to_chrs;

# Lock only the old assembly as the new one is supposed to be hidden
my ($cb,$author_obj);

my $locked_slices = [];
eval {
	$cb = Bio::Vega::ContigLockBroker->new(-hostname => hostname);
	$author_obj = Bio::Vega::Author->new(-name => $author, -email => $email);
	foreach my $slice (@from_chrs) {
		$support->log_verbose("Locking ".$slice->seq_region_name."\n");
		$cb->lock_clones_by_slice($slice,$author_obj,$vega_db);
		push @$locked_slices, $slice;
	}
};
if($@){
	# remove the stored locks on the locked slices
	foreach my $slice (@$locked_slices) {
		$cb->remove_by_slice($slice,$author_obj,$vega_db);
	}
	throw("Problem locking assembly $from_assembly chromosomes with author name $author\n$@\n");
}

# set the top level attribute on old and new assembly
my $is_top_level ;
foreach my $slice (@from_chrs, @to_chrs) {
	my $srid = $slice->get_seq_region_id;
	$sth_attrib_type_select->execute($srid);
	$is_top_level->{$srid} = $sth_attrib_type_select->fetchrow_array();
	$sth_attrib_type_insert->execute($srid) unless $is_top_level->{$srid};
}

# change the coord_system of the new assembly if not different
my $sth_cs = $dbh->prepare($sql_sr_update);
if($cs_change) {
	foreach my $A_chr (@A_chr_list) {
		$sth_cs->execute( $A_chr, $altassembly );
	}
}

# Loop over the old assembly chromosomes and transfer the genes/simple_features to the new assembly
$dbh->begin_work;
my $R_chr;

eval {
	SET: for my $i ( 0 .. scalar(@R_chr_list) - 1 ) {
		$R_chr = $R_chr_list[$i];

		# for stats
		my $total_genes      = 0;
		my $transfered_genes = 0;
		my $created_genes	 = 0;
		my $total_sf         = 0;
		my $transfered_sf    = 0;
		my $skipped_sf       = 0;
		my $skipped_g       = 0;

		my $ref_sl =
	  		$sliceAd->fetch_by_region( 'chromosome', $R_chr, undef, undef, undef,$assembly );

		$support->log_verbose("Annotation transfer $R_chr => $to_assembly assembly...\n");


		my @genes;
		@genes       = @{ $ref_sl->get_all_Genes };
		$total_genes = scalar @genes;
		@genes = sort { $a->start <=> $b->start } @genes;

		# transfer simple features (PolyA_site/_signal)
		my @proj_feat;
		my @simple_features = @{ $sfeat_Ad->fetch_all_by_Slice($ref_sl) };
		$total_sf = scalar @simple_features;
		foreach my $f (@simple_features) {
			my $tf = $f->transform('chromosome',$altassembly);
			if ( defined $tf ) {
				$tf->dbID(undef);
				$tf->adaptor(undef);
				my $existing_sf = $sfeat_Ad->fetch_all_by_Slice( $tf->feature_Slice, $tf->analysis->logic_name );
				if ( ! @$existing_sf ) {
					push @proj_feat, $tf;
					$transfered_sf++;
				}
				else {
					$support->log_verbose(
						sprintf(
							"SKIP %s %d %d %d, already saved on %s\n",
							$tf->display_label, $tf->start,
							$tf->end,           $tf->strand,
							$tf->seq_region_name
						)
					);
					$skipped_sf++;
				}

			}
			else {
				$support->log_verbose(
					sprintf(
						" %s %d %d %d cannot be transfered on $to_assembly assembly\n",
						$f->display_label, $f->start, $f->end, $f->strand
					)
				);
				$skipped_sf++;

			}

		}

	  GENE: foreach my $g (@genes) {
			my $tg = Bio::Vega::Gene->new;
			$tg->analysis( $g->analysis );
			$tg->biotype( $g->biotype );
			$tg->status( $g->status );
			$tg->source( $g->source );
			$tg->stable_id( $g->stable_id );

			$tg->gene_author( $g->gene_author );
			$tg->description( $g->description );

			$support->log_verbose(
				sprintf(
					"Transferring gene %s %s %d %d\n",
					$g->stable_id, $g->biotype, $g->start, $g->end
				)
			);

			my @proj_trans;
			my $transcript = $g->get_all_Transcripts;

		  TRANSCRIPT: foreach my $t ( @{ $transcript } ) {
				my $tt = &transform($t);
				if ( defined $tt ) {
					push @proj_trans, $tt;
					$support->log_verbose(
						sprintf(
							"\t%s %s %d %d transferred successfully on %s:\n",
							$t->stable_id, $t->biotype, $t->start, $t->end,
							$tt->seq_region_name
						)
					);
					$tt->status( $t->status );
					&remove_all_db_ids($tt);
					&log_compare_transcripts( $t, $tt );
				} else {
					$support->log_verbose(
						sprintf(
							"\t%s %s %d %d cannot be transfered on $to_assembly assembly\n",
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
					scalar( @{ $transcript } )
				)
			);

			my $missing_transcript = scalar(@{ $transcript }) - scalar(@proj_trans);

			if (  $missing_transcript == 0 ) {
			   # essential for loading attribute list so that get_all_Attributes
			   # won't return empty list
				$tg->add_Attributes( @{ $g->get_all_Attributes() } );

				my $genes = &transcripts2genes($tg,\@proj_trans);

				foreach my $gene (@$genes) {
					&write_gene($gene);
					my $existing_gene;
					eval {
						$existing_gene = $geneAd->fetch_all_by_Slice( $gene->feature_Slice,
									$gene->analysis->logic_name );
					};
					if($@) {
						$support->log_verbose($@);
						next GENE;
					}
					my $exist = 0;
					my $tg_key = join(":",$gene->start,$gene->end,$gene->biotype,$gene->status,$gene->source);
					foreach (@$existing_gene) {
						my $existing_gkey = join(":",$_->start,$_->end,$_->biotype,$_->status,$_->source);
						$exist = 1 if $tg_key eq $existing_gkey;
					}

					if($exist) {
						$support->log_verbose(
							sprintf(
								"SKIP GENE %s %s (%s:%d-%d => %s:%d-%d) already transfered\n",
								$gene->stable_id, $gene->biotype, $R_chr,
								$g->start,      $g->end,      $gene->seq_region_name,
								$gene->start,     $gene->end
							));

						$support->log_verbose("Transfered gene: ".join("\t",
						$gene->analysis->logic_name,
						$gene->biotype,
						$gene->status,
						$gene->source,
						$gene->stable_id,
						$gene->gene_author->name,
						$gene->description,
						scalar(@{$gene->get_all_Transcripts}))."\n");
						foreach my $eg (@$existing_gene){
							$support->log_verbose("Existing gene: ".join("\t",
							$eg->analysis->logic_name,
							$eg->biotype,
							$eg->status,
							$eg->source,
							$eg->stable_id,
							$eg->gene_author->name,
							$eg->description,
							scalar(@{$eg->get_all_Transcripts}))."\n");
						}

						$skipped_g++;
						next GENE;
					}
				}

				$created_genes += (scalar(@$genes) -1);

				foreach my $gene (@$genes) {
					if ( $geneAd->store($gene) ) {
						$transfered_genes++;
						$support->log_verbose(
							sprintf(
							"GENE %s %s successfully TRANSFERED (%s:%d-%d => %s:%d-%d)\n",
								$gene->stable_id, $gene->biotype, $R_chr,
								$g->start,      $g->end,      $gene->seq_region_name,
								$gene->start,     $gene->end
							)
						);
					} else {
						throw(
							sprintf(
								"GENE %s %s cannot be saved (%s:%d-%d => %s:%d-%d)",
								$gene->stable_id, $gene->biotype, $R_chr,
								$g->start,      $g->end,      $gene->seq_region_name,
								$gene->start,     $gene->end
							)
						);
					}
				}
				if(scalar(@$genes) -1) {
					# print info about splitted gene for the annotators
					$support->log_verbose(sprintf("WARNING: Check Gene %s, it has been splitted into %s\n",$g->stable_id,join(",",map($_->stable_id,@$genes))));
				}
			} else {
				$support->log_verbose(
					sprintf(
								"SKIP GENE %s %s (%s:%d-%d => %s assembly) with %d missing transcripts\n",
								$tg->stable_id, $tg->biotype, $R_chr,
								$g->start,      $g->end,      $to_assembly,
								$missing_transcript
							)
				);
			}
		}

		# save polyA features if gene transfered
		$sfeat_Ad->store(@proj_feat) if ( @proj_feat && $transfered_genes );

		$support->log_verbose(
sprintf(
"INFO: Annotation transfer %s:%s => %s
INFO: transfered Gene: %d/%d
INFO: created Gene:	%d
INFO: skipped Gene: %d/%d
INFO: transfered PolyA features: %d/%d
INFO: skipped PolyA features: %d/%d\n",
				$R_chr,            $assembly,
				$to_assembly,
				$transfered_genes, ($total_genes+$created_genes),
				$created_genes,
				$skipped_g, ($total_genes+$created_genes),
				$transfered_sf,    $total_sf,
				$skipped_sf,       $total_sf
			)
		);
		$support->log_stamped("INFO: Done with $R_chr\n");

	}

	$write_db ? $dbh->commit : $dbh->rollback;
};

if ($@) {
	$dbh->rollback;
	$support->log_verbose(
		    "UNABLE TO TRANSFER ANNOTATION FROM $from_assembly assembly ($R_chr:$assembly) to $to_assembly assembly\n[$@]\n"
	);
}


# set the alt-assemblies to the previous coord_system and
# remove the temporary 'MAPPING' coordinate system
if($cs_change) {
	foreach my $A_chr (@A_chr_list) {
		$sth_cs->execute( $A_chr, $assembly );
	}
	$dbh->do($sql_meta_delete);
	$dbh->do($sql_cs_delete);
	$dbh->do($sql_mc_delete);
	$dbh->do($sql_meta_bl_delete);
}

# remove the assemblies locks
eval {
	foreach my $slice (@$locked_slices) {
		$support->log_verbose("Unlocking ".$slice->seq_region_name."\n");
		$cb->remove_by_slice($slice,$author_obj,$vega_db);
	}
};
if($@){
	warning("Cannot remove locks from assembly $from_assembly chromosomes with author name $author\n$@\n");
}

# remove the top level attribute on old and new assembly
# leave existing attributes
foreach my $slice (@from_chrs, @to_chrs) {
	my $srid = $slice->get_seq_region_id;
	$sth_attrib_type_delete->execute($srid) unless $is_top_level->{$srid};
}

# finish logfile
$support->finish_log;
### end main

sub transform {
	my $self = shift;
	my @coords = ('chromosome',$altassembly);
	my $new_transcript;
	if ( !defined $new_transcript ) {
		my @segments = $self->project(@coords);

		# if it projects, maybe the exons transform well?
		# lazy load them here
		return undef if ( !@segments );
		$self->get_all_Exons();
	}
	if ( exists $self->{'_trans_exon_array'} ) {
		my @new_exons;
		my ( $low_start, $hi_end, $slice );

		# we want to check whether the transform preserved 5prime 3prime
		# ordering. This assumes 5->3 order. No complaints on transsplicing.
		my ( $last_new_start, $last_old_strand, $last_new_strand, $start_exon,
			 $end_exon, $last_seq_region_name );
		my $first        = 1;
		my $ignore_order = 1;
		my $order_broken = 0;
		for my $old_exon ( @{ $self->{'_trans_exon_array'} } ) {
			my $new_exon = $old_exon->transform(@coords);
			return undef if ( !defined $new_exon );
			if ( !defined $new_transcript ) {
				if ( !$first ) {
					if ( $old_exon->strand() != $last_old_strand ) {

						# transsplicing, ignore ordering
						$ignore_order = 1;
					}
					if ( $new_exon->slice()->seq_region_name() ne
						 $last_seq_region_name )
					{
						return undef;
					}
					if (     $last_new_strand == 1
						 and $new_exon->start() < $last_new_start )
					{
						$order_broken = 1;
					}
					if (     $last_new_strand == -1
						 and $new_exon->start() > $last_new_start )
					{
						$order_broken = 1;
					}
					if ( $new_exon->start() < $low_start ) {
						$low_start = $new_exon->start();
					}
					if ( $new_exon->end() > $hi_end ) {
						$hi_end = $new_exon->end();
					}
				}
				else {
					$first     = 0;
					$low_start = $new_exon->start();
					$hi_end    = $new_exon->end();
				}
				$last_seq_region_name = $new_exon->slice()->seq_region_name();
				$last_old_strand      = $old_exon->strand();
				$last_new_start       = $new_exon->start();
				$last_new_strand      = $new_exon->strand();
			}
			if ( defined $self->{'translation'} ) {
				if ( $self->translation()->start_Exon() == $old_exon ) {
					$start_exon = $new_exon;
				}
				if ( $self->translation()->end_Exon() == $old_exon ) {
					$end_exon = $new_exon;
				}
			}
			push( @new_exons, $new_exon );
		}
		if ( $order_broken && !$ignore_order ) {
			warning( "Order of exons broken in transform of " . $self->dbID() );
			return undef;
		}
		if ( !defined $new_transcript ) {
			%$new_transcript = %$self;
			bless $new_transcript, ref($self);
			$new_transcript->start($low_start);
			$new_transcript->end($hi_end);
			$new_transcript->slice( $new_exons[0]->slice() );
			$new_transcript->strand( $new_exons[0]->strand() );
		}
		$new_transcript->{'_trans_exon_array'} = \@new_exons;

		# should be ok to do inside exon array loop
		# translations only exist together with the exons ...
		if ( defined $self->{'translation'} ) {
			my $new_translation;
			%$new_translation = %{ $self->{'translation'} };
			bless $new_translation, ref( $self->{'translation'} );
			$new_transcript->{'translation'} = $new_translation;
			$new_translation->start_Exon($start_exon);
			$new_translation->end_Exon($end_exon);
		}
	}
	if ( exists $self->{'_supporting_evidence'} ) {
		my @new_features;
		for my $old_feature ( @{ $self->{'_supporting_evidence'} } ) {
			my $new_feature = $old_feature->transform(@coords);
			if ( defined $new_feature ) {
				push @new_features, $new_feature;
			}
		}
		$new_transcript->{'_supporting_evidence'} = \@new_features;
	}

	# flush cached internal values that depend on the exon coords
	$new_transcript->{'transcript_mapper'}   = undef;
	$new_transcript->{'coding_region_start'} = undef;
	$new_transcript->{'coding_region_end'}   = undef;
	$new_transcript->{'cdna_coding_start'}   = undef;
	$new_transcript->{'cdna_coding_end'}     = undef;
	return $new_transcript;
}


sub transcripts2genes {
	my ($tg,$trans) = @_;
	my @genes;
	# group the transcripts by seq_region_name and minus/plus strand
	my $hash_trans;
	for (@$trans) {
		$hash_trans->{$_->seq_region_name} ||= {};
		&parsenpush_trans($hash_trans->{$_->seq_region_name},$_);
	}

	# Use a hash to split genes that have:
	# 1. transcripts on different chromosomes
	# 2. transcripts on different strands
	# 3. transcript with exons on different strands


	# sort the sets by translatable transcripts number, total exons length and transcript number
	my @gene_transcripts;
	map (push(@gene_transcripts, values %$_), values(%$hash_trans));
	@gene_transcripts = sort sort_transcripts @gene_transcripts;
	my $number = scalar @gene_transcripts;

	my $gt = shift @gene_transcripts;

	my $set = 2;
	# loop over the remaining transcript sets and create new gene without stable id
	foreach my $nt (@gene_transcripts) {
		my $ng = Bio::Vega::Gene->new;
		$ng->analysis( $tg->analysis );
		$ng->biotype( $tg->biotype );
		$ng->status( $tg->status );
		$ng->source( $tg->source );
		$ng->gene_author( $tg->gene_author );
		$ng->description( $tg->description );
		# add gene attributes and change gene name
		foreach (@{ $tg->get_all_Attributes() }) {
			if($_->code eq 'name') {
				my $attrib = new Bio::EnsEMBL::Attribute->new
			       (-CODE => $_->code,
			        -NAME => $_->name,
			        -DESCRIPTION => $_->description,
			        -VALUE => $_->value."_$set");
				$ng->add_Attributes($attrib);
			} else {
				$ng->add_Attributes($_);
			}

		}
		my $remark = "Gene $set/$number (".$tg->stable_id.") automatically splitted by the annotation transfer script";
		my $attribute = Bio::EnsEMBL::Attribute->new
	       (-CODE => 'hidden_remark',
	        -NAME => 'Hidden Remark',
	        -VALUE => $remark);
		$ng->add_Attributes($attribute);

		map( $ng->add_Transcript($_), @$nt );
		push @genes,$ng;
		$set++;
	}

	# create the main gene here with same stable id
	map( $tg->add_Transcript($_), @$gt );
	my $remark = "Gene 1/$number (".$tg->stable_id.") automatically splitted by the annotation transfer script";
	my $attribute = Bio::EnsEMBL::Attribute->new
       (-CODE => 'hidden_remark',
        -NAME => 'Hidden Remark',
        -VALUE => $remark);
	$tg->add_Attributes($attribute) if(@gene_transcripts);
	my ($name_attrib) = @{ $tg->get_all_Attributes('name') };
	# change main gene name
	$name_attrib->value($name_attrib->value."_1") if(@gene_transcripts);;
	push @genes,$tg;

	return \@genes;
}

sub sort_transcripts {

	&get_translatable_trans($b) <=> &get_translatable_trans($a)
								||
		  &get_exons_length($b) <=> &get_exons_length($a)
								||
					scalar(@$b) <=> scalar(@$a);
}

sub sort_exons {

		  &get_exons_length($b) <=> &get_exons_length($a)
								||
					scalar(@$b) <=> scalar(@$a);
}

sub get_translatable_trans {
	my ($trans_arr) = @_;
	my $i;
	foreach (@$trans_arr) { $i++ if($_->translate); }

	return $i;
}

sub get_exons_length {
	my ($object_arr) = @_;
	my $l;
	foreach (@$object_arr) { $l += $_->length; }

	return $l;
}

sub parsenpush_trans{
	my ($hash,$trans) = @_;
	my $exon_hash = {};
	foreach my $exon (@{$trans->get_all_Exons}) {
		$exon_hash->{$exon->strand} ||= [];
		push @{$exon_hash->{$exon->strand}}, $exon;
	}

	if(scalar(keys %{$exon_hash}) == 1) {
		$hash->{$trans->strand} ||= [];
		push @{$hash->{$trans->strand}}, $trans;
	} else {
		# split transcript with trans splice exons
		my $translation_stable_id;
		if ( defined $trans->translation ) {
			$translation_stable_id = $trans->translation->stable_id;
			my $remark = "Lost translation $translation_stable_id after transfer of transcript ".$trans->stable_id;
			my $attribute = Bio::EnsEMBL::Attribute->new
		       (-CODE => 'hidden_remark',
		        -NAME => 'Hidden Remark',
		        -VALUE => $remark);
		    $trans->add_Attributes($attribute);
		}
		my @transcript_exons = sort sort_exons values(%$exon_hash);

		# create two transcripts, 1 on each strand
		my @attribs = @{$trans->get_all_Attributes()};

		# 2nd transcript
		my $new_trans;
		%$new_trans = %$trans;
		bless $new_trans, ref($trans);
		$new_trans->{'attributes'} = undef;
		$new_trans->{'_trans_exon_array'} = pop @transcript_exons;
		$new_trans->{'stable_id'} = undef;
		$new_trans->{'translation'} = undef;
		# Add transcript attribs and change transcript name
		foreach (@attribs) {
			if($_->code eq 'name') {
				my $attrib = new Bio::EnsEMBL::Attribute->new
			       (-CODE => $_->code,
			        -NAME => $_->name,
			        -DESCRIPTION => $_->description,
			        -VALUE => $_->value.'_2');
				$new_trans->add_Attributes($attrib);
			} else {
				$new_trans->add_Attributes($_);
			}
		}

		# 1st transcript
		$trans->{'_trans_exon_array'} = pop @transcript_exons;
		$trans->{'translation'} = undef;
		# Change transcript name
		foreach (@attribs) { $_->value($_->value.'_1') if($_->code eq 'name'); }

		$trans->recalculate_coordinates();
		$new_trans->recalculate_coordinates();

		foreach ($trans,$new_trans) {
			$hash->{$_->strand} ||= [];
			push @{$hash->{$_->strand}}, $_;
		}
	}
}


sub remove_all_db_ids {
	my $t = shift;

	$t->dbID(undef);
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
#	foreach my $a (@{ $g->get_all_Attributes })	{
#		printf "code %s\tname %s\tvalue %s\n",$a->code,$a->name,$a->value;
#	}

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
	}

	foreach my $t ( @{ $g->get_all_Transcripts } ) {
		printf( "\t$seq_id\ttranscript\t%s\t%s\t%s\t%d\t%d\t%d\n",
			$t->stable_id, $t->biotype, $t->status, $t->start, $t->end, $t->strand );

#		foreach my $a (@{ $t->get_all_Attributes })	{
#			printf "code %s\tname %s\tvalue %s\n",$a->code,$a->name,$a->value;
#		}

		foreach my $e ( @{ $t->get_all_Exons } ) { #sort { $a->start cmp $b->start }
			printf( "\t$seq_id\texon_transcript\t%s\t%d\t%d\t%d\t%d\t%d\n",
				$e->stable_id, $e->start, $e->end, $e->strand, $e->phase,
				$e->end_phase );
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

