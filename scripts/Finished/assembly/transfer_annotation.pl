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
  --write							  to write the changes
  --author							  author login used to lock the assemblies (must be set to write the changes)
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
use Bio::Vega::DBSQL::DBAdaptor;
use Bio::Vega::ContigLockBroker;
use Bio::Vega::Author;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
	'assembly=s',         'altassembly=s',
	'chromosomes|chr=s@', 'altchromosomes|altchr=s@', 'write!',
	'reverse!', 'haplotype!', 'ref_start=i', 'ref_end=i', 'alt_start=i',
	'alt_end=i', 'author=s', 'email=s'
);
$support->allowed_params( $support->get_common_params, 'assembly',
	'altassembly', 'chromosomes', 'altchromosomes', );

if ( $support->param('help') or $support->error ) {
	print STDOUT $support->error if $support->error;
	pod2usage(1);
}

$support->param('verbose', 1);
$support->param('interactive', 0);

$support->comma_to_list( 'chromosomes', 'altchromosomes' );

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params( 'assembly', 'altassembly' );

# database connection
my $vega_db = $support->get_database('ensembl');
bless $vega_db, "Bio::Vega::DBSQL::DBAdaptor";
my $dbh = $vega_db->dbc->db_handle;

my $assembly    = $support->param('assembly');
my $altassembly = $support->param('altassembly');
my $write 		= $support->param('write');
my $author 		= $support->param('author');
my $email 		= $support->param('email') || $author;
my $haplo 		= $support->param('haplotype');

my $R_start = $support->param('ref_start') || undef;
my $R_end = $support->param('ref_end') || undef;
my $A_start = $support->param('alt_start') || undef;
my $A_end = $support->param('alt_end') || undef;

throw("must set author name to lock the assemblies if you want to write the changes") if (!$author && $write);


my $geneAd   = $vega_db->get_GeneAdaptor;
my $sliceAd  = $vega_db->get_SliceAdaptor;
my $sfeat_Ad = $vega_db->get_SimpleFeatureAdaptor();

# make sure that the coordinate system versions are different
# if they are the same (i.e. both Otter) create a temporary cs version MAPPING
my $sql_meta_insert = qq{
	insert ignore into meta (meta_key, meta_value) values
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
my $sql = 0;

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
	my $created_genes	 = 0;
	my $transfered_sf    = 0;
	my $skipped_sf       = 0;
	my $skipped_g        = 0;

	my $sth_cs = $dbh->prepare($sql_sr_update);
	$sth_cs->execute( $A_chr, $altassembly ) unless ( !$cs_change );

	my $ref_sl =
  		$sliceAd->fetch_by_region( 'chromosome', $R_chr, $R_start, $R_end, undef,$assembly );
	my $alt_sl =
  		$sliceAd->fetch_by_region( 'chromosome', $A_chr, $A_start, $A_end, undef,$altassembly );

  	$vega_db->get_AssemblyMapperAdaptor()->delete_cache();

	my $alt_seq_region_id = $alt_sl->get_seq_region_id;
	my $ref_seq_region_id = $ref_sl->get_seq_region_id;

	$sth_attrib_type_select->execute($alt_seq_region_id);
	my $alt_attrib_set = $sth_attrib_type_select->fetchrow_array;
	$sth_attrib_type_insert->execute($alt_seq_region_id) unless $alt_attrib_set;

	$sth_attrib_type_select->execute($ref_seq_region_id);
	my $ref_attrib_set = $sth_attrib_type_select->fetchrow_array;
	$sth_attrib_type_insert->execute($ref_seq_region_id) unless $ref_attrib_set;

	# Lock the reference and alternative assemblies
	my ($cb,$author_obj);
	my $refslice_locked = 0;
	eval {
		$cb = Bio::Vega::ContigLockBroker->new(-hostname => hostname);
		$author_obj = Bio::Vega::Author->new(-name => $author, -email => $email);
		$support->log_verbose("Locking $R_chr\n");
		$cb->lock_clones_by_slice($ref_sl,$author_obj,$vega_db);
		$refslice_locked = 1;
	};
	if($@){
		warning("Cannot lock assemblies $R_chr and $A_chr with author name $author\n$@\n");
		# remove the stored locks on the reference slice in case alt slice locking fails
		$cb->remove_by_slice($ref_sl,$author_obj,$vega_db) if $refslice_locked;

			$sth_cs->execute( $A_chr, $support->param('altassembly') )
		  unless ( !$cs_change );
		$sth_attrib_type_delete->execute($alt_seq_region_id) unless $alt_attrib_set;
		$sth_attrib_type_delete->execute($ref_seq_region_id) unless $ref_attrib_set;

			next SET;
	}

	$dbh->begin_work;

	eval {
		$support->log_verbose("Annotation transfer $R_chr => $A_chr...\n");


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
				my $existing_sf = $sfeat_Ad->fetch_all_by_Slice( $tf->feature_Slice, $tf->analysis->logic_name );
				if ( ! @$existing_sf ) {
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
				$support->log_verbose(
					sprintf(
						" %s %d %d %d cannot be transfered on $A_chr\n",
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
			my $transcript = $g->get_all_Transcripts;

		  TRANSCRIPT: foreach my $t ( @{ $transcript } ) {
				my $tt = $t->transfer($alt_sl);
				if ( defined $tt ) {
					#my $status = 'CLEAN_TRANSFER';
        			#add_trans_remark($R_chr, $status, $tt);
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
				} else {
					if(&transcript_is_missed($t,$alt_sl)) {
						$support->log_verbose(
							sprintf(
								"\t%s %s %d %d cannot be transfered on $A_chr\n",
								$t->stable_id, $t->biotype, $t->start, $t->end
							)
						);
						&add_hidden_remark(	$t->slice->seq_region_name,
						    				sprintf("transcript %s cannot be transfered",$t->stable_id),
						    				$tg );

			        } else {
						if ($t->translation and $haplo) {
							my $cds_t = &get_coding_part_of_transcript($t);
				            my $tt = $cds_t->transfer($alt_sl);
				            if (defined $tt) {
				            	$tt->status( $t->status );
				              	#my $status = 'CLEAN_CDS_ONLY_TRANSFER';
				              	#add_trans_remark($R_chr, $status, $tt);
				              	push @proj_trans, $tt;
								$support->log_verbose(
					              	sprintf("\t%s %s %d %d CDS transferred cleanly:\n",
					                	$t->stable_id,
					                    $t->biotype,
					                    $t->start,
					                    $t->end)
								);
				            	&log_compare_transcripts($cds_t, $tt);
				            	&remove_all_db_ids($tt);
				            } else {
				            	my $new_t = &project_transcript_the_hard_way($g, $t, $alt_sl);
				              	$new_t->status( $t->status );
				              	#my $status = 'COMPLEX_CODING_WITHOUT_CDS';
				              	#add_trans_remark($R_chr, $status, $new_t);
				              	push @proj_trans, $new_t;
								$support->log_verbose(
				              		sprintf("\t%s %s %d %d complex transfer of coding without CDS\n",
				                    	$t->stable_id,
				                     	$t->biotype,
				                     	$t->start,
				                     	$t->end)
				                );
				                &remove_all_db_ids($new_t);
				              	&log_summarise_projection_by_exon($t, $alt_sl);
				            }
						} else {
							my $new_t = &project_transcript_the_hard_way($g, $t, $alt_sl );
							$new_t->status( $t->status );
				            #my $status = 'COMPLEX_NONCODING';
				            #add_trans_remark($ref_set, $status, $new_t);
							push @proj_trans, $new_t;
							$support->log_verbose(
			              		sprintf("\t%s %s %d %d complex transfer of non-coding\n",
			                     $t->stable_id,
			                     $t->biotype,
			                     $t->start,
			                     $t->end)
							);
							&remove_all_db_ids($new_t);
			              	&log_summarise_projection_by_exon($t, $alt_sl);
						}
			        }
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

			if (  scalar(@proj_trans) ) {
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
						$support->log_verbose(sprintf("WARNING: Check Gene %s with %d missing transcripts\n",$gene->stable_id,$missing_transcript)) if $missing_transcript;
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
								"SKIP GENE %s %s (%s:%d-%d => %s:%d-%d) with %d missing transcripts\n",
								$tg->stable_id, $tg->biotype, $R_chr,
								$g->start,      $g->end,      $A_chr,
								$tg->start,     $tg->end,
								$missing_transcript
							)
				);
			}
		}

		# save polyA features if gene transfered
		$sfeat_Ad->store(@proj_feat) if ( @proj_feat && $transfered_genes );

		$write ? $dbh->commit : $dbh->rollback;
	};

	if ($@) {
		$dbh->rollback;
		$support->log_verbose(
			    "UNABLE TO TRANSFER ANNOTATION FROM $R_chr:$assembly to $A_chr:"
			  . $support->param('altassembly')
			  . "\n[$@]\n" );
	} else {
		# print annotation transfer stats
		$support->log_verbose(
			sprintf(
"INFO: Annotation transfer %s:%s => %s/%s
INFO: transfered Gene: %d/%d
INFO: created Gene:	%d
INFO: skipped Gene: %d/%d
INFO: transfered PolyA features: %d/%d
INFO: skipped PolyA features: %d/%d\n",
				$R_chr,            $assembly,
				$A_chr,            $support->param('altassembly'),
				$transfered_genes, ($total_genes+$created_genes),
				$created_genes,
				$skipped_g, ($total_genes+$created_genes),
				$transfered_sf,    $total_sf,
				$skipped_sf,       $total_sf
			)
		);
	}

	# remove the assemblies locks
	eval {
		$support->log_verbose("Removing $R_chr Locks\n");
		$cb->remove_by_slice($ref_sl,$author_obj,$vega_db);
	};
	if($@){
		warning("Cannot remove locks from assemblies $R_chr and $A_chr with author name $author\n$@\n");
	}

	$sth_cs->execute( $A_chr, $support->param('altassembly') )
	  unless ( !$cs_change );
	$sth_attrib_type_delete->execute($alt_seq_region_id) unless $alt_attrib_set;
	$sth_attrib_type_delete->execute($ref_seq_region_id) unless $ref_attrib_set;
}

# remove temporary 'MAPPING' coordinate system from meta and coord_system table
if ($cs_change) {
	$dbh->do($sql_meta_delete);
	$dbh->do($sql_cs_delete);
	$dbh->do($sql_mc_delete);
	$dbh->do($sql_meta_bl_delete);
}

$support->log_stamped("\nDone.\n");

# finish logfile
$support->finish_log;

sub print_sql {
	my ($g,$tg, $ref_sid, $alt_sid, $R_chr, $A_chr) = @_;

	my $update_filename = "update_${R_chr}_${A_chr}.sql";
	my $update_back_filename = "update_${R_chr}_${A_chr}_back.sql";

	open(OUT, ">>$update_filename") or die "cannot create file $update_filename\n";
	open(BACK, ">>$update_back_filename") or die "cannot create file $update_back_filename\n";

	# Gene SQL updates
	print OUT sprintf qq{
 UPDATE gene SET seq_region_id = %d,seq_region_start = %d,seq_region_end = %d,seq_region_strand = %d
 WHERE gene_id = %d AND seq_region_id = %d;
    },$alt_sid,$tg->start,$tg->end,$tg->strand,$g->dbID,$ref_sid;
    print BACK sprintf qq{
 UPDATE gene SET seq_region_id = %d,seq_region_start = %d,seq_region_end = %d,seq_region_strand = %d
 WHERE gene_id = %d AND seq_region_id = %d;
    },$ref_sid,$g->start,$g->end,$g->strand,$g->dbID,$alt_sid;

    # Transcrit SQL updates
    my ($trans,$ttrans);
	map($trans->{$_->stable_id} = $_ , @{ $g->get_all_Transcripts });
	map($ttrans->{$_->stable_id} = $_ , @{ $tg->get_all_Transcripts });
	foreach my $stable_id (sort keys %$trans) {
		my $t = $trans->{$stable_id};
		my $tt = $ttrans->{$stable_id};
		print OUT sprintf qq{
   UPDATE transcript SET seq_region_id = %d,seq_region_start = %d,seq_region_end = %d,seq_region_strand = %d
   WHERE transcript_id = %d AND seq_region_id = %d;
        },$alt_sid,$tt->start,$tt->end,$tt->strand,$t->dbID,$ref_sid;
        print BACK sprintf qq{
   UPDATE transcript SET seq_region_id = %d,seq_region_start = %d,seq_region_end = %d,seq_region_strand = %d
   WHERE transcript_id = %d AND seq_region_id = %d;
        },$ref_sid,$t->start,$t->end,$t->strand,$t->dbID,$alt_sid;

        # Exon SQL updates
		my ($exon,$texon);
		map($exon->{$_->stable_id} = $_ , @{ $t->get_all_Exons() });
		map($texon->{$_->stable_id} = $_ , @{ $tt->get_all_Exons() });
		foreach my $sid (sort keys %$exon) {
			my $e = $exon->{$sid};
			my $te = $texon->{$sid};
			print OUT sprintf qq{
     UPDATE exon SET seq_region_id = %d,seq_region_start = %d,seq_region_end = %d,seq_region_strand = %d, phase = %d, end_phase = %d
     WHERE exon_id = %d AND seq_region_id = %d;
            },$alt_sid,$te->start,$te->end,$te->strand,$te->phase,$te->end_phase,$e->dbID,$ref_sid;
            print BACK sprintf qq{
     UPDATE exon SET seq_region_id = %d,seq_region_start = %d,seq_region_end = %d,seq_region_strand = %d, phase = %d, end_phase = %d
     WHERE exon_id = %d AND seq_region_id = %d;
            },$ref_sid,$e->start,$e->end,$e->strand,$e->phase,$e->end_phase,$e->dbID,$alt_sid;
		}
	}
}

sub remove_all_db_ids {
	my $t = shift;

	$t->dbID(undef);
	$t->stable_id(undef) if $haplo;
	$t->adaptor(undef);
	if ( $t->translation ) {
		$t->translation->dbID(undef);
		$t->translation->adaptor(undef);
		$t->translation->stable_id(undef) if $haplo;
	}
	foreach my $e ( @{ $t->get_all_Exons } ) {
		$e->dbID(undef);
		$e->adaptor(undef);
		$e->stable_id(undef) if $haplo;
	}

	return $t;
}

sub transcript_is_missed {
	my ( $t, $alt_sl ) = @_;

	foreach my $e ( @{ $t->get_all_Exons } ) {
		my $alt_e = $e->transfer($alt_sl);
		if ( defined $alt_e ) {
			return 0;
		}
	}

	return 1;
}

sub get_coding_part_of_transcript {
	my ($t) = @_;

	my @cds = @{ $t->get_all_translateable_Exons };
	$t->{_trans_exon_array} = [];
	foreach my $e (@cds) { $t->add_Exon($e); }
	my $tr = $t->translation;
	$tr->start_Exon( $cds[0] );
	$tr->start(1);
	$tr->end_Exon( $cds[-1] );
	$tr->end( $tr->end_Exon->length );
	$t->translation($tr);
	$t->stable_id( $t->stable_id );

	return $t;
}

sub project_transcript_the_hard_way {
	my ( $g, $t, $alt_sl ) = @_;

	my $new_trans;
	%$new_trans = %$t;
	bless $new_trans, ref($t);
	$new_trans->{'translation'} = undef;
	if ( defined $t->translation ) {
		&add_hidden_remark(	$t->slice->seq_region_name,
						    sprintf("transcript %s has lost its translation %s",$t->stable_id,$t->translation->stable_id),
						    $new_trans);
	}

	my @new_e;
	my $cs = $alt_sl->coord_system;
	foreach my $e ( @{ $new_trans->get_all_Exons } ) {
		my $te = $e->transfer($alt_sl);

		if ( defined $te ) {
			push @new_e, $te;
		}
		else {
			my @bits = @{ $e->project( $cs->name, $cs->version) };
			if (@bits) {
				# need to make a new exon from the bits
				my ( $coord, $start, $end, $strand );
				foreach my $bit (@bits) {
					$strand = $bit->to_Slice->strand;
					$start = $bit->to_Slice->start;
					$end = $bit->to_Slice->end;
					$coord->{$strand} ||= {};
					$coord->{$strand}->{length} +=  ($end-$start);
					if ( not defined $coord->{$strand}->{start}
					or $start < $coord->{$strand}->{start} )
					{
						$coord->{$strand}->{start} = $start;
					}
					if ( not defined $coord->{$strand}->{end}
						or $end > $coord->{$strand}->{end} )
					{
						$coord->{$strand}->{end} = $end;
					}
				}
				my $flag = 0;
				foreach ( sort { $coord->{$b}->{length} <=> $coord->{$a}->{length} } keys  %$coord ) {
					my $new_e = Bio::Vega::Exon->new;
					$new_e->start($coord->{$_}->{start});
					$new_e->end($coord->{$_}->{end});
					$new_e->strand($_);
					$new_e->phase(-1);
					$new_e->end_phase(-1);
					$new_e->slice($alt_sl);
					$new_e->stable_id( $e->stable_id ) unless $flag;
					push @new_e, $new_e;
					$flag = 1;
				}
				&add_hidden_remark($e->slice->seq_region_name,
									    sprintf("exon %s has been altered",$e->stable_id),
									    $new_trans);
				$support->log_verbose(sprintf("WARNING: Check altered Exon %s in Transcript %s of Gene %s\n",$e->stable_id,$t->stable_id,$g->stable_id));


			} else {
				&add_hidden_remark($e->slice->seq_region_name,
								    sprintf("exon %s %d-%d:%d cannot be transfered",$e->stable_id,$e->start,$e->end,$e->strand),
								    $new_trans);
			}
		}
	}

	$new_trans->{_trans_exon_array} = [];
	foreach my $e (@new_e) {
		$new_trans->add_Exon($e);
	}

	return $new_trans;
}

sub add_hidden_remark {
  my ($ref_set, $status, $obj) = @_;

  my $remark = "Automatic annotation transfer from $ref_set: $status";
  my $t_at = Bio::EnsEMBL::Attribute->new(-VALUE => $remark,
                                          -CODE  => 'hidden_remark',
                                          -NAME  => 'Hidden Remark',
                                         );
  $obj->add_Attributes($t_at);
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


sub log_summarise_projection_by_exon {
	my ( $t, $alt_sl ) = @_;

	my $cs = $alt_sl->coord_system;
	my $exons_total       = 0;
	my $exons_transferred = 0;
	my $exons_missed      = 0;
	my $exons_part_missed = 0;
	my $exons_over_indel  = 0;

	foreach my $e ( @{ $t->get_all_Exons } ) {
		$exons_total++;

		my $te = $e->transfer( $alt_sl );

		if ( defined $te ) {
			$exons_transferred++;
		}
		else {
			my @bits = @{ $e->project( $cs->name, $cs->version ) };

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
			"\t%d exons; %d transferred, %d missed, %d part missed, %d over indel\n",
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
#	foreach my $a (@{ $g->get_all_Attributes })	{
#		printf "code %s\tname %s\tvalue %s\n",$a->code,$a->name,$a->value;
#	}

	foreach my $t ( @{ $g->get_all_Transcripts } ) {
		printf( "\t$seq_id\ttranscript\t%s\t%s\t%s\n",
			$t->stable_id, $t->biotype, $t->status );

#		foreach my $a (@{ $t->get_all_Attributes })	{
#			printf "code %s\tname %s\tvalue %s\n",$a->code,$a->name,$a->value;
#		}

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

