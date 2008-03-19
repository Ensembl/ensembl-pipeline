#!/software/bin/perl

=head1 NAME

align_by_component_identity_zfish.pl - create a whole genome alignment between two closely
related set of zebrafish assemblies, step 1

=head1 SYNOPSIS

align_by_component_identity_zfish.pl [arguments]

Required arguments:

    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --dbname, db_name=NAME              database name NAME
    --from_assembly                     old assembly date
    --to_assembly                       new assembly date

Optional arguments:

    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)
    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)

    -v, --verbose=0|1                   verbose logging (default: true)
    -i, --interactive=0|1               run script interactively (default: false)
    -n, --dry=0|1, --dry_run=0|1        don't write results to database (default: false)
    -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script is part of a series of scripts to create a mapping between two
sets of zebrafish assemblies. It assembles the chromosome coordinate systems of two different
assemblies of a genome by creating a whole genome alignment between the two.

The process handles major rearrangements or components moved from one chromosome to another.

See "Related scripts" below for an overview of the whole process.

This particular script creates a whole genome alignment between two closely
related sets of assemblies. You will need a database containing the reference assembly
and the alternative chromosomes.

The alignment is created in two steps:

    1. Match components with same name and version directly using the sdiff method
       in Algorithm::Diff and create alignment blocks for these regions.

       The result is stored in the assembly table as an assembly between the
       chromosomes of both genome assemblies.

    2. Store non-aligned blocks in a temporary table (tmp_align). They can
       later be aligned using lastz by align_nonident_zfish.pl.

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
    unshift(@INC, "$SERVERROOT/extra");
	unshift(@INC, "$SERVERROOT/bioperl-0.7.2");
    unshift(@INC, "$SERVERROOT/bioperl-1.2.3-patched");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Algorithm::Diff qw(sdiff);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options( 'from_assembly=s', 'to_assembly=s' );
$support->allowed_params( $support->get_common_params, 'from_assembly',
						  'to_assembly' );
if ( $support->param('help') or $support->error ) {
	warn $support->error if $support->error;
	pod2usage(1);
}
$support->param( 'verbose',     1 );
$support->param( 'interactive', 0 );
my $write_db      = not $support->param('dry_run');
my $from_assembly = $support->param('from_assembly');
my $to_assembly   = $support->param('to_assembly');

throw("Must set from/to_assembly!\n") unless($from_assembly && $to_assembly);

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# first set connection parameters for alternative db
# both databases have to be on the same host, so we don't need to configure
# them separately
for my $prm (qw(host port user pass dbname)) {
	$support->param( "alt$prm", $support->param($prm) )
	  unless ( $support->param("alt$prm") );
}

# database connection
my $dba = $support->get_database( 'ensembl', '' );
my $dbh = $dba->dbc->db_handle;
my $sa  = $dba->get_SliceAdaptor;
#####
# create temporary table for storing non-aligned blocks
#
if ($write_db) {
	$dbh->do(
		qq(
        CREATE TABLE IF NOT EXISTS tmp_align (
          tmp_align_id int(10) unsigned NOT NULL auto_increment,
          alt_seq_region_name varchar(20) NOT NULL,
          alt_start int(10) UNSIGNED NOT NULL,
          alt_end int(10) UNSIGNED NOT NULL,
          ref_seq_region_name varchar(20) NOT NULL,
          ref_start int(10) UNSIGNED NOT NULL,
          ref_end int(10) UNSIGNED NOT NULL,

          PRIMARY KEY (tmp_align_id)
          ) ENGINE=InnoDB;
  )
	);

	# clear tmp_align table of entries from previous runs
	$dbh->do(qq(DELETE FROM tmp_align));
}
else {
	$support->log(
		 "\nHere I would create and empty 'tmp_align' table, if you'd let me\n",
		 1
	);
}

my $component_cs = 'contig';    # 'contig' or 'clone'
my $match         = {};
my $nomatch       = {};
my $seq_region_id = {};
my %stats_total;
my %stats_chr;
my @block_length;
my $fmt1 = "%-40s%10.0f\n";
my $fmt2 = "%-40s%9.2f%%\n";
my $fmt3 = "%-12s%-12s%-12s%-12s%-12s%-9s%-12s%-9s%-12s\n";
my $fmt4 = "%10.0f  %10.0f    %7.0f   %10.0f  %10.0f  %7.0f %-20s %-9s %-20s\n";
my $fmt4b = "%10.0f  %10.0f    %7.0f   %10.0f  %10.0f  %7.0f %-20s %-20s\n";
my $fmt5  = "%-40s%10s\n";
my $fmt6  = "%-10s%-12s%-10s%-12s\n";
my $sth1 = $dbh->prepare(
	qq{
    INSERT IGNORE INTO assembly (asm_seq_region_id, cmp_seq_region_id,
                                 asm_start, asm_end, cmp_start, cmp_end, ori)
    VALUES (?, ?, ?, ?, ?, ?, ?)
}
);
my $sth2 = $dbh->prepare(
	qq{
  INSERT INTO tmp_align values(NULL, ?, ?, ?, ?, ?, ?)
}
);
my $sth3 = $dbh->prepare(
	qq{
		SELECT sa.name AS R_chr,  sc.name AS A_chr, a.*
		FROM assembly a, seq_region sa, seq_region sc, coord_system ca, coord_system cc
		WHERE a.asm_seq_region_id = sa.seq_region_id
		AND a.cmp_seq_region_id = sc.seq_region_id
		AND sa.coord_system_id = ca.coord_system_id
		AND sc.coord_system_id = cc.coord_system_id
		AND cc.name = 'chromosome'
		AND ca.name = 'chromosome'
		AND sa.name = ?
		ORDER BY a.asm_start ASC
	}
);
$support->log_stamped("Looping over chromosomes...\n");
my $chr_list  = $sa->fetch_all('chromosome');
my @from_chrs =
  sort { lc( $a->seq_region_name ) cmp lc( $b->seq_region_name ) }
  grep ( $_->seq_region_name =~ /$from_assembly/, @$chr_list );
my @to_chrs =
  sort { lc( $a->seq_region_name ) cmp lc( $b->seq_region_name ) }
  grep ( $_->seq_region_name =~ /$to_assembly/, @$chr_list );

if ( scalar(@from_chrs) != scalar(@to_chrs) ) {
	throw(   "Chromosome lists do not match by length:\n["
		   . join( " ", map( $_->seq_region_name, @from_chrs ) ) . "]\n["
		   . join( " ", map( $_->seq_region_name, @to_chrs ) )
		   . "]\n" );
}

# Check that the chromosome names match
for my $i ( 0 .. scalar(@from_chrs) - 1 ) {
	my $R_sr_name = $from_chrs[$i]->seq_region_name;
	my $A_sr_name = $to_chrs[$i]->seq_region_name;
	my ($R_chr) = $R_sr_name =~ /(.*)_$from_assembly/;
	my ($A_chr) = $A_sr_name =~ /(.*)_$to_assembly/;
	throw(   "chromosome names don't match $R_chr != $A_chr\n["
		   . join( " , ", map( $_->seq_region_name, @from_chrs ) ) . "]\n["
		   . join( " , ", map( $_->seq_region_name, @to_chrs ) )
		   . "]\n" )
	  unless $R_chr eq $A_chr;
	$support->log_verbose("$R_sr_name	=>	$A_sr_name\n");
	$seq_region_id->{$R_sr_name} = $sa->get_seq_region_id( $from_chrs[$i] );
	$seq_region_id->{$A_sr_name} = $sa->get_seq_region_id( $to_chrs[$i] );
}
CHR: for my $i ( 0 .. scalar(@from_chrs) - 1 ) {

	# get the chromosome slices
	my $R_slice = $from_chrs[$i];
	my $R_chr   = $R_slice->seq_region_name;

	my $A_slice = $to_chrs[$i];
	my $A_chr   = $A_slice->seq_region_name;

	$support->log_stamped( "Chromosome $R_chr/$A_chr ...\n", 1 );
	my @R_components = @{ $R_slice->project($component_cs) };
	$dba->get_AssemblyMapperAdaptor()->delete_cache();
	my @A_components = @{ $A_slice->project($component_cs) };
	$stats_chr{$R_chr}->{'A_cmp_number'} = scalar(@A_components);
	$stats_chr{$R_chr}->{'R_cmp_number'} = scalar(@R_components);
	my $i = 0;
	my $j = 0;
	my ( $left, $right, $type, $tag );
	my $match_flag     = 0;
	my @assembly_diffs = sdiff( \@R_components, \@A_components, \&get_cmp_key );

	# loop over sdiff results
  DIFF: foreach my $diff (@assembly_diffs) {
		$type = $diff->[0];
		( $left, $right, $tag ) = ( '-', '-', '' );
		( $i, $j ) = ( 0, 0 );
		if ( $type eq '+' ) {
			$right = &get_cmp_key( $diff->[2], 1 );
			$tag = '>';
		}
		elsif ( $type eq '-' ) {
			$left = &get_cmp_key( $diff->[1], 1 );
			$tag = '<';
		}
		elsif ( $type eq 'u' ) {
			$left  = &get_cmp_key( $diff->[1], 1 );
			$right = &get_cmp_key( $diff->[2], 1 );
		}
		elsif ( $type eq 'c' ) {
			$left  = &get_cmp_key( $diff->[1], 1 );
			$right = &get_cmp_key( $diff->[2], 1 );
			$tag   = '|';
		}
		if ( $type eq 'u' ) {
			( $i, $j ) = ( 1, 1 );
			found_match(
						 $R_chr,     $A_chr,   $diff->[1],
						 $diff->[2], $i,       $j,
						 $match,     $nomatch, $match_flag
			);
			$stats_chr{$R_chr}->{'identical'}++;
			$match_flag = 1;
		}
		else {
			$i = 1 if $type eq '-';
			$j = 1 if $type eq '+';
			if ( $type eq 'c' ) {
				( $i, $j ) = ( 1, 1 );
				my ( $R_acc, $R_sv ) = split /\./,
				  $diff->[1]->to_Slice->seq_region_name;
				my ( $A_acc, $A_sv ) = split /\./,
				  $diff->[2]->to_Slice->seq_region_name;
				if ( $R_acc eq $A_acc && $R_sv eq $A_sv ) {
					$stats_chr{$R_chr}->{'mismatch'}++;

					my ( $R_chr_start, $R_chr_end ) =
					  ( $diff->[1]->from_start, $diff->[1]->from_end );
					my ( $R_ctg_start, $R_ctg_end ) =
					  ( $diff->[1]->to_Slice->start,
						$diff->[1]->to_Slice->end );
					my $R_ctg_strand = $diff->[1]->to_Slice->strand;
					my $R_coords = [
									 $R_ctg_start, $R_ctg_end, $R_chr_start,
									 $R_chr_end,   $R_ctg_strand
					];

					my ( $A_chr_start, $A_chr_end ) =
					  ( $diff->[2]->from_start, $diff->[2]->from_end );
					my ( $A_ctg_start, $A_ctg_end ) =
					  ( $diff->[2]->to_Slice->start,
						$diff->[2]->to_Slice->end );
					my $A_ctg_strand = $diff->[2]->to_Slice->strand;
					my $A_coords = [
									 $A_ctg_start, $A_ctg_end, $A_chr_start,
									 $A_chr_end,   $A_ctg_strand
					];

					my $blocks = &parse_projections( $R_coords, $A_coords );

					foreach my $block (@$blocks) {
						my ( $R_s, $R_e, $A_s, $A_e, $tag, $ori ) = @$block;
						if ($ori) {
							$support->log("start a new align block\n");
							push @{ $match->{$R_chr} },
							  [ $A_s, $A_e, $j, $R_s, $R_e, $i, $A_chr, $ori ];
						}
						else {
							if ( $tag == -1 ) {
								$support->log(
									   "start a new overlap non-align block\n");
								push @{ $nomatch->{$R_chr} },
								  [ $A_s, $A_e, $j, $R_s, $R_e, $i, $A_chr ];
							}
						}
					}
					$match_flag = 2;
					$support->log(       sprintf( "%-40s\t%-10s %-40s\n", $left, $tag, $right )
					);
					next DIFF;
				}
				elsif ( $R_acc ne $A_acc ) {
					# Project the reference slice to the chromosome cs in case it
					# matches another chromosome
					my $s = $diff->[1]->to_Slice->seq_region_Slice;
					my @A_chrs = @{ $s->project( "chromosome", "Otter" ) };
					foreach $A_chr (@A_chrs) {
						my ($A_chr_name) =
						  $A_chr->to_Slice->seq_region_name =~
						  /(.*)_$to_assembly/;
						if ($A_chr_name) {
							$dba->get_AssemblyMapperAdaptor()->delete_cache();
							my ($proj_chr) =
							  @{ $s->project_to_slice( $A_chr->to_Slice ) };

							my ( $R_chr_start, $R_chr_end ) =
							  ( $diff->[1]->from_start, $diff->[1]->from_end );
							my ( $R_ctg_start, $R_ctg_end ) = (
													$diff->[1]->to_Slice->start,
													$diff->[1]->to_Slice->end
							);
							my $R_ctg_strand = $diff->[1]->to_Slice->strand;
							my $R_coords = [
										 $R_ctg_start, $R_ctg_end, $R_chr_start,
										 $R_chr_end,   $R_ctg_strand
							];

							my ( $A_chr_start, $A_chr_end ) = (
													 $proj_chr->to_Slice->start,
													 $proj_chr->to_Slice->end
							);
							my ( $A_ctg_start, $A_ctg_end ) =
							  ( $proj_chr->from_start, $proj_chr->from_end );
							my $A_ctg_strand = $proj_chr->to_Slice->strand;
							my $A_coords = [
										 $A_ctg_start, $A_ctg_end, $A_chr_start,
										 $A_chr_end,   $A_ctg_strand
							];
							my $blocks =
							  &parse_projections( $R_coords, $A_coords );

							foreach my $block (@$blocks) {
								my ( $R_s, $R_e, $A_s, $A_e, $tag, $ori ) =
								  @$block;
								if ($ori) {
									$support->log("start a new align block\n");
									push @{ $match->{$R_chr} },
									  [
										$A_s, $A_e, $j, $R_s, $R_e, $i,
										$proj_chr->to_Slice->seq_region_name,
										$ori
									  ];
								}
								else {
									if ( $tag == -1 ) {
										$support->log("start a new overlap non-align block\n");
										push @{ $nomatch->{$R_chr} },
										  [
											$A_s,
											$A_e,
											$j,
											$R_s,
											$R_e,
											$i,
											$proj_chr->to_Slice->seq_region_name
										  ];
									}
								}
							}
							$match_flag = 2;
							$support->log(
										   sprintf(
													"%-40s\t%-10s %-40s\n",
													$left, $tag, $right
										   )
							);
							$support->log(
										   sprintf(
													"%-40s\t%-10s %-40s\n",
													"", "", &get_cmp_key( $proj_chr, 1 )
										   )
							);
							next DIFF;
						}
					}
				}
			}
			elsif ( $type eq '-' ) {
				my ( $R_acc, $R_sv ) = split /\./,
				  $diff->[1]->to_Slice->seq_region_name;
				# Project the reference slice to the chromosome cs in case
				# it matches another chromosome
				my $s = $diff->[1]->to_Slice->seq_region_Slice;
				my @A_chrs = @{ $s->project( "chromosome", "Otter" ) };
				foreach $A_chr (@A_chrs) {
					my ($A_chr_name) =
					  $A_chr->to_Slice->seq_region_name =~ /(.*)_$to_assembly/;
					if ($A_chr_name) {
						$dba->get_AssemblyMapperAdaptor()->delete_cache();
						my ($proj_chr) =
						  @{ $s->project_to_slice( $A_chr->to_Slice ) };

						my ( $R_chr_start, $R_chr_end ) =
						  ( $diff->[1]->from_start, $diff->[1]->from_end );
						my ( $R_ctg_start, $R_ctg_end ) = (
													$diff->[1]->to_Slice->start,
													$diff->[1]->to_Slice->end
						);
						my $R_ctg_strand = $diff->[1]->to_Slice->strand;
						my $R_coords = [
										 $R_ctg_start, $R_ctg_end, $R_chr_start,
										 $R_chr_end,   $R_ctg_strand
						];

						my ( $A_chr_start, $A_chr_end ) = (
													 $proj_chr->to_Slice->start,
													 $proj_chr->to_Slice->end
						);
						my ( $A_ctg_start, $A_ctg_end ) =
						  ( $proj_chr->from_start, $proj_chr->from_end );
						my $A_ctg_strand = $proj_chr->to_Slice->strand;
						my $A_coords = [
										 $A_ctg_start, $A_ctg_end, $A_chr_start,
										 $A_chr_end,   $A_ctg_strand
						];

						my $blocks = &parse_projections( $R_coords, $A_coords );

						foreach my $block (@$blocks) {
							my ( $R_s, $R_e, $A_s, $A_e, $tag, $ori ) = @$block;
							if ($ori) {
								$support->log("start a new align block\n");
								push @{ $match->{$R_chr} },
								  [
									$A_s, $A_e, $j, $R_s, $R_e, $i,
									$proj_chr->to_Slice->seq_region_name, $ori
								  ];
							}
							else {
								if ( $tag == -1 ) {
									$support->log("start a new overlap non-align block\n");
									push @{ $nomatch->{$R_chr} },
									  [
										$A_s, $A_e, $j, $R_s, $R_e, $i,
										$proj_chr->to_Slice->seq_region_name
									  ];
								}
							}
						}
						$match_flag = 2;
						$right = &get_cmp_key( $proj_chr, 1 );
						$support->log(
									   sprintf(
												"%-40s\t%-10s %-40s\n",
												$left, $tag, $right
									   )
						);
						next DIFF;
					}
				}
			}
			$match_flag = 0;
		}
		$support->log( sprintf( "%-40s\t%-10s %-40s\n", $left, $tag, $right ) );
	}
}

sub parse_projections {
	my ( $R, $A ) = @_;
	my ( $R_ctg_start, $R_ctg_end, $R_chr_start, $R_chr_end, $R_ctg_strand ) =
	  @$R;
	my ( $A_ctg_start, $A_ctg_end, $A_chr_start, $A_chr_end, $A_ctg_strand ) =
	  @$A;
	my $blocks = [];
	my $ori    = $R_ctg_strand * $A_ctg_strand;

	# check that contig slices overlap
	if ( $A_ctg_end < $R_ctg_start || $R_chr_end < $A_ctg_start ) {
		return $blocks;
	}
	my $start_offset = $R_ctg_start - $A_ctg_start;
	my $end_offset   = $R_ctg_end - $A_ctg_end;
	my ( $R_chr_s, $R_chr_e, $A_chr_s, $A_chr_e );

	# create the 1st non-align block if exists
	if ($start_offset) {
		if ( $R_ctg_strand == 1 ) {
			if ( $start_offset < 0 ) {
				$R_chr_s = $R_chr_start;
				$R_chr_e = $R_chr_start - $start_offset - 1;
			}
			else {
				$R_chr_s = $R_chr_start - $start_offset;
				$R_chr_e = $R_chr_start - 1;
			}
		}
		else {
			if ( $start_offset < 0 ) {
				$R_chr_s = $R_chr_end + $start_offset + 1;
				$R_chr_e = $R_chr_end;
			}
			else {
				$R_chr_s = $R_chr_end + 1;
				$R_chr_e = $R_chr_end + $start_offset;
			}
		}
		if ( $A_ctg_strand == 1 ) {
			if ( $start_offset < 0 ) {
				$A_chr_s = $A_chr_start + $start_offset;
				$A_chr_e = $A_chr_start - 1;
			}
			else {
				$A_chr_s = $A_chr_start;
				$A_chr_e = $A_chr_start + $start_offset - 1;
			}
		}
		else {
			if ( $start_offset < 0 ) {
				$A_chr_s = $A_chr_end + 1;
				$A_chr_e = $A_chr_end - $start_offset;
			}
			else {
				$A_chr_s = $A_chr_end - $start_offset + 1;
				$A_chr_e = $A_chr_end;
			}
		}
		push @$blocks,
		  [
			$R_chr_s, $R_chr_e,
			$A_chr_s, $A_chr_e,
			$start_offset < 0 ? -1 : 1, 0
		  ];
	}

	# create the overlapping block
	( $R_chr_s, $R_chr_e, $A_chr_s, $A_chr_e ) =
	  ( $R_chr_start, $R_chr_end, $A_chr_start, $A_chr_end );
	if ($start_offset) {
		if ( $start_offset < 0 ) {
			if ( $R_ctg_strand == 1 ) {
				$R_chr_s = $R_chr_start - $start_offset;
			}
			else {
				$R_chr_e = $R_chr_end + $start_offset;
			}
		}
		else {
			if ( $A_ctg_strand == 1 ) {
				$A_chr_s = $A_chr_start + $start_offset;
			}
			else {
				$A_chr_e = $A_chr_end - $start_offset;
			}
		}
	}
	if ($end_offset) {
		if ( $end_offset < 0 ) {
			if ( $A_ctg_strand == 1 ) {
				$A_chr_e = $A_chr_end + $end_offset;
			}
			else {
				$A_chr_s = $A_chr_start - $end_offset;
			}
		}
		else {
			if ( $R_ctg_strand == 1 ) {
				$R_chr_e = $R_chr_end - $end_offset;
			}
			else {
				$R_chr_s = $R_chr_start + $end_offset;
			}
		}
	}
	push @$blocks, [ $R_chr_s, $R_chr_e, $A_chr_s, $A_chr_e, 0, $ori ];

	# create the 2nd non-align block if exists
	if ($end_offset) {
		if ( $R_ctg_strand == 1 ) {
			if ( $end_offset > 0 ) {
				$R_chr_s = $R_chr_end - $end_offset + 1;
				$R_chr_e = $R_chr_end;
			}
			else {
				$R_chr_s = $R_chr_end + 1;
				$R_chr_e = $R_chr_end - $end_offset;
			}
		}
		else {
			if ( $end_offset > 0 ) {
				$R_chr_s = $R_chr_start;
				$R_chr_e = $R_chr_start + $end_offset - 1;
			}
			else {
				$R_chr_s = $R_chr_start + $end_offset;
				$R_chr_e = $R_chr_start - 1;
			}
		}
		if ( $A_ctg_strand == 1 ) {
			if ( $end_offset > 0 ) {
				$A_chr_s = $A_chr_end + 1;
				$A_chr_e = $A_chr_end + $end_offset;
			}
			else {
				$A_chr_s = $A_chr_end + $end_offset + 1;
				$A_chr_e = $A_chr_end;
			}
		}
		else {
			if ( $end_offset > 0 ) {
				$A_chr_s = $A_chr_start - $end_offset;
				$A_chr_e = $A_chr_start - 1;
			}
			else {
				$A_chr_s = $A_chr_start;
				$A_chr_e = $A_chr_start - $end_offset - 1;
			}
		}
		push @$blocks,
		  [
			$R_chr_s, $R_chr_e,
			$A_chr_s, $A_chr_e,
			$end_offset > 0 ? -1 : 1, 0
		  ];
	}
	return $blocks;
}
for my $i ( 0 .. scalar(@from_chrs) - 1 ) {

	# get the chromosome slices
	my $R_slice  = $from_chrs[$i];
	my $R_chr    = $R_slice->seq_region_name;
	my $R_length = $R_slice->length;
	my $A_slice  = $to_chrs[$i];
	my $A_chr    = $A_slice->seq_region_name;
	my $A_length = $A_slice->length;
	my $c;
	next unless $match->{$R_chr};

	# store directly aligned blocks in assembly table
	my $number_aligned_blocks = scalar( @{ $match->{$R_chr} || [] } );
	if ($write_db) {
		$support->log(
					 "Adding assembly entries for directly aligned blocks...\n",
					 1 );
		foreach $c ( 0 .. $number_aligned_blocks - 1 ) {
			$sth1->execute(
							$seq_region_id->{$R_chr},
							$seq_region_id->{ $match->{$R_chr}->[$c]->[6] },
							$match->{$R_chr}->[$c]->[3],
							$match->{$R_chr}->[$c]->[4],
							$match->{$R_chr}->[$c]->[0],
							$match->{$R_chr}->[$c]->[1],
							$match->{$R_chr}->[$c]->[7]
			);
		}
		$support->log( "Done inserting $number_aligned_blocks entries.\n", 1 );
	}
	else {
		$support->log(
"\nHere I would insert $number_aligned_blocks rows into 'assembly' table, if you'd let me\n",
			1
		);
	}

# loop through the directly aligned blocks and fill in the non-aligned block hash
# sort the match blocks by ref chromosome start
	my ( $r_start, $r_end, $a_start, $a_end, $a_chr, $ref_gap, $alt_gap );
	my ( $sr_start, $sr_end, $sa_start, $sa_end ) = ( 0, 0, 0, 0 );
	my $last = [ $A_length + 1, 0, 0, $R_length + 1, 0, 0, $A_chr ];
	my $sa_chr = $A_chr;
	foreach ( sort { $a->[3] <=> $b->[3] } @{ $match->{$R_chr} }, $last ) {
		$a_start = $_->[0];
		$a_end   = $_->[1];
		$r_start = $_->[3];
		$r_end   = $_->[4];
		$a_chr   = $_->[6];
		if ( $a_chr eq $sa_chr ) {
			$ref_gap = $r_start - $sr_end - 1;
			$alt_gap = $a_start - $sa_end - 1;
			if (    ( $ref_gap > 0 && $alt_gap > 0 )
				 && ( $ref_gap < 1000000 && $alt_gap < 1000000 ) )
			{
				push @{ $nomatch->{$R_chr} },
				  [
					$sa_end + 1,
					$a_start - 1,
					1,
					$sr_end + 1,
					$r_start - 1,
					1,
					$a_chr,
				  ];
			}
		}
		$sa_start = $_->[0];
		$sa_end   = $_->[1];
		$sr_start = $_->[3];
		$sr_end   = $_->[4];
		$sa_chr   = $_->[6];
	}

	# filter single assembly inserts from non-aligned blocks (i.e. cases where
	# a block has components only in one assembly, not in the other) - there is
	# nothing to be done with them
	@{ $nomatch->{$R_chr} } =
	  grep { $_->[2] > 0 and $_->[5] > 0 } @{ $nomatch->{$R_chr} }
	  if ( $nomatch->{$R_chr} );

	# store non-aligned blocks in tmp_align table
	my $number_nonaligned_blocks = scalar( @{ $nomatch->{$R_chr} || [] } );
	if ($write_db) {
		if ( $nomatch->{$R_chr} ) {
			$support->log( "Storing non-aligned blocks in tmp_align table...\n",
						   1 );
			foreach $c ( 0 .. $number_nonaligned_blocks - 1 ) {
				$sth2->execute(
								$nomatch->{$R_chr}->[$c]->[6],
								$nomatch->{$R_chr}->[$c]->[0],
								$nomatch->{$R_chr}->[$c]->[1],
								$R_chr,
								$nomatch->{$R_chr}->[$c]->[3],
								$nomatch->{$R_chr}->[$c]->[4],
				);
			}
			$support->log(
						  "Done inserting $number_nonaligned_blocks entries.\n",
						  1 );
		}
	}
	else {
		$support->log(
"\nHere I would insert $number_nonaligned_blocks rows into 'tmp_align' table, if you'd let me\n",
			1
		);
	}

	# stats for this chromosome
	$stats_chr{$R_chr}->{'A_only'} =
	  $stats_chr{$R_chr}->{'A_cmp_number'} - $stats_chr{$R_chr}->{'identical'} -
	  $stats_chr{$R_chr}->{'mismatch'};
	$stats_chr{$R_chr}->{'R_only'} =
	  $stats_chr{$R_chr}->{'R_cmp_number'} - $stats_chr{$R_chr}->{'identical'} -
	  $stats_chr{$R_chr}->{'mismatch'};
	for ( $c = 0 ; $c < scalar( @{ $match->{$R_chr} || [] } ) ; $c++ ) {
		$stats_chr{$R_chr}->{'A_matchlength'} +=
		  $match->{$R_chr}->[$c]->[1] - $match->{$R_chr}->[$c]->[0];
		$stats_chr{$R_chr}->{'R_matchlength'} +=
		  $match->{$R_chr}->[$c]->[4] - $match->{$R_chr}->[$c]->[3];
	}
	$stats_chr{$R_chr}->{'A_coverage'} =
	  100 * $stats_chr{$R_chr}->{'A_matchlength'} / $A_length;
	$stats_chr{$R_chr}->{'R_coverage'} =
	  100 * $stats_chr{$R_chr}->{'R_matchlength'} / $R_length;
	map { $stats_total{$_} += $stats_chr{$R_chr}->{$_} } keys %stats_chr;
	$support->log( "\nStats for chromosome $R_chr:\n\n", 1 );
	$support->log( sprintf( $fmt5, "Alternative chromosome name:", $A_chr ),
				   2 );
	$support->log( sprintf( $fmt1, "Length (alternative):", $A_length ), 2 );
	$support->log( sprintf( $fmt1, "Length (reference):",   $R_length ), 2 );
	$support->log(
				   sprintf( $fmt1,
							"Identical components:",
							$stats_chr{$R_chr}->{'identical'} ),
				   2
	);
	$support->log(
				   sprintf( $fmt1,
							"Identical components that were skipped:",
							$stats_chr{$R_chr}->{'skipped'} ),
				   2
	);
	$support->log(
				   sprintf( $fmt1,
							"Components with start/end mismatch:",
							$stats_chr{$R_chr}->{'mismatch'} ),
				   2
	);
	$support->log(
				   sprintf( $fmt1,
							"Components only in alternative assembly:",
							$stats_chr{$R_chr}->{'A_only'} ),
				   2
	);
	$support->log(
				   sprintf( $fmt1,
							"Components only in reference assembly:",
							$stats_chr{$R_chr}->{'R_only'} ),
				   2
	);
	$support->log(
				   sprintf( $fmt2,
							"Direct match coverage (alternative):",
							$stats_chr{$R_chr}->{'A_coverage'} ),
				   2
	);
	$support->log(
				   sprintf( $fmt2,
							"Direct match coverage (reference):",
							$stats_chr{$R_chr}->{'R_coverage'} ),
				   2
	);

	# Aligned blocks
	if ( $match->{$R_chr} ) {
		$support->log( "\nDirectly aligned blocks:\n\n", 1 );
		$support->log(
			sprintf(
				$fmt3,
				qw(ALT_START ALT_END ALT_COMPONENTS REF_START REF_END REF_COMPONENTS ALT_CHR ORI REF_CHR)
			),
			2
		);
		$support->log( ( '-' x 71 ) . "\n", 2 );
		for ( $c = 0 ; $c < scalar( @{ $match->{$R_chr} } ) ; $c++ ) {
			$support->log(
						   sprintf(
									$fmt4, @{ $match->{$R_chr}->[$c] }, $R_chr
						   ),
						   2
			);

			# sanity check: aligned region pairs must have same length
			my $e_len =
			  $match->{$R_chr}->[$c]->[1] - $match->{$R_chr}->[$c]->[0] + 1;
			my $v_len =
			  $match->{$R_chr}->[$c]->[4] - $match->{$R_chr}->[$c]->[3] + 1;
			$support->log_warning( "Length mismatch: $e_len <> $v_len\n", 2 )
			  unless ( $e_len == $v_len );
		}
	}

	# Non-aligned blocks
	if ( $nomatch->{$R_chr} ) {
		$support->log( "\nNon-aligned blocks:\n\n", 1 );
		$support->log(
			sprintf(
				$fmt3,
				qw(ALT_START ALT_END ALT_COMPONENTS REF_START REF_END REF_COMPONENTS ALT_CHR REF_CHR)
			),
			2
		);
		$support->log( ( '-' x 71 ) . "\n", 2 );
		for ( $c = 0 ; $c < scalar( @{ $nomatch->{$R_chr} } ) ; $c++ ) {
			$support->log(
						   sprintf( $fmt4b,
									@{ $nomatch->{$R_chr}->[$c] }, $R_chr ),
						   2
			);

			# find longest non-aligned block
			my $A_length =
			  $nomatch->{$R_chr}->[$c]->[1] - $nomatch->{$R_chr}->[$c]->[0] + 1;
			my $R_length =
			  $nomatch->{$R_chr}->[$c]->[4] - $nomatch->{$R_chr}->[$c]->[3] + 1;
			push @block_length,
			  [ $nomatch->{$R_chr}->[$c]->[6], $A_length, $R_chr, $R_length ];
		}
	}
	$support->log_stamped( "\nDone with chromosome $R_chr.\n", 1 );

	#last CHR;
}

# overall stats
$support->log("\nOverall stats:\n");
$support->log(
			   sprintf( $fmt1,
						"Identical components:",
						$stats_total{'identical'} ),
			   1
);
$support->log(
			   sprintf( $fmt1,
						"Identical components that were skipped:",
						$stats_total{'skipped'} ),
			   1
);
$support->log(
			   sprintf( $fmt1,
						"Components with start/end mismatch:",
						$stats_total{'mismatch'} ),
			   1
);
$support->log(
			   sprintf( $fmt1,
						"Components only in alternative assembly:",
						$stats_total{'A_only'} ),
			   1
);
$support->log(
			   sprintf( $fmt1,
						"Components only in reference assembly:",
						$stats_total{'R_only'} ),
			   1
);
$support->log("\nNon-match block lengths:\n");
$support->log( sprintf( $fmt6, qw(ALT_CHR ALT_LENGTH REF_CHR REF_LENGTH) ), 1 );
$support->log( ( '-' x 42 ) . "\n", 1 );

foreach my $block ( sort { $a->[1] <=> $b->[1] } @block_length ) {
	$support->log( sprintf( "%-10s%10.0f  %-10s%10.0f\n", @{$block} ), 1 );
}
$support->log_stamped("\nDone.\n");

# finish logfile
$support->finish_log;
### end main
sub get_cmp_key {
	my ( $cmp, $flag ) = @_;
	my $slice = $cmp->to_Slice;
	my $key   =
	    $slice->seq_region_name . ":"
	  . $slice->start . "-"
	  . $slice->end . ":"
	  . $slice->strand;
	$key .= ":"
	  . $cmp->from_start . "-"
	  . $cmp->from_end . ":"
	  . ( $cmp->from_end - $cmp->from_start + 1 )
	  if $flag;
	return $key;
}

=head2 found_match

  Arg[1]      : String $R_chr - reference chromosome name
  Arg[2]      : String $A_chr - alternative chromosome name
  Arg[3]      : Bio::EnsEMBL::ProjectionSegment $R_seg - current reference
                segment
  Arg[4]      : Bio::EnsEMBL::ProjectionSegment $A_seg - current alternative
                segment
  Arg[5]      : Boolean $i - indicates if there is a reference component
  Arg[6]      : Boolean $j - indicates if there is an alternative component
  Arg[7]      : Hashref $match - datastructure to store aligned blocks
  Arg[8]      : Hashref $nomatch - datastructure to store non-aligned blocks
  Arg[9]      : Boolean $match_flag - flag indicating if last component was a match

  Description : This function is called when two components match (i.e. have the
                same name.version and seq_reg start and end in both assemblies). Depending on the state
                of the last component (match or nomatch), it extends aligned blocks
                or creates a new aligned block.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_match {
	my ( $R_chr, $A_chr, $R_seg, $A_seg, $i, $j, $match, $nomatch, $match_flag )
	  = @_;
	my $R_start = $R_seg->from_start;
	my $R_end   = $R_seg->from_end;
	my $A_start = $A_seg->from_start;
	my $A_end   = $A_seg->from_end;

	# last component was a match
	if ( $match_flag == 1 ) {

		# adjust align block end
		if ( $match->{$R_chr} ) {

			# check that there is no single gap in the block
			my $c     = scalar( @{ $match->{$R_chr} } ) - 1;
			my $A_gap = $A_start - $match->{$R_chr}->[$c]->[1];
			my $R_gap = $R_start - $match->{$R_chr}->[$c]->[4];
			if ( $A_gap == $R_gap ) {
				$support->log("adjust align block end\n");
				$match->{$R_chr}->[$c]->[1] = $A_end;
				$match->{$R_chr}->[$c]->[2]++ if $j;
				$match->{$R_chr}->[$c]->[4] = $R_end;
				$match->{$R_chr}->[$c]->[5]++ if $i;
			}
			else {
				$support->log("start a new align block (because of gap)\n");
				push @{ $match->{$R_chr} },
				  [ $A_start, $A_end, $j, $R_start, $R_end, $i, $A_chr, 1 ];
			}
		}
	}

	# last component was a non-match
	else {

		# start a new align block
		$support->log("start a new align block\n");
		push @{ $match->{$R_chr} },
		  [ $A_start, $A_end, $j, $R_start, $R_end, $i, $A_chr, 1 ];
	}
}

=head2 found_nomatch
  Arg[1]      : String $R_chr - reference chromosome name
  Arg[2]      : String $A_chr - alternative chromosome name
  Arg[3]      : Bio::EnsEMBL::ProjectionSegment $R_seg - current reference
                segment
  Arg[4]      : Bio::EnsEMBL::ProjectionSegment $A_seg - current alternative
                segment
  Arg[5]      : Boolean $i - indicates if there is a reference component
  Arg[6]      : Boolean $j - indicates if there is an alternative component
  Arg[7]      : Hashref $match - datastructure to store aligned blocks
  Arg[8]      : Hashref $nomatch - datastructure to store non-aligned blocks
  Arg[9]      : Boolean $match_flag - flag indicating if last component was a match
  Description : This function is called when two components don't match (either
                different name.version or length mismatch in the two
                assemblies). Depending on the state of the last component (nomatch
                or match), it extends non-aligned blocks or creates a new non-aligned block.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_nomatch {
	my ( $R_chr, $A_chr, $R_seg, $A_seg, $i, $j, $match, $nomatch, $match_flag )
	  = @_;
	my ( $R_start, $R_end, $A_start, $A_end ) = ( 0, 0, 0, 0 );
	if ($i) {
		$R_start = $R_seg->from_start;
		$R_end   = $R_seg->from_end;
	}
	if ($j) {
		$A_start = $A_seg->from_start;
		$A_end   = $A_seg->from_end;
	}

	# last component was a match
	if ($match_flag) {

		# start a new non-align block
		$support->log("start a new non-align block\n");
		push @{ $nomatch->{$R_chr} },
		  [ $A_start, $A_end, $j, $R_start, $R_end, $i, $A_chr, ];
	}

	# last component was a non-match
	else {

		# adjust non-align block end
		if ( $nomatch->{$R_chr} ) {
			$support->log("adjust non-align block end\n");
			my $c = scalar( @{ $nomatch->{$R_chr} || [] } ) - 1;
			$nomatch->{$R_chr}->[$c]->[0] ||= $A_start;
			$nomatch->{$R_chr}->[$c]->[1] = $A_end if $A_end;
			$nomatch->{$R_chr}->[$c]->[2]++ if $j;
			$nomatch->{$R_chr}->[$c]->[3] ||= $R_start;
			$nomatch->{$R_chr}->[$c]->[4] = $R_end if $R_end;
			$nomatch->{$R_chr}->[$c]->[5]++ if $i;
		}
		else {
			$support->log("start a new non-align block\n");
			push @{ $nomatch->{$R_chr} },
			  [ $A_start, $A_end, $j, $R_start, $R_end, $i, $A_chr, ];
		}
	}
}
