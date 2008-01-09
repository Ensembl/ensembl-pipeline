#!/usr/local/bin/perl


my $component_cs = 'contig';    # 'contig' or 'clone'

=head1 NAME

align_by_component_identity.pl - create a whole genome alignment between two closely
related assemblies, step 1

=head1 SYNOPSIS

align_by_component_identity.pl [arguments]

Required arguments:

    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --dbname, db_name=NAME              database name NAME
    --assembly=ASSEMBLY                 assembly version ASSEMBLY
    --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional arguments:

    --altdbname=NAME                    alternative database NAME

    --chromosomes, --chr=LIST           only process LIST chromosomes
    --altchromosomes, --altchr=LIST     supply alternative chromosome names (the two lists must agree)
    --ref_start                         start coordinate on reference chromosomes
    --ref_end                           end coordinate on reference chromosomes
    --alt_start                         start coordinate on alternative chromosomes
    --alt_end                           end coordinate on alternative chromosomes

    --skipcomponents=FILE               read list of components to skip from FILE

    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)

    -v, --verbose=0|1                   verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry=0|1, --dry_run=0|1        don't write results to database
    -h, --help, -?                      print help (this message)
    --exctype=PAR|HAP                   if defined, will modify 'assembly_exception' table instead of 'assembly'
                                        ( encode regions are to be considered 'PAR' regions for the time being )

=head1 DESCRIPTION

This script is part of a series of scripts to create a mapping between two
assemblies. It assembles the chromosome coordinate systems of two different
assemblies of a genome by creating a whole genome alignment between the two.

The process assumes that the two assemblies are reasonably similar, i.e. there
are no major rearrangements or components moved from one chromosome to another.

See "Related scripts" below for an overview of the whole process.

This particular script creates a whole genome alignment between two closely
related assemblies. You will need a database containing the reference assembly
and the alternative chromosomes which can be created using
load_alternative_assembly.pl.

The alignment is created in two steps:

    1. Match components with same name and version directly and create alignment
       blocks for these regions. Components can be tagged manually to be excluded
       from these direct matches by listing them in a file of components to skip
       (--skipcomponents argument). This can be useful to get better results in
       regions with major assembly differences.

       The result is stored in the assembly table as an assembly between the
       chromosomes of both genome assemblies.

    2. Store non-aligned blocks in a temporary table (tmp_align). They can
       later be aligned using blastz by align_nonident_regions.pl.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two assemblies
is done by a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

modified by Leo Gordon <lg4@sanger.ac.uk>

and Mustapha Larbaoui <ml6@sanger.ac.uk>

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
	$SERVERROOT = "$Bin/../../..";
	unshift( @INC, "$SERVERROOT/ensembl_HEAD/modules" );
	unshift( @INC, "$SERVERROOT/bioperl-live" );
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Attribute;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
	'assembly=s',               'altdbname=s',
	'altassembly=s',            'chromosomes|chr=s@',
	'altchromosomes|altchr=s@', 'skipcomponents|skip_components=s',
	'exctype=s', 'ref_start=i', 'ref_end=i', 'alt_start=i', 'alt_end=i'
);
$support->allowed_params( $support->get_common_params, 'assembly', 'altdbname',
	'altassembly', 'chromosomes', 'altchromosomes', 'skipcomponents', 'exctype',
);

if ( $support->param('help') or $support->error ) {
	warn $support->error if $support->error;
	pod2usage(1);
}

$support->comma_to_list( 'chromosomes', 'altchromosomes' );

$support->param( 'verbose',     0 );    # throw away all that garbage
$support->param( 'interactive', 0 );    # stop that garbage from coming up

my $write_db = not $support->param('dry_run');

my $exctype = $support->param('exctype') || '';

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params( 'assembly', 'altassembly' );

# first set connection parameters for alternative db
# both databases have to be on the same host, so we don't need to configure
# them separately
for my $prm (qw(host port user pass dbname)) {
	$support->param( "alt$prm", $support->param($prm) )
	  unless ( $support->param("alt$prm") );
}

# reference database
my $R_dba = $support->get_database( 'ensembl', '' );
my $R_dbh = $R_dba->dbc->db_handle;
my $R_sa  = $R_dba->get_SliceAdaptor;
my $R_asm = $support->param('assembly');
my $R_start = $support->param('ref_start') || undef;
my $R_end = $support->param('ref_end') || undef;

# database containing the alternative assembly
my $A_dba = $support->get_database( 'ensembl', 'alt' );
my $A_sa  = $A_dba->get_SliceAdaptor;
my $A_asm = $support->param('altassembly');
my $A_start = $support->param('alt_start') || undef;
my $A_end = $support->param('alt_end') || undef;


#####
# create temporary table for storing non-aligned blocks
#
if ($write_db) {
	$R_dbh->do(
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
	$R_dbh->do(qq(DELETE FROM tmp_align));
}
else {
	$support->log(
		"\nHere I would create and empty 'tmp_align' table, if you'd let me\n",
		1
	);
}


#####
# read list of components to skip from file
#
$support->log("Reading list of components to skip from file...\n");
my %skip = ();
if ( $support->param('skipcomponents') ) {
	my $infh = $support->filehandle( '<', $support->param('skipcomponents') );
	while (<$infh>) {
		chomp;
		$skip{$_} = 1;
	}
}
$support->log("Done.\n\n");

my $match   = {};
my $nomatch = {};
my %stats_total;
my @block_length;

my $fmt1 = "%-40s%10.0f\n";
my $fmt2 = "%-40s%9.2f%%\n";
my $fmt3 = "%-12s%-12s%-12s%-12s%-12s%-9s\n";
my $fmt4 = "%10.0f  %10.0f    %7.0f   %10.0f  %10.0f  %7.0f\n";
my $fmt5 = "%-40s%10s\n";
my $fmt6 = "%-10s%-12s%-10s%-12s\n";

my $sth1 = $R_dbh->prepare(
	$exctype
	? qq{
    INSERT IGNORE INTO assembly_exception (exc_seq_region_id, seq_region_id,
                                           exc_seq_region_start, exc_seq_region_end,
                                           seq_region_start, seq_region_end,
                                           ori, exc_type)
    VALUES (?, ?, ?, ?, ?, ?, 1, '$exctype')
}
	: qq{
    INSERT IGNORE INTO assembly (asm_seq_region_id, cmp_seq_region_id,
                                 asm_start, asm_end, cmp_start, cmp_end, ori)
    VALUES (?, ?, ?, ?, ?, ?, 1)
}
);

my $sth2 = $R_dbh->prepare(
	qq{
  INSERT INTO tmp_align values(NULL, ?, ?, ?, ?, ?, ?)
}
);

$support->log_stamped("Looping over chromosomes...\n");

my @R_chr_list = $support->param('chromosomes');
if ( !scalar(@R_chr_list) ) {
	@R_chr_list = $support->sort_chromosomes;

	if ( scalar( $support->param('altchromosomes') ) ) {
		die "AltChromosomes list is defined while Chromosomes list is not!";
	}
}

my @A_chr_list = $support->param('altchromosomes');
if ( !scalar(@A_chr_list) ) {
	@A_chr_list = @R_chr_list;
}
elsif ( scalar(@R_chr_list) != scalar(@A_chr_list) ) {
	die "Chromosome lists do not match by length";
}

for my $i ( 0 .. scalar(@R_chr_list) - 1 ) {
	my $R_chr = $R_chr_list[$i];
	my $A_chr = $A_chr_list[$i];

	$support->log_stamped( "Chromosome $R_chr/$A_chr ...\n", 1 );

	# fetch chromosome slices
	my $R_slice =
	  $R_sa->fetch_by_region( 'chromosome', $R_chr, $R_start, $R_end, undef,
		$R_asm );
	print STDOUT $R_slice->seq_region_name." ".$R_slice->start." -> ".$R_slice->end."\n";
	my $A_slice =
	  $A_sa->fetch_by_region( 'chromosome', $A_chr, $A_start, $A_end, undef,
		$A_asm );
	my $A_slice_ref =
	  $A_sa->fetch_by_region( 'chromosome', $A_chr, $A_start, $A_end, undef,
		$A_asm );


	my @A_components = @{ $A_slice->project($component_cs) };
	my @R_components = @{ $R_slice->project($component_cs) };

	if($A_start){
		map ($_->[0] += ($A_start-1) , @A_components);
		map ($_->[1] += ($A_start-1) , @A_components);
	}

	if($R_start){
		map ($_->[0] += ($R_start-1) , @R_components);
		map ($_->[1] += ($R_start-1) , @R_components);
	}


	# loop over alternative components
	my $last       = 0;
	my $j          = 0;
	my $match_flag = 0;
	my $last_A_seg;
	my %stats_chr;

	foreach my $A_seg (@A_components) {

		my $A_component = $A_seg->to_Slice;

		print STDOUT
			"Alternative component ($j) "
			  . $A_component->seq_region_name . ":"
			  . $A_component->start . "-"
			  . $A_component->end . ":"
			  . $A_component->strand
			  . " $A_chr:"
			  . $A_seg->from_start . "-"
			  . $A_seg->from_end . "\n";

		# walk reference components
	  REFCOMPONENTS:
		for ( my $i = $last ; $i < scalar(@R_components) ; $i++ ) {

			my $R_component = $R_components[$i]->to_Slice;

			print STDOUT "Reference component ($i) "
			  . $R_component->seq_region_name . ":"
			  . $R_component->start . "-"
			  . $R_component->end . ":"
			  . $R_component->strand
			  . " $R_chr:"
			  . $R_components[$i]->from_start . "-"
			  . $R_components[$i]->from_end . "\n";

			# same name.version and strand found
			if ( $A_component->seq_region_name eq $R_component->seq_region_name
				and $A_component->strand == $R_component->strand )
			{

				# same component start/end -> identical assembly
				if (    $A_component->start == $R_component->start
					and $A_component->end == $R_component->end )
				{

					# check if component is tagged to be skipped
					# this can be used to resolve some odd assembly differences
					if ( $skip{ $A_component->seq_region_name } ) {

						$support->log_verbose(
							"Skipping matching reference component ($i)"
							  . $R_component->seq_region_name . ":"
							  . $R_component->start . "-"
							  . $R_component->end . ":"
							  . $R_component->strand
							  . "$R_chr:"
							  . $R_components[$i]->from_start . "-"
							  . $R_components[$i]->from_end . "\n",
							2
						);

						&found_nomatch(
							$R_chr,            $A_chr,
							$match,            $nomatch,
							$A_seg,            $last_A_seg,
							$R_components[$i], $R_components[ $i - 1 ],
							$match_flag,       $i,
							$j
						);

						$stats_chr{'skipped'}++;
						$match_flag = 0;

					}
					else {

						$support->log_verbose(
							"Found matching reference component ($i)"
							  . $R_component->seq_region_name . ":"
							  . $R_component->start . "-"
							  . $R_component->end . ":"
							  . $R_component->strand
							  . "$R_chr:"
							  . $R_components[$i]->from_start . "-"
							  . $R_components[$i]->from_end . "\n",
							2
						);
						print STDOUT "Found matching reference component ($i)"
						  . $R_component->seq_region_name . ":"
						  . $R_component->start . "-"
						  . $R_component->end . ":"
						  . $R_component->strand
						  . "$R_chr:"
						  . $R_components[$i]->from_start . "-"
						  . $R_components[$i]->from_end . "\n";
						&found_match(
							$R_chr,            $A_chr,
							$match,            $nomatch,
							$A_seg,            $last_A_seg,
							$R_components[$i], $R_components[ $i - 1 ],
							$match_flag,       $i,
							$j
						);

						$stats_chr{'identical'}++;
						$match_flag = 1;
					}

					# start/end mismatch
				}
				else {

					$support->log_verbose(
						"Start/end mismatch for component ($i) "
						  . $R_component->seq_region_name . ":"
						  . $R_component->start . "-"
						  . $R_component->end . ":"
						  . $R_component->strand
						  . " $R_chr:"
						  . $R_components[$i]->from_start . "-"
						  . $R_components[$i]->from_end . "\n",
						2
					);

					print STDOUT "Start/end mismatch for component ($i) "
						  . $R_component->seq_region_name . ":"
						  . $R_component->start . "-"
						  . $R_component->end . ":"
						  . $R_component->strand
						  . " $R_chr:"
						  . $R_components[$i]->from_start . "-"
						  . $R_components[$i]->from_end . "\n";

					&found_nomatch(
						$R_chr,            $A_chr,
						$match,            $nomatch,
						$A_seg,            $last_A_seg,
						$R_components[$i], $R_components[ $i - 1 ],
						$match_flag,       $i,
						$j
					);

					$stats_chr{'mismatch'}++;
					$match_flag = 0;
				}
				$i++;
				$last = $i;
				last REFCOMPONENTS;

				# different components
			}
			else {

				$support->log_verbose(
					"Skipping component ($i)"
					  . $R_component->seq_region_name . ":"
					  . $R_component->start . "-"
					  . $R_component->end . ":"
					  . $R_component->strand
					  . " $R_chr:"
					  . $R_components[$i]->from_start . "-"
					  . $R_components[$i]->from_end . "\n",
					2
				);

				print STDOUT "Skipping component ($i)"
					  . $R_component->seq_region_name . ":"
					  . $R_component->start . "-"
					  . $R_component->end . ":"
					  . $R_component->strand
					  . " $R_chr:"
					  . $R_components[$i]->from_start . "-"
					  . $R_components[$i]->from_end . "\n";

				&found_nomatch(
					$R_chr,            $A_chr,
					$match,            $nomatch,
					$A_seg,            $last_A_seg,
					$R_components[$i], $R_components[ $i - 1 ],
					$match_flag,       $i,
					$j
				);

				$match_flag = 0;

			}
		}

		$last_A_seg = $A_seg;
		$j++;
	}

	# adjust the final component count
	if ($match_flag) {

		# last component was a match, adjust matching component count
		if ( $match->{$R_chr} ) {

			my $c = scalar( @{ $match->{$R_chr} } ) - 1;
			$match->{$R_chr}->[$c]->[2] =
			  scalar(@A_components) - $match->{$R_chr}->[$c]->[2];
			$match->{$R_chr}->[$c]->[5] =
			  scalar(@R_components) - $match->{$R_chr}->[$c]->[5];

		}

	}
	else {

		# last component was a non-match, adjust non-matching component count
		if ( $nomatch->{$R_chr} ) {

			my $c = scalar( @{ $nomatch->{$R_chr} } ) - 1;
			$nomatch->{$R_chr}->[$c]->[2] =
			  scalar(@A_components) - $nomatch->{$R_chr}->[$c]->[2];
			$nomatch->{$R_chr}->[$c]->[5] =
			  scalar(@R_components) - $nomatch->{$R_chr}->[$c]->[5];

		}

	}

	# filter single assembly inserts from non-aligned blocks (i.e. cases where
	# a block has components only in one assembly, not in the other) - there is
	# nothing to be done with them
	@{ $nomatch->{$R_chr} } =
	  grep { $_->[2] > 0 and $_->[5] > 0 } @{ $nomatch->{$R_chr} }
	  if ( $nomatch->{$R_chr} );

	# store directly aligned blocks in assembly table
	my $number_aligned_blocks = scalar( @{ $match->{$R_chr} || [] } );
	if ($write_db) {

		$support->log(
			"Adding assembly entries for directly aligned blocks...\n", 1 );

		foreach my $c ( 0 .. $number_aligned_blocks - 1 ) {
			$sth1->execute(
				$R_sa->get_seq_region_id($R_slice),
				$R_sa->get_seq_region_id($A_slice_ref),
				$match->{$R_chr}->[$c]->[3],
				$match->{$R_chr}->[$c]->[4],
				$match->{$R_chr}->[$c]->[0],
				$match->{$R_chr}->[$c]->[1]
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

	# store non-aligned blocks in tmp_align table
	my $number_nonaligned_blocks = scalar( @{ $nomatch->{$R_chr} || [] } );
	if ($write_db) {

		if ( $nomatch->{$R_chr} ) {

			$support->log( "Storing non-aligned blocks in tmp_align table...\n",
				1 );

			foreach my $c ( 0 .. $number_nonaligned_blocks - 1 ) {
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
				"Done inserting $number_nonaligned_blocks entries.\n", 1 );
		}
	}
	else {
		$support->log(
"\nHere I would insert $number_nonaligned_blocks rows into 'tmp_align' table, if you'd let me\n",
			1
		);
	}

	# stats for this chromosome
	$stats_chr{'A_only'} =
	  scalar(@A_components) - $stats_chr{'identical'} - $stats_chr{'mismatch'};
	$stats_chr{'R_only'} =
	  scalar(@R_components) - $stats_chr{'identical'} - $stats_chr{'mismatch'};
	for ( my $c = 0 ; $c < scalar( @{ $match->{$R_chr} || [] } ) ; $c++ ) {
		$stats_chr{'A_matchlength'} +=
		  $match->{$R_chr}->[$c]->[1] - $match->{$R_chr}->[$c]->[0];
		$stats_chr{'R_matchlength'} +=
		  $match->{$R_chr}->[$c]->[4] - $match->{$R_chr}->[$c]->[3];
	}
	$stats_chr{'A_coverage'} =
	  100 * $stats_chr{'A_matchlength'} / $A_slice->length;
	$stats_chr{'R_coverage'} =
	  100 * $stats_chr{'R_matchlength'} / $R_slice->length;
	map { $stats_total{$_} += $stats_chr{$_} } keys %stats_chr;

	$support->log( "\nStats for chromosome $R_chr:\n\n", 1 );
	$support->log( sprintf( $fmt5, "Alternative chromosome name:", $A_chr ),
		2 );
	$support->log( sprintf( $fmt1, "Length (alternative):", $A_slice->length ),
		2 );
	$support->log( sprintf( $fmt1, "Length (reference):", $R_slice->length ),
		2 );
	$support->log(
		sprintf( $fmt1, "Identical components:", $stats_chr{'identical'} ), 2 );
	$support->log(
		sprintf( $fmt1,
			"Identical components that were skipped:",
			$stats_chr{'skipped'} ),
		2
	);
	$support->log(
		sprintf( $fmt1,
			"Components with start/end mismatch:",
			$stats_chr{'mismatch'} ),
		2
	);
	$support->log(
		sprintf( $fmt1,
			"Components only in alternative assembly:",
			$stats_chr{'A_only'} ),
		2
	);
	$support->log(
		sprintf( $fmt1,
			"Components only in reference assembly:",
			$stats_chr{'R_only'} ),
		2
	);
	$support->log(
		sprintf( $fmt2,
			"Direct match coverage (alternative):",
			$stats_chr{'A_coverage'} ),
		2
	);
	$support->log(
		sprintf( $fmt2,
			"Direct match coverage (reference):",
			$stats_chr{'R_coverage'} ),
		2
	);

	# Aligned blocks
	if ( $match->{$R_chr} ) {

		$support->log( "\nDirectly aligned blocks:\n\n", 1 );
		$support->log(
			sprintf( $fmt3,
				qw(ALT_START ALT_END ALT_COMPONENTS REF_START REF_END REF_COMPONENTS)
			),
			2
		);
		$support->log( ( '-' x 71 ) . "\n", 2 );

		for ( my $c = 0 ; $c < scalar( @{ $match->{$R_chr} } ) ; $c++ ) {

			$support->log( sprintf( $fmt4, @{ $match->{$R_chr}->[$c] } ), 2 );

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
			sprintf( $fmt3,
				qw(ALT_START ALT_END ALT_COMPONENTS REF_START REF_END REF_COMPONENTS)
			),
			2
		);
		$support->log( ( '-' x 71 ) . "\n", 2 );

		for ( my $c = 0 ; $c < scalar( @{ $nomatch->{$R_chr} } ) ; $c++ ) {

			$support->log( sprintf( $fmt4, @{ $nomatch->{$R_chr}->[$c] } ), 2 );

			# find longest non-aligned block
			my $A_length =
			  $nomatch->{$R_chr}->[$c]->[1] - $nomatch->{$R_chr}->[$c]->[0] + 1;
			my $R_length =
			  $nomatch->{$R_chr}->[$c]->[4] - $nomatch->{$R_chr}->[$c]->[3] + 1;
			push @block_length, [ $A_chr, $A_length, $R_chr, $R_length ];

		}
	}

	$support->log_stamped( "\nDone with chromosome $R_chr.\n", 1 );
}

# overall stats
$support->log("\nOverall stats:\n");
$support->log(
	sprintf( $fmt1, "Identical components:", $stats_total{'identical'} ), 1 );
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

=head2 found_match

  Arg[1]      : String $R_chr - reference chromosome name
  Arg[2]      : String $A_chr - alternative chromosome name
  Arg[3]      : Hashref $match - datastructure to store aligned blocks
  Arg[4]      : Hashref $nomatch - datastructure to store non-aligned blocks
  Arg[5]      : Bio::EnsEMBL::ProjectionSegment $A_seg - current alternative
                segment
  Arg[6]      : Bio::EnsEMBL::ProjectionSegment $last_A_seg - last alternative
                segment
  Arg[7]      : Bio::EnsEMBL::ProjectionSegment $R_seg - current reference
                segment
  Arg[8]      : Bio::EnsEMBL::ProjectionSegment $last_R_seg - last reference
                segment
  Arg[9]      : Boolean $match_flag - flag indicating if last component was a match
  Arg[10]     : Int $i - reference component count
  Arg[11]     : Int $j - alternative component count
  Description : This function is called when two components match (i.e. have the
                same name.version in both assemblies). Depending on the state
                of the last component (match or nomatch), it extends aligned blocks
                or finishes the non-aligned block and creates a new aligned
                block.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_match {
	my (
		$R_chr,      $A_chr,      $match, $nomatch,
		$A_seg,      $last_A_seg, $R_seg, $last_R_seg,
		$match_flag, $i,          $j
	  )
	  = @_;

	# last component was a match
	if ($match_flag) {

		# adjust align block end
		if ( $match->{$R_chr} ) {

			my $c = scalar( @{ $match->{$R_chr} } ) - 1;

		  # if the gaps between this component and the last are different, start
		  # a new block
			if ( ( $A_seg->from_start - $match->{$R_chr}->[$c]->[1] ) !=
				( $R_seg->from_start - $match->{$R_chr}->[$c]->[4] ) )
			{

				$support->log(
					"Gap size mismatch at A:$A_chr:"
					  . $match->{$R_chr}->[$c]->[1] . '-'
					  . $A_seg->from_start
					  . ", R:$R_chr:"
					  . $match->{$R_chr}->[$c]->[4] . '-'
					  . $R_seg->from_start . "\n",
					2
				);

				# finish the last align block
				$match->{$R_chr}->[$c]->[1] = $last_A_seg->from_end;
				$match->{$R_chr}->[$c]->[2] = $j - $match->{$R_chr}->[$c]->[2];
				$match->{$R_chr}->[$c]->[4] = $last_R_seg->from_end;
				$match->{$R_chr}->[$c]->[5] = $i - $match->{$R_chr}->[$c]->[5];

				# start a new align block
				push @{ $match->{$R_chr} },
				  [
					$A_seg->from_start, $A_seg->from_end,
					$j,                 $R_seg->from_start,
					$R_seg->from_end,   $i,
					$A_chr,
				  ];

				# adjust align block end
			}
			else {
				$match->{$R_chr}->[$c]->[1] = $A_seg->from_end;
				$match->{$R_chr}->[$c]->[4] = $R_seg->from_end;
			}
		}

		# last component was a non-match
	}
	else {

		# start a new align block
		push @{ $match->{$R_chr} },
		  [
			$A_seg->from_start, $A_seg->from_end, $j,
			$R_seg->from_start, $R_seg->from_end, $i,
			$A_chr,
		  ];

		# finish the last non-align block
		if ( $nomatch->{$R_chr} ) {
			my $c = scalar( @{ $nomatch->{$R_chr} } ) - 1;
			$nomatch->{$R_chr}->[$c]->[1] = $last_A_seg->from_end;
			$nomatch->{$R_chr}->[$c]->[2] = $j - $nomatch->{$R_chr}->[$c]->[2];
			$nomatch->{$R_chr}->[$c]->[4] = $last_R_seg->from_end;
			$nomatch->{$R_chr}->[$c]->[5] = $i - $nomatch->{$R_chr}->[$c]->[5];
		}
	}
}

=head2 found_nomatch

  Arg[1]      : String $R_chr - reference chromosome name
  Arg[2]      : String $A_chr - alternative chromosome name
  Arg[3]      : Hashref $match - datastructure to store aligned blocks
  Arg[4]      : Hashref $nomatch - datastructure to store non-aligned blocks
  Arg[5]      : Bio::EnsEMBL::ProjectionSegment $A_seg - current alternative
                segment
  Arg[6]      : Bio::EnsEMBL::ProjectionSegment $last_A_seg - last alternative
                segment
  Arg[7]      : Bio::EnsEMBL::ProjectionSegment $R_seg - current reference
                segment
  Arg[8]      : Bio::EnsEMBL::ProjectionSegment $last_R_seg - last reference
                segment
  Arg[9]      : Boolean $match_flag - flag indicating if last component was a match
  Arg[10]     : Int $i - reference component count
  Arg[11]     : Int $j - alternative component count
  Description : This function is called when two components don't match (either
                different name.version or length mismatch in the two
                assemblies). Depending on the state of the last component (nomatch
                or match), it extends non-aligned blocks or finishes the
                aligned block and creates a new non-aligned block.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_nomatch {
	my (
		$R_chr,      $A_chr,      $match, $nomatch,
		$A_seg,      $last_A_seg, $R_seg, $last_R_seg,
		$match_flag, $i,          $j
	  )
	  = @_;

	# last component was a match
	if ($match_flag) {
		print STDOUT "last component was a match $R_chr ".$A_seg->from_start." ".$A_seg->from_end." ".$j." ".
			$R_seg->from_start." ".$R_seg->from_end." ".$i." ".
			$A_chr."\n";
		# start a new non-align block
		push @{ $nomatch->{$R_chr} },
		  [
			$A_seg->from_start, $A_seg->from_end, $j,
			$R_seg->from_start, $R_seg->from_end, $i,
			$A_chr,
		  ];

		# finish the last align block
		if ( $nomatch->{$R_chr} ) {
			my $c = scalar( @{ $match->{$R_chr} || [] } ) - 1;
			$match->{$R_chr}->[$c]->[1] = $last_A_seg->from_end;
			$match->{$R_chr}->[$c]->[2] = $j - $match->{$R_chr}->[$c]->[2];
			$match->{$R_chr}->[$c]->[4] = $last_R_seg->from_end;
			$match->{$R_chr}->[$c]->[5] = $i - $match->{$R_chr}->[$c]->[5];
		}

		# last component was a non-match
	}
	else {
		print STDOUT "adjust non-align block end $nomatch->{$R_chr}\n";
		# adjust non-align block end
		if ( $nomatch->{$R_chr} ) {
			my $c = scalar( @{ $nomatch->{$R_chr} || [] } ) - 1;
			$nomatch->{$R_chr}->[$c]->[1] = $A_seg->from_end;
			$nomatch->{$R_chr}->[$c]->[4] = $R_seg->from_end;
		} else {
					push @{ $nomatch->{$R_chr} },
		  [
			$A_seg->from_start, $A_seg->from_end, $j,
			$R_seg->from_start, $R_seg->from_end, $i,
			$A_chr,
		  ];
		}
	}
}

