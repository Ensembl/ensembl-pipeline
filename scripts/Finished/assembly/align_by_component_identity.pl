#!/software/bin/perl


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

    1. Match components with same name and version directly using the sdiff method
       in Algorithm::Diff and create alignment blocks for these regions. Components
       can be tagged manually to be excluded from these direct matches by listing
       them in a file of components to skip (--skipcomponents argument). This can be
       useful to get better results in regions with major assembly differences.

       The result is stored in the assembly table as an assembly between the
       chromosomes of both genome assemblies.

    2. Store non-aligned blocks in a temporary table (tmp_align). They can
       later be aligned using lastz by align_nonident.pl.

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

$support->param( 'verbose',     1 );    # throw away all that garbage
$support->param( 'interactive', 0 );    # stop that garbage from coming up

my $write_db = not $support->param('dry_run');

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

my %stats_total;
my @block_length;

my $fmt1 = "%-40s%10.0f\n";
my $fmt2 = "%-40s%9.2f%%\n";
my $fmt3 = "%-12s%-12s%-12s%-12s%-12s%-9s\n";
my $fmt4 = "%10.0f  %10.0f    %7.0f   %10.0f  %10.0f  %7.0f\n";
my $fmt5 = "%-40s%10s\n";
my $fmt6 = "%-10s%-12s%-10s%-12s\n";

my $sth1 = $R_dbh->prepare(
	qq{
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
	my $match   = {};
	my $nomatch = {};
	my ($i,$j)  = (0,0);
	my ($left,$right,$type,$tag);
	my $match_flag = 0;
	my %stats_chr;

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
	$R_dba->get_AssemblyMapperAdaptor()->delete_cache();
	my @R_components = @{ $R_slice->project($component_cs) };

	if($A_start){
		map ($_->[0] += ($A_start-1) , @A_components);
		map ($_->[1] += ($A_start-1) , @A_components);
	}

	if($R_start){
		map ($_->[0] += ($R_start-1) , @R_components);
		map ($_->[1] += ($R_start-1) , @R_components);
	}

	my @assembly_diffs = sdiff(\@R_components,\@A_components,\&get_cmp_key);

	# loop over sdiff results
	DIFF: foreach my $diff (@assembly_diffs)  {
		$type = $diff->[0];
		($left,$right,$tag) = ('-','-','');
		($i,$j) = (0,0);
		if($type eq '+'){
			$right = &get_cmp_key($diff->[2],1);
			$tag = '>';
		}elsif($type eq '-'){
			$left = &get_cmp_key($diff->[1],1);
			$tag = '<';
		}elsif($type eq 'u'){
			$left = &get_cmp_key($diff->[1],1);
			$right = &get_cmp_key($diff->[2],1);
		}elsif($type eq 'c'){
			$left = &get_cmp_key($diff->[1],1);
			$right = &get_cmp_key($diff->[2],1);
			$tag = '|';
		}

		if($type eq 'u') {
			($i,$j) = (1,1);
			if( $skip{ $diff->[1]->to_Slice->seq_region_name } ) {
				found_nomatch($R_chr,$A_chr,$diff->[1],$diff->[2],$i,$j,$match, $nomatch,$match_flag);
				$stats_chr{'skipped'}++;
				$match_flag = 0;
			} else {
				found_match($R_chr,$A_chr,$diff->[1],$diff->[2],$i,$j,$match, $nomatch,$match_flag);
				$stats_chr{'identical'}++;
				$match_flag = 1;
			}
		}else{
			$i = 1 if $type eq '-';
			$j = 1 if $type eq '+';
			if($type eq 'c') {
				($i,$j) = (1,1);
				my ($R_acc, $R_sv) = split /\./,$diff->[1]->to_Slice->seq_region_name;
				my ($A_acc, $A_sv) = split /\./,$diff->[2]->to_Slice->seq_region_name;

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

			}
			found_nomatch($R_chr,$A_chr,$diff->[1],$diff->[2],$i,$j,$match, $nomatch,$match_flag);
			$match_flag = 0;
		}

		$support->log(sprintf("%-40s\t%-10s %-40s\n",$left,$tag,$right));

	}

	my $c;


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

		foreach $c ( 0 .. $number_aligned_blocks - 1 ) {
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
	for ( $c = 0 ; $c < scalar( @{ $match->{$R_chr} || [] } ) ; $c++ ) {
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

		for ( $c = 0 ; $c < scalar( @{ $match->{$R_chr} } ) ; $c++ ) {

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

		for ( $c = 0 ; $c < scalar( @{ $nomatch->{$R_chr} } ) ; $c++ ) {

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

sub get_cmp_key {
	my ($cmp, $flag) = @_;
	my $slice = $cmp->to_Slice;
	my $key = $slice->seq_region_name . ":". $slice->start . "-". $slice->end . ":". $slice->strand;
	$key.= ":".$cmp->from_start."-".$cmp->from_end.":".($cmp->from_end-$cmp->from_start+1) if $flag;
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
	my (
		$R_chr, $A_chr, $R_seg, $A_seg,
		$i, $j, $match, $nomatch, $match_flag
	  )
	  = @_;

	my $R_start = $R_seg->from_start;
	my $R_end = $R_seg->from_end;
	my $A_start = $A_seg->from_start;
	my $A_end = $A_seg->from_end;

	# last component was a match
	if ( $match_flag == 1 ) {

		# adjust align block end
		if ( $match->{$R_chr} ) {
			# check that there is no single gap in the block
			my $c = scalar( @{ $match->{$R_chr} } ) - 1;
			my $A_gap = $A_start-$match->{$R_chr}->[$c]->[1];
			my $R_gap = $R_start-$match->{$R_chr}->[$c]->[4];

			if($A_gap == $R_gap){
				$support->log("adjust align block end\n");

				$match->{$R_chr}->[$c]->[1] = $A_end;
				$match->{$R_chr}->[$c]->[2]++ if $j;
				$match->{$R_chr}->[$c]->[4] = $R_end;
				$match->{$R_chr}->[$c]->[5]++ if $i;
			} else {
				$support->log("start a new align block (because of gap)\n");
				push @{ $match->{$R_chr} },
				  [
					$A_start, $A_end, $j,
					$R_start, $R_end, $i,
					$A_chr,
				  ];
			}
		}

	}
	# last component was a non-match
	else {

		# start a new align block
		$support->log("start a new align block\n");
		push @{ $match->{$R_chr} },
		  [
			$A_start, $A_end, $j,
			$R_start, $R_end, $i,
			$A_chr,
		  ];
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
	my (
		$R_chr, $A_chr, $R_seg, $A_seg,
		$i, $j, $match, $nomatch, $match_flag
	  )
	  = @_;
	my ($R_start, $R_end, $A_start, $A_end) = (0,0,0,0);

	if($i) {
		$R_start = $R_seg->from_start;
		$R_end = $R_seg->from_end;
	}
	if($j) {
		$A_start = $A_seg->from_start;
		$A_end = $A_seg->from_end;
	}

	# last component was a match
	if ($match_flag) {
		# start a new non-align block
		$support->log("start a new non-align block\n");
		push @{ $nomatch->{$R_chr} },
		  [
			$A_start, $A_end, $j,
			$R_start, $R_end, $i,
			$A_chr,
		  ];
	}
	# last component was a non-match
	else {
		# adjust non-align block end
		if ( $nomatch->{$R_chr} ) {
			$support->log("adjust non-align block end\n");
			my $c = scalar( @{ $nomatch->{$R_chr} || [] } ) - 1;
			$nomatch->{$R_chr}->[$c]->[0] ||=  $A_start;
			$nomatch->{$R_chr}->[$c]->[1] = $A_end if $A_end;
			$nomatch->{$R_chr}->[$c]->[2]++ if $j;
			$nomatch->{$R_chr}->[$c]->[3] ||= $R_start;
			$nomatch->{$R_chr}->[$c]->[4] = $R_end if $R_end;
			$nomatch->{$R_chr}->[$c]->[5]++ if $i;
		} else {
			$support->log("start a new non-align block\n");
					push @{ $nomatch->{$R_chr} },
		  [
			$A_start, $A_end, $j,
			$R_start, $R_end, $i,
			$A_chr,
		  ];
		}
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


