#!/software/bin/perl

=head1 NAME

align_nonident_regions.pl - create whole genome alignment between two closely
related assemblies for non-identical regions

=head1 SYNOPSIS

align_nonident_regions.pl [arguments]

Required arguments:

    --dbname, db_name=NAME              database name NAME
    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --assembly=ASSEMBLY                 assembly version ASSEMBLY
    --altdbname=NAME                    alternative database NAME
    --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional arguments:

    --chromosomes, --chr=LIST           only process LIST chromosomes
    --bindir=DIR                        look for program binaries in DIR
    --tmpfir=DIR                        use DIR for temporary files (useful for
                                        re-runs after failure)

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

This script is part of a series of scripts to create a mapping between two
assemblies. It assembles the chromosome coordinate systems of two different
assemblies of a genome by creating a whole genome alignment between the two.

The process assumes that the two assemblies are reasonably similar, i.e. there
are no major rearrangements or clones moved from one chromosome to another.

See "Related files" below for an overview of the whole process.

This particular script creates a whole genome alignment between two closely
related assemblies for non-identical regions. These regions are identified by
another script (align_by_clone_identity.pl) and stored in a temporary database
table (tmp_align).

Alignments are calculated by this algorithm:

    1. fetch region from tmp_align
    2. write soft-masked sequences to temporary files
    3. align using lastz
    4. filter best hits (for query sequences, i.e. alternative regions) using
       axtBest
    5. parse lastz output to create blocks of exact matches only
    6. remove overlapping target (reference) alignments
    7. write alignments to assembly table

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
	unshift(@INC, "$SERVERROOT/bioperl-0.7.2");
    unshift(@INC, "$SERVERROOT/bioperl-1.2.3-patched");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use AssemblyMapper::BlastzAligner;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'assembly=s',
    'altdbname=s',
    'altassembly=s',
    'bindir=s',
    'tmpdir=s',
    'chromosomes|chr=s@',
    'altchromosomes|altchr=s@',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'altdbname',
    'altassembly',
    'bindir',
    'tmpdir',
    'chromosomes',
    'altchromosomes',
);

if ($support->param('help') or $support->error) {
  warn $support->error if $support->error;
  pod2usage(1);
}

$support->param('verbose', 1);
$support->param('interactive', 0);




# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
    'assembly',
    'altassembly'
);

#####
# connect to database and get adaptors
#
my ($dba, $dbh, $sql, $sth);

# first set connection parameters for alternative db
# both databases have to be on the same host, so we don't need to configure
# them separately
for my $prm ( qw(host port user pass dbname) ) {
    $support->param("alt$prm", $support->param($prm)) unless ($support->param("alt$prm"));
}

# reference database
my $R_dba = $support->get_database('ensembl', '');
my $R_dbh = $R_dba->dbc->db_handle;
my $R_sa = $R_dba->get_SliceAdaptor;
my $R_asm = $support->param('assembly');

# database containing the alternative assembly
my $A_dba = $support->get_database('core', 'alt');
my $A_sa = $A_dba->get_SliceAdaptor;
my $A_asm = $support->param('altassembly');


# loop over non-aligned regions in tmp_align table
$support->log_stamped("Looping over non-aligned blocks...\n");

# parse the strings into lists:
$support->comma_to_list('chromosomes', 'altchromosomes');

my @R_chr_list = $support->param('chromosomes');
if(! scalar(@R_chr_list)) {
    @R_chr_list = $support->sort_chromosomes;

    if(scalar($support->param('altchromosomes')) ) {
        die "AltChromosomes list is defined while Chromosomes list is not!";
    }
}

my @A_chr_list = $support->param('altchromosomes');
if(! scalar(@A_chr_list) ) {
    @A_chr_list = @R_chr_list;
} elsif(scalar(@R_chr_list) != scalar(@A_chr_list)) {
    die "Chromosome lists do not match by length";
}


$sql = qq(SELECT * FROM tmp_align);
my @where = ();
if(@R_chr_list) {
  my $chr_string = join("', '", @R_chr_list);
  push @where, "ref_seq_region_name IN ('$chr_string')";
}
if(@A_chr_list) {
  my $altchr_string = join("', '", @A_chr_list);
  push @where, "alt_seq_region_name IN ('$altchr_string')";
}
if(scalar(@where)) {
    $sql .= ' WHERE '.join(' AND ', @where);
}
$sth = $R_dbh->prepare($sql);
$sth->execute;

while (my $row = $sth->fetchrow_hashref) {

  # create BlastzAligner object
	my $aligner = AssemblyMapper::BlastzAligner->new(-SUPPORT => $support);

	# create tmpdir to store input and output
	$aligner->create_tempdir($support->param('tmpdir'));

  my $id = $row->{'tmp_align_id'};
  $aligner->id($id);
  $aligner->seq_region_name($row->{'ref_seq_region_name'});

  $support->log_stamped("Block with tmp_align_id = $id\n", 1);

  my $A_slice = $A_sa->fetch_by_region(
      'chromosome',
      $row->{'alt_seq_region_name'},
      $row->{'alt_start'},
      $row->{'alt_end'},
      1,
      $support->param('altassembly'),
  );

  my $R_slice = $R_sa->fetch_by_region(
      'chromosome',
      $row->{'ref_seq_region_name'},
      $row->{'ref_start'},
      $row->{'ref_end'},
      1,
      $support->param('assembly'),
  );

  # write sequences to file
  my $A_basename = "alt_seq.$id";
  my $R_basename = "ref_seq.$id";

  $support->log("Writing sequences to fasta...\n", 2);

  $aligner->write_sequence(
      $A_slice,
      $support->param('altassembly'),
      $A_basename
  );

  $aligner->write_sequence(
      $R_slice,
      $support->param('assembly'),
      $R_basename
  );

  $support->log("Done.\n", 2);

  # align using blastz
  $support->log("Running lastz...\n", 2);
  $aligner->run_lastz($A_basename, $R_basename);
  $support->log("Done.\n", 2);


  # find best alignment with axtBest
  $support->log("Finding best alignment with axtBest...\n", 2);
  $aligner->find_best_alignment;
  $support->log("Done.\n", 2);

  # parse lastz output, and convert relative alignment coordinates to
  # chromosomal coords
  $support->log("Parsing lastz output...\n", 2);

  $aligner->parse_lastz_output;

  $aligner->adjust_coords(
      $row->{'alt_start'},
      $row->{'alt_end'},
      { $id => [ $row->{'ref_start'}, $row->{'ref_end'} ] }
  );

  $support->log("Done.\n", 2);

  # cleanup temp files
  $support->log("Cleaning up temp files...\n", 2);
  $aligner->cleanup_tmpfiles(
    "$A_basename.fa",
    "$R_basename.fa",
  );
  $support->log("Done.\n", 2);

  # log alignment stats
  $aligner->log_block_stats(2);

  $support->log_stamped("Done with block $id.\n", 1);


	# write alignments to assembly table
	$aligner->write_assembly($R_dba, [$row->{'ref_seq_region_name'}], [$row->{'alt_seq_region_name'}]);

	# overall stats
	$aligner->log_overall_stats;

	# cleanup
	$support->log_stamped("\nRemoving tmpdir...\n");
	$aligner->remove_tempdir;

} # while ($row = fetchrow...)



$support->log_stamped("Done.\n\n");

# remind to drop tmp_align
$support->log("\nDon't forget to drop the tmp_align table when all is done!\n\n");

# finish logfile
$support->finish_log;


