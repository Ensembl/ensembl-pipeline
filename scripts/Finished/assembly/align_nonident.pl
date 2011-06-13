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
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use AssemblyMapper::BlastzAligner;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

$support->param('verbose', 1);
$support->param('interactive', 0);

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
my ($block_sql, $block_sth);

# first set connection parameters for alternative db
# both databases have to be on the same host, so we don't need to configure
# them separately
for my $prm ( qw(host port user pass dbname) ) {
    $support->param("alt$prm", $support->param($prm)) unless ($support->param("alt$prm"));
}

# reference database
my $R_dba = $support->get_database('ensembl', '');
my $R_pipe_dba = &get_pipe_db($R_dba);
my $R_dbc = $R_dba->dbc;
my $R_asm = $support->param('assembly');

# database containing the alternative assembly
my $A_dba = $support->get_database('ensembl', 'alt');
my $A_pipe_dba = &get_pipe_db($A_dba);
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


$block_sql = qq(SELECT * FROM tmp_align);
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
    $block_sql .= ' WHERE '.join(' AND ', @where);
}

$block_sth = $R_dbc->prepare($block_sql);
$block_sth->execute;

# Pre-prepare per-block mask queries
my $ref_mask_sql = qq(SELECT ref_mask_start AS mask_start, ref_mask_end AS mask_end
                        FROM tmp_mask
                       WHERE tmp_align_id = ? AND ref_mask_start IS NOT NULL);
my $ref_mask_sth = $R_dbc->prepare($ref_mask_sql);

my $alt_mask_sql = qq(SELECT alt_mask_start AS mask_start, alt_mask_end AS mask_end
                        FROM tmp_mask
                       WHERE tmp_align_id = ? AND alt_mask_start IS NOT NULL);
my $alt_mask_sth = $R_dbc->prepare($alt_mask_sql);

BLOCK: while (my $row = $block_sth->fetchrow_hashref) {

  my $id = $row->{'tmp_align_id'};

  # create BlastzAligner object
    my $aligner = AssemblyMapper::BlastzAligner->new(-SUPPORT => $support);

    # create tmpdir to store input and output
    $aligner->create_tempdir($support->param('tmpdir'));


  $aligner->id($id);
  $aligner->seq_region_name($row->{'ref_seq_region_name'});

  $support->log_stamped("Block with tmp_align_id = $id\n", 1);

  $ref_mask_sth->execute($id);
  my $ref_masks = $ref_mask_sth->fetchall_arrayref({}); # get array of hashrefs
  my $n_ref_masks = scalar(@$ref_masks);

  $alt_mask_sth->execute($id);
  my $alt_masks = $alt_mask_sth->fetchall_arrayref({}); # get array of hashrefs
  my $n_alt_masks = scalar(@$alt_masks);

  $support->log("Fetched $n_ref_masks ref masks for this block\n", 2);
  $support->log("Fetched $n_alt_masks alt masks for this block\n", 2);

   # write sequences to file
  my $A_basename = "alt_seq.$id";
  my $R_basename = "ref_seq.$id";

  $support->log("Writing sequences to fasta...\n", 2);

  # This is needed otherwise will get a sequence of N's for the ref slice
  ($R_pipe_dba ? $R_pipe_dba : $R_dba)->get_AssemblyMapperAdaptor()->delete_cache();

    my $R_slice;
  if ($R_pipe_dba) {
    eval {
    $R_slice = $R_pipe_dba->get_SliceAdaptor->fetch_by_region(
      'chromosome',
      $row->{'ref_seq_region_name'},
      $row->{'ref_start'},
      $row->{'ref_end'},
      1,
      $support->param('assembly'),
    );
    };
  }
  $R_slice = $R_dba->get_SliceAdaptor->fetch_by_region(
      'chromosome',
      $row->{'ref_seq_region_name'},
      $row->{'ref_start'},
      $row->{'ref_end'},
      1,
      $support->param('assembly'),
   ) unless $R_slice;

  $aligner->write_sequence(
      $R_slice,
      $support->param('assembly'),
      $R_basename,
      undef,
      $ref_masks,
  );

  ($A_pipe_dba ? $A_pipe_dba : $A_dba)->get_AssemblyMapperAdaptor()->delete_cache();


  my $A_slice;
  if($A_pipe_dba){
      eval {
    $A_slice = $A_pipe_dba->get_SliceAdaptor->fetch_by_region(
      'chromosome',
      $row->{'alt_seq_region_name'},
      $row->{'alt_start'},
      $row->{'alt_end'},
      1,
      $support->param('altassembly'),
    );
    };
  }
  $A_slice = $A_dba->get_SliceAdaptor->fetch_by_region(
      'chromosome',
      $row->{'alt_seq_region_name'},
      $row->{'alt_start'},
      $row->{'alt_end'},
      1,
      $support->param('altassembly'),
   ) unless $A_slice;

   $aligner->write_sequence(
      $A_slice,
      $support->param('altassembly'),
      $A_basename,
      undef,
      $alt_masks,
  );


  $support->log("Done.\n", 2);

  # skip unmasked ref/alt sequences longer than 1.1MB
  # This will avoid everlasting alignment...
   $support->log("Checking sequences...\n", 2);
   if($aligner->bad_sequences($A_basename, $R_basename)){
        $support->log_warning("Skip block $id (not soft-masked and too long)...\n", 2);
        next BLOCK;
   }

  # align using lastz
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
$support->log("\nDon't forget to drop the tmp_align & tmp_mask tables when all is done!\n\n");

# finish logfile
$support->finish_log;

sub get_pipe_db {
    my ($dba) = @_;
    my $metakey = 'pipeline_db_head';
    my ($opt_str) = @{ $dba->get_MetaContainer()->list_value_by_key($metakey) };

    return undef unless $opt_str;

    my %anycase_options = (
        eval $opt_str,
    );
    if ($@) {
        throw("Error evaluating '$opt_str' : $@");
    }
    my %uppercased_options = ();
    while( my ($k,$v) = each %anycase_options) {
        $uppercased_options{uc($k)} = $v;
    }

    return Bio::EnsEMBL::DBSQL::DBAdaptor->new(%uppercased_options);
}


