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
    --chromosomes, --chr=LIST           only process LIST chromosomes
    --altchromosomes, --altchr=LIST     supply alternative chromosome names (the two lists must agree)

Optional arguments:

    --multiple                          produce a "many to one" mapping (default: off)
                                        many regions on the reference assembly can be mapped to one region
                                        on the alternative assembly (e.g. in new clone overlaps)
    --altdbname=NAME                    alternative database NAME

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
and Michael Gray <mg13@sanger.ac.uk>

=head1 CONTACT

Please post comments/questions to Anacode
<anacode-people@sanger.ac.uk>

=cut

use strict;
use warnings;

# Supress warnings from stats - but I want to see them for now
# no warnings 'uninitialized';

use FindBin qw($Bin);
use lib "$Bin";

use AssemblyMapper::Support;

use Pod::Usage;
use Readonly;
use Switch;

use Bio::EnsEMBL::Utils::Exception qw( throw );
use Algorithm::Diff qw(sdiff);

# project() returns triplets [start,end,slice]
#           where start,end are in source slice coords
#           and slice is in requested coord system
#           (the triplets are blessed as Bio::EnsEMBL::ProjectionSegment objects too)
#
Readonly my $PROJ_START => 0;      # $projection_segment->from_start()
Readonly my $PROJ_END   => 1;      # $projection_segment->from_end()
Readonly my $PROJ_SLICE => 2;      # $projection_segment->to_Slice()

# sdiff returns triplets [modifier,old_ele,new_ele]
#
Readonly my $SDIFF_MOD_TYPE  => 0;
Readonly my $SDIFF_LEFT_ELE  => 1;
Readonly my $SDIFF_RIGHT_ELE => 2;

$| = 1;

my $support = AssemblyMapper::Support->new(
    extra_options => [
        'skipcomponents|skip_components=s',
        'exctype=s',
        'ref_start=i',
        'ref_end=i',
        'alt_start=i',
        'alt_end=i',
        'multiple!'
    ],
    );

unless ($support->parse_arguments(@_)) {
    warn $support->error if $support->error;
    pod2usage(1);
}

my $write_db = not $support->param('dry_run');
my $multiple = $support->param('multiple');

if ($multiple) {
    $support->log_error("--multiple mode disabled for now following May 2011 rewrite\n");
}

$support->connect_dbs;

#####
# create temporary tables for storing non-aligned blocks and their masks
#
if ($write_db) {
    $support->ref_dbc->do(
        qq(
        CREATE TABLE IF NOT EXISTS tmp_align (
          tmp_align_id int(10) unsigned NOT NULL auto_increment,
          align_session_id int(10) unsigned NOT NULL,
          alt_start int(10) UNSIGNED NOT NULL,
          alt_end int(10) UNSIGNED NOT NULL,
          ref_start int(10) UNSIGNED NOT NULL,
          ref_end int(10) UNSIGNED NOT NULL,

          PRIMARY KEY (tmp_align_id),
          KEY align_session_id_idx (align_session_id)
          ) ENGINE=InnoDB;
  )
    );

    $support->ref_dbc->do(
        qq(
        CREATE TABLE IF NOT EXISTS tmp_mask (
          tmp_mask_id int(10) unsigned NOT NULL auto_increment,
          tmp_align_id int(10) unsigned NOT NULL,
          alt_mask_start int(10) UNSIGNED,
          alt_mask_end int(10) UNSIGNED,
          ref_mask_start int(10) UNSIGNED,
          ref_mask_end int(10) UNSIGNED,

          PRIMARY KEY (tmp_mask_id),
          KEY tmp_align_idx (tmp_align_id)
          ) ENGINE=InnoDB;
  )
    );

}
else {
    $support->log("\nHere I might create 'tmp_align' and 'tmp_mask' tables, if you'd let me\n");
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
my $fmt3 = "%-12s%-12s%-15s%-12s%-12s%-12s\n";
my $fmt33 = "%-12s%-12s%-15s%-12s%-12s%-15s%-12s\n";
my $fmt4 = " %10.0f  %10.0f    %7.0f   %10.0f  %10.0f  %7.0f\n";
my $fmt_rmask = "                                    [%10.0f  %10.0f mask ]\n";
my $fmt_amask = "[%10.0f  %10.0f mask ]\n";
my $fmt44 = "%10.0f  %10.0f    %7.0f   %10.0f  %10.0f  %7.0f %7.0f\n";
my $fmt5 = "%-40s%10s\n";
my $fmt6 = "%-10s%-12s%-10s%-12s\n";

my $sth1 = $support->ref_dbc->prepare(
    qq{
    INSERT IGNORE INTO assembly (asm_seq_region_id, cmp_seq_region_id,
                                 asm_start, asm_end, cmp_start, cmp_end, ori)
    VALUES (?, ?, ?, ?, ?, ?, ? )
}
    );

my $sth2 = $support->ref_dbc->prepare(
    qq{
  INSERT INTO tmp_align
    (align_session_id, alt_start, alt_end, ref_start, ref_end)
    VALUES (?, ?, ?, ?, ?)
}
    );

my $sth_mask = $support->ref_dbc->prepare(
    qq{
  INSERT INTO tmp_mask
      (tmp_align_id, alt_mask_start, alt_mask_end, ref_mask_start, ref_mask_end)
      values (?, ?, ?, ?, ?)
}
    );

my $ok = $support->iterate_chromosomes(
    prev_stage =>     '00-new_align_session',
    this_stage =>     '10-align_identical',
    worker =>         \&do_align,
    );

sub do_align {
    my ($asp, $data) = @_;

    my $match   = {};
    my $nomatch = {};
    my $match_flag = 0;
    my %stats_chr = ( identical => 0, mismatch => 0, skipped => 0 );

    my $R_chr   = $asp->ref_chr;
    my $A_chr   = $asp->alt_chr;

    my $R_slice = $asp->ref_slice;
    my $A_slice = $asp->alt_slice;

    # we need to fetch the alternative slice from the reference db explicitely by
    # coord_system, since toplevel attribute is not set there
    my $cs_name = $A_slice->coord_system_name;
    my $A_slice_ref = $R_slice->adaptor->fetch_by_region($cs_name, $A_chr, undef, undef, undef, $support->param('altassembly'));

    my $alt_seq_region_id_in_ref_db = $A_slice_ref->get_seq_region_id();

    my $R_length = $R_slice->length;
    my $A_length = $A_slice->length;

    my $Ref_start = $asp->ref_start;
    my $Alt_start = $asp->alt_start;

    my $R_asm = $asp->ref_asm;
    my $A_asm = $asp->alt_asm;

    my $R_dba = $asp->ref_dba;

    # project() returns triplets [start,end,slice]
    #           where start,end are from start of source slice
    #           and slice is in requested coord system
    #           (the triplets are blessed as Bio::EnsEMBL::ProjectionSegment objects too)
    #
    # NB base for triplet start,end coords is 1 => start of slice
    # They are NOT in the coord system of the slice...

    my @A_components = @{ $A_slice->project($component_cs) };
    $R_dba->get_AssemblyMapperAdaptor()->delete_cache();
    my @R_components = @{ $R_slice->project($component_cs) };

    # ...which should get taken care of here.
    #
    if($Alt_start){
        map ($_->[$PROJ_START] += ($Alt_start-1) , @A_components);
        map ($_->[$PROJ_END]   += ($Alt_start-1) , @A_components);
    }

    if($Ref_start){
        map ($_->[$PROJ_START] += ($Ref_start-1) , @R_components);
        map ($_->[$PROJ_END]   += ($Ref_start-1) , @R_components);
    }

    my @assembly_diffs = sdiff(\@R_components,\@A_components,\&get_cmp_key);

    # loop over sdiff results
  DIFF: foreach my $diff (@assembly_diffs)  {

      # sdiff returns triplets [modifier,old_ele,new_ele]

      my $type  = $diff->[$SDIFF_MOD_TYPE];
      my $R_ele = $diff->[$SDIFF_LEFT_ELE];
      my $A_ele = $diff->[$SDIFF_RIGHT_ELE];

      my $R_exists = $R_ele ? 1 : 0;
      my $A_exists = $A_ele ? 1 : 0;

      my ($R_ele_slice, $A_ele_slice);
      my ($R_key, $A_key) = ('<UNDEF>', '<UNDEF>');

      if ($R_exists) {
          $R_ele_slice = $R_ele->to_Slice;
          $R_key = &get_cmp_key($R_ele,1);
      }
      if ($A_exists) {
          $A_ele_slice = $A_ele->to_Slice;
          $A_key = &get_cmp_key($A_ele,1);
      }

      my ($left,$right,$tag) = ('-','-','');

      switch ($type) {

          case '+' { $left  = '-';    $tag   = '>'; $right = $A_key; }
          case '-' { $left  = $R_key; $tag   = '<'; $right = '-';    }
          case 'u' { $left  = $R_key; $tag   = '';  $right = $A_key; }
          case 'c' { $left  = $R_key; $tag   = '|'; $right = $A_key; }

          else { throw("Not expecting type '$type'"); }
      }

    MATCH_TYPE: {

        if($type eq 'u') {

            found_match($R_chr,$A_chr,$R_ele,$A_ele,$R_exists,$A_exists,$match,$match_flag);
            $stats_chr{'identical'}++;
            $match_flag = 1;

        }else{

            my $R_coords;
            if ($R_exists) {
                $R_coords = {
                    R_ctg_start  => $R_ele_slice->start,
                    R_ctg_end    => $R_ele_slice->end,
                    R_ctg_strand => $R_ele_slice->strand,
                    R_chr_start  => $R_ele->from_start,
                    R_chr_end    => $R_ele->from_end,
                };
            }

            if($type eq 'c') {

                my ($R_acc, $R_sv) = split m{\.},$R_ele_slice->seq_region_name;
                my ($A_acc, $A_sv) = split m{\.},$A_ele_slice->seq_region_name;

                if ( $R_acc eq $A_acc && $R_sv eq $A_sv ) {

                    $stats_chr{$R_chr}->{'mismatch'}++;

                    my $A_coords = {
                        A_ctg_start  => $A_ele_slice->start,
                        A_ctg_end    => $A_ele_slice->end,
                        A_ctg_strand => $A_ele_slice->strand,
                        A_chr_start  => $A_ele->from_start,
                        A_chr_end    => $A_ele->from_end,
                    };

                    my $blocks = &parse_projections( $R_coords, $A_coords );
                    store_blocks( $blocks, $R_chr, $A_chr, $A_exists, $match, $nomatch, $multiple );

                    $match_flag = 2;
                    last MATCH_TYPE;
                }
                elsif ( $R_acc ne $A_acc ) {

                    # Project the reference clone to the chromosome cs in case it
                    # matches another part in the alt. chromosome
                    my $s = $R_ele_slice->seq_region_Slice;
                    $R_dba->get_AssemblyMapperAdaptor()->delete_cache();
                    my @A_chrs = @{ $s->project( "chromosome", $A_asm ) };

                    foreach my $A_chr_proj (@A_chrs) {

                        my ($A_chr_name) = $A_chr_proj->to_Slice->seq_region_name;
                        if ($A_chr_name eq $A_chr) {

                            $R_dba->get_AssemblyMapperAdaptor()->delete_cache();

                            my ($proj_chr) =
                                @{ $s->project_to_slice( $A_chr_proj->to_Slice ) };

                            my $A_coords = {
                                A_ctg_start  => $proj_chr->from_start,
                                A_ctg_end    => $proj_chr->from_end,
                                A_ctg_strand => $proj_chr->to_Slice->strand,
                                A_chr_start  => $proj_chr->to_Slice->start,
                                A_chr_end    => $proj_chr->to_Slice->end,
                            };

                            my $blocks = &parse_projections( $R_coords, $A_coords );
                            store_blocks( $blocks, $R_chr, $proj_chr->to_Slice->seq_region_name, $A_exists,
                                          $match, $nomatch, $multiple );

                            $match_flag = 2;
                            $support->log(
                                sprintf("%-40s\t%-10s %-40s\n", $left, $tag, $right)
                                );
                            $support->log(
                                sprintf("%-40s\t%-10s %-40s\n", "", "", &get_cmp_key( $proj_chr, 1 ))
                                );
                            next DIFF;
                        }
                    }
                }
            }
            elsif ( $type eq '-' ) {

                my ( $R_acc, $R_sv ) = split m{\.}, $R_ele_slice->seq_region_name;

                # Project the reference clone slice to the chromosome cs in case
                # it matches another part of the alt. chromosome
                my $s = $R_ele_slice->seq_region_Slice;
                $R_dba->get_AssemblyMapperAdaptor()->delete_cache();
                my @A_chrs = @{ $s->project( "chromosome", $A_asm ) };

                foreach my $A_chr_proj (@A_chrs) {

                    my ($A_chr_name) = $A_chr_proj->to_Slice->seq_region_name;

                    if ($A_chr_name eq $A_chr) {

                        $R_dba->get_AssemblyMapperAdaptor()->delete_cache();
                        my ($proj_chr) = @{ $s->project_to_slice( $A_chr_proj->to_Slice ) };

                        my $A_coords = {
                            A_ctg_start  => $proj_chr->from_start,
                            A_ctg_end    => $proj_chr->from_end,
                            A_ctg_strand => $proj_chr->to_Slice->strand,
                            A_chr_start  => $proj_chr->to_Slice->start,
                            A_chr_end    => $proj_chr->to_Slice->end,
                        };

                        my $blocks = &parse_projections( $R_coords, $A_coords );
                        store_blocks( $blocks, $R_chr, $proj_chr->to_Slice->seq_region_name, $A_exists,
                                      $match, $nomatch, $multiple );

                        $match_flag = 2;
                        $right = &get_cmp_key( $proj_chr, 1 );
                        last MATCH_TYPE;
                    }
                }
            }
            $match_flag = 0;
        }
      } # MATCH_TYPE

      $support->log(sprintf("%-40s\t%-10s %-40s\n",$left,$tag,$right));

  } # DIFF

    # coalesce directly aligned blocks to improve handling (below) of
    # out-of-sequence alt blocks, sequence gaps, and zero-length gaps on one side

    my $number_aligned_blocks = scalar( @{ $match->{$R_chr} || [] } );
    $support->log("\nFound $number_aligned_blocks directly aligned blocks before coalescing\n", 1);

    my $current = $match->{$R_chr}->[0];
    if ($current) {
        my $c = 1;

      COALESCE: while (my $next_block = $match->{$R_chr}->[$c]) {

            my $ref_abutt = ( ($next_block->{R_start} - $current->{R_end}) == 1 );
            my $alt_abutt = ( ($next_block->{A_start} - $current->{A_end}) == 1 );
            my $same_ori  = ( ($next_block->{ori} * $current->{ori}) > 0 );

            if ($ref_abutt and $alt_abutt and $same_ori) {

                my $c_block = get_block_loc_str($current);
                my $n_block = get_block_loc_str($next_block);
                $support->log( "Merging adjacent blocks: $c_block : $n_block\n", 2);

                $current->{A_end} = $next_block->{A_end};
                $current->{R_end} = $next_block->{R_end};
                $current->{A_count} += $next_block->{A_count};
                $current->{R_count} += $next_block->{R_count};

                splice @{$match->{$R_chr}}, $c, 1; # drop $next_block from match array
                                                   # $c will now point at new $next_block
                next COALESCE;
            }

            $current = $next_block;
            $c++;
        }
    }

    my $c;

    # store directly aligned blocks in assembly table
    $number_aligned_blocks = scalar( @{ $match->{$R_chr} || [] } );
    if ($write_db) {

        $support->log(
            "Adding assembly entries for directly aligned blocks...\n", 1 );

        foreach $c ( 0 .. $number_aligned_blocks - 1 ) {
            $sth1->execute(
                $asp->ref_seq_region_id,
                $alt_seq_region_id_in_ref_db,
                $match->{$R_chr}->[$c]->{R_start},
                $match->{$R_chr}->[$c]->{R_end},
                $match->{$R_chr}->[$c]->{A_start},
                $match->{$R_chr}->[$c]->{A_end},
                $match->{$R_chr}->[$c]->{ori},
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

    # start/end coords of previous match, start with dummy zero
    my $A_dummy_prev_end = $Alt_start ? $Alt_start - 1 : 0;
    my $R_dummy_prev_end = $Ref_start ? $Ref_start - 1 : 0;
    my $prev_match = {
        A_start => $A_dummy_prev_end,
        A_end   => $A_dummy_prev_end,
        R_start => $R_dummy_prev_end,
        R_end   => $R_dummy_prev_end,
    };

    # dummy match at end of chromosome
    my $A_dummy_last_start = $Alt_start ? $Alt_start + $A_length : 1 + $A_length;
    my $R_dummy_last_start = $Ref_start ? $Ref_start + $R_length : 1 + $R_length;
    my $last = {
        A_start => $A_dummy_last_start,
        A_end   => $A_dummy_last_start,
        A_count => 0,
        R_start => $R_dummy_last_start,
        R_end   => $R_dummy_last_start,
        R_count => 0,
        A_name  => $A_chr
    };

    my @aug_match_list;
    @aug_match_list = @{ $match->{$R_chr} } if $match->{$R_chr};
    push @aug_match_list, $last;

    # copies of augmented match list, sorted by ref start and alt start respectively.
    my @matches_by_ref = sort { $a->{R_start} <=> $b->{R_start} } @aug_match_list;
    my @matches_by_alt = sort { $a->{A_start} <=> $b->{A_start} } @aug_match_list;

    # these cannot be initialised in-loop as we may wish to restart a block without clearing mask lists.
    my @ref_masks = ();
    my @alt_masks = ();

    my ($match_by_ref, $match_by_alt); # declare here or redo won't work as anticipated

  DIR_ALIGN_BLOCK: while ( $match_by_ref = shift @matches_by_ref, $match_by_alt = shift @matches_by_alt ) {

        unless ($match_by_ref and $match_by_alt) {
            my $what;
            $what = 'ref' unless $match_by_ref;
            $what = 'alt' unless $match_by_alt;
            die "Oops, reached the end and ran out of $what blocks";
        }

        my $a_start = $match_by_ref->{A_start};
        my $r_start = $match_by_ref->{R_start};
        my $a_chr   = $match_by_ref->{A_name}; # NB not the same as $A_chr !!

        my $sa_end  = $prev_match->{A_end};
        my $sr_end  = $prev_match->{R_end};

        my $ref_abutt = ( ($r_start - $sr_end) == 1 );
        my $alt_abutt = ( ($a_start - $sa_end) == 1 );

        if ( $ref_abutt and $alt_abutt ) {
            $support->log( "Both ref and alt blocks abutt, skipping (ref:$sr_end-$r_start)\n", 1);
            next DIR_ALIGN_BLOCK;
        }

        my $gap = get_gap_loc_str($prev_match, $match_by_ref);
        my $a_block = get_block_loc_str($match_by_alt);

        # Short-cut pointer comparison: both elements come from same source array, @aug_match_list.
        my $in_step = ($match_by_ref == $match_by_alt);
        unless ($in_step) {

            my $alt_a_start = $match_by_alt->{A_start};
            if ($alt_a_start < $a_start) {

                $support->log("There's an extra alt match block in this gap\n", 1);
                $support->log("  gap: $gap\n", 1);
                $support->log("  alt: $a_block\n", 1);

                # Alt block from a later match occurs in this gap, so:

                # - queue the intruding alt block for masking and move onto the next
                push @alt_masks, $match_by_alt unless $match_by_alt == $last;
                $match_by_alt = shift @matches_by_alt;

                # - restart comparison to use the new alt block, and in case there's another alt insert
                redo DIR_ALIGN_BLOCK if $match_by_alt;

            } else {

                $support->log("There's an extra ref match block, extending gap\n", 1);
                $support->log("  gap: $gap\n", 1);
                $support->log("  alt: $a_block\n", 1);

                # Alt block from this match occured earlier, thus this ref block is an 'extra'
                # ocurring in the gap between alt blocks, so:

                # - queue the intruding ref block for masking and move onto the next
                push @ref_masks, $match_by_ref unless $match_by_ref == $last;
                $match_by_ref = shift @matches_by_ref;

                # - restart comparision to use the new ref block, and in case there's another ref insert
                redo DIR_ALIGN_BLOCK if $match_by_ref;
            }

        }

        if ( ($r_start - $sr_end) < 1 ) {
            # This should never get hit
            $support->log_warning( "We shouldn't be going backwards!\n" );
            next DIR_ALIGN_BLOCK;
        }

        my ($chr_r_start, $chr_r_end) = ($sr_end+1, $r_start-1);

        # We filter out the null case where both abutt above.
        if ( $ref_abutt or $alt_abutt ) {

            $support->log("Zero-length ref gap (ref:$chr_r_start-$chr_r_end)\n", 1) if $ref_abutt;
            $support->log("Zero-length alt gap (ref:$chr_r_start-$chr_r_end)\n", 1) if $alt_abutt;

            # We put both ref and alt blocks onto the mask lists and extend the gap to the next one
            push @alt_masks, $match_by_alt unless $match_by_alt == $last;
            push @ref_masks, $match_by_ref unless $match_by_ref == $last;

            $match_by_alt = shift @matches_by_alt;
            $match_by_ref = shift @matches_by_ref;

            redo DIR_ALIGN_BLOCK if ($match_by_alt and $match_by_ref);
        }

        if ( ($a_start - $sa_end) < 1 ) {
            # This should never get hit now
            $support->log_warning("Negative alt gap for ref: $chr_r_start - $chr_r_end\n");
        }

        # We should apply masks before seeing if there's sequence...

        my $ref_seq = $asp->ref_sa->fetch_by_region(
            'chromosome', $R_chr, $chr_r_start, $chr_r_end, undef, $R_asm)->seq;
        my $ref_count = $ref_seq =~ s/([atgc])/$1/ig;

        my ($chr_a_start, $chr_a_end) =    $sa_end+1 < $a_start-1 ?
            ($sa_end+1, $a_start-1) :
            ($a_start-1, $sa_end+1);

        my $alt_seq = $asp->alt_sa->fetch_by_region(
            'chromosome', $A_chr, $chr_a_start, $chr_a_end, undef, $A_asm)->seq;
        my $alt_count = $alt_seq =~ s/([atgc])/$1/ig;

        if ( $ref_count and $alt_count ) {

            # We have a pair of unmatched blocks, both of which have sequence. Save them!
            push @{ $nomatch->{$R_chr} },
            {
                A_start => $chr_a_start, # 0
                A_end   => $chr_a_end,   # 1
                A_count => 1,            # 2
                R_start => $chr_r_start, # 3
                R_end   => $chr_r_end,   # 4
                R_count => 1,            # 5
                A_name  => $a_chr,       # 6
                A_masks => [ @alt_masks ],
                R_masks => [ @ref_masks ],
            };

        } elsif ( $ref_count or $alt_count ) {

            $support->log("Unmatched block: gap in ref (ref: $gap)\n", 1) unless $ref_count;
            $support->log("Unmatched block: gap in alt (ref: $gap)\n", 1) unless $alt_count;

            unless ($match_by_alt and $match_by_ref) {
                $support->log("Ran out of sequence (may be new extension on one side)\n", 1);
                last DIR_ALIGN_BLOCK;
            }

            # We add masks, and extend the unmatched block to the next one
            push @alt_masks, $match_by_alt unless $match_by_alt == $last;
            push @ref_masks, $match_by_ref unless $match_by_ref == $last;

            $match_by_alt = shift @matches_by_alt;
            $match_by_ref = shift @matches_by_ref;

            redo DIR_ALIGN_BLOCK if ($match_by_alt and $match_by_ref);

            $support->log_warning("Reached end of list with gap on one side of unmatched block\n");

        } else {
            $support->log("Skipping unmatched block: gaps in both ref and alt (ref: $gap)\n", 1);
        }

    } continue { # DIR_ALIGN_BLOCK

        $prev_match->{A_end} = $match_by_ref->{A_end};
        $prev_match->{R_end} = $match_by_ref->{R_end};

        # clear the mask lists before starting on next gap
        @ref_masks = ();
        @alt_masks = ();
    }

    die "Oops, reached the end and still had ref blocks" if @matches_by_ref;
    die "Oops, reached the end and still had alt blocks" if @matches_by_alt;

    # filter single assembly inserts from non-aligned blocks (i.e. cases where
    # a block has components only in one assembly, not in the other) - there is
    # nothing to be done with them
    @{ $nomatch->{$R_chr} } =
        grep { $_->{A_count} > 0 and $_->{R_count} > 0 } @{ $nomatch->{$R_chr} }
    if ( $nomatch->{$R_chr} );



    # store non-aligned blocks in tmp_align table
    my $number_nonaligned_blocks = scalar( @{ $nomatch->{$R_chr} || [] } );
    if ($write_db) {

        if ( $nomatch->{$R_chr} ) {

            $support->log( "Storing non-aligned blocks & masks in tmp_align & tmp_mask tables...\n",
                           1 );

            my $n_masks = 0;

            foreach $c ( 0 .. $number_nonaligned_blocks - 1 ) {

                my $block = $nomatch->{$R_chr}->[$c];

                $sth2->execute(
                    $asp->session_id,
                    $block->{A_start},
                    $block->{A_end},
                    $block->{R_start},
                    $block->{R_end},
                    );

                my $tmp_align_id = $sth2->{mysql_insertid}; # retrieve auto-incr insert id

                foreach my $mask ( @{$block->{A_masks}} )
                {
                    $sth_mask->execute($tmp_align_id, $mask->{A_start}, $mask->{A_end}, undef, undef);
                    $n_masks++;
                }

                foreach my $mask ( @{$block->{R_masks}} )
                {
                    $sth_mask->execute($tmp_align_id, undef, undef, $mask->{R_start}, $mask->{R_end});
                    $n_masks++;
                }
            }

            $support->log( "Done inserting $number_nonaligned_blocks tmp_align entries,\n", 1 );
            $support->log( " and $n_masks tmp_mask entries.\n", 1 );
        }
    }
    else {
        my $n_masks = 0;
        foreach my $block ( @{ $nomatch->{$R_chr} || [] } ) {
            map { $n_masks++ } @{$block->{A_masks}};
            map { $n_masks++ } @{$block->{R_masks}};
        }
        $support->log( "\n" );
        $support->log( "Here I would insert $number_nonaligned_blocks rows into 'tmp_align' table,\n", 1 );
        $support->log( "and $n_masks rows into 'tmp_mask' table, if you'd let me\n", 1 );
    }

    # stats for this chromosome
    $stats_chr{'A_only'} =
        scalar(@A_components) - $stats_chr{'identical'} - $stats_chr{'mismatch'};
    $stats_chr{'R_only'} =
        scalar(@R_components) - $stats_chr{'identical'} - $stats_chr{'mismatch'};
    for ( $c = 0 ; $c < scalar( @{ $match->{$R_chr} || [] } ) ; $c++ ) {
        $stats_chr{'A_matchlength'} +=
            $match->{$R_chr}->[$c]->{A_end} - $match->{$R_chr}->[$c]->{A_start};
        $stats_chr{'R_matchlength'} +=
            $match->{$R_chr}->[$c]->{R_end} - $match->{$R_chr}->[$c]->{R_start};
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
            sprintf( $fmt33,
                     qw(ALT_START ALT_END ALT_COMPONENTS REF_START REF_END REF_COMPONENTS ORIENTATION)
            ),
            2
            );
        $support->log( ( '-' x 80 ) . "\n", 2 );

        for ( $c = 0 ; $c < scalar( @{ $match->{$R_chr} } ) ; $c++ ) {

            $support->log( sprintf( $fmt44, @{$match->{$R_chr}->[$c]}{qw(A_start A_end A_count R_start R_end R_count ori)} ), 2 );

            # sanity check: aligned region pairs must have same length
            my $e_len =
                $match->{$R_chr}->[$c]->{A_end} - $match->{$R_chr}->[$c]->{A_start} + 1;
            my $v_len =
                $match->{$R_chr}->[$c]->{R_end} - $match->{$R_chr}->[$c]->{R_start} + 1;

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

            my $nmb = $nomatch->{$R_chr}->[$c];
            $support->log(sprintf( $fmt4, @{ $nmb }{qw(A_start A_end A_count R_start R_end R_count A_name)} ), 2);

            if ($nmb->{R_masks}) {
                foreach my $mask (@{$nmb->{R_masks}}) {
                    $support->log(sprintf( $fmt_rmask, $mask->{R_start}, $mask->{R_end} ), 2);
                }
            }

            if ($nmb->{A_masks}) {
                foreach my $mask (@{$nmb->{A_masks}}) {
                    $support->log(sprintf( $fmt_amask, $mask->{A_start}, $mask->{A_end} ), 2);
                }
            }

            # find longest non-aligned block
            my $A_length = $nmb->{A_end} - $nmb->{A_start} + 1;
            my $R_length = $nmb->{R_end} - $nmb->{R_start} + 1;
            push @block_length, [ $A_chr, $A_length, $R_chr, $R_length ];

        }
    }

    $support->log_stamped( "\nDone with chromosome $R_chr.\n", 1 );

    return 1;
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

sub get_block_loc_str {
    my ($block) = @_;
    return sprintf('R:%d-%d,A:%d-%d',
                   $block->{R_start}, $block->{R_end},
                   $block->{A_start}, $block->{A_end},
        );
}

sub get_gap_loc_str {
    my ($prev, $block) = @_;
    return sprintf('R:%d-%d,A:%d-%d',
                   $prev->{R_end}+1, $block->{R_start}-1,
                   $prev->{A_end}+1, $block->{A_start}-1,
        );
}

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
    Arg[8]      : Boolean $match_flag - flag indicating if last component was a match

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
        $i, $j, $match, $match_flag
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
            my $A_gap = $A_start-$match->{$R_chr}->[$c]->{A_end};
            my $R_gap = $R_start-$match->{$R_chr}->[$c]->{R_end};

            if($A_gap == $R_gap){
                $support->log("adjust align block end\n");

                $match->{$R_chr}->[$c]->{A_end} = $A_end;
                $match->{$R_chr}->[$c]->{A_count}++ if $j;
                $match->{$R_chr}->[$c]->{R_end} = $R_end;
                $match->{$R_chr}->[$c]->{R_count}++ if $i;
            } else {
                $support->log("start a new align block (because of gap)\n");
                push @{ $match->{$R_chr} },
                {
                    A_start => $A_start, # 0
                    A_end   => $A_end,   # 1
                    A_count => $j,       # 2
                    R_start => $R_start, # 3
                    R_end   => $R_end,   # 4
                    R_count => $i,       # 5
                    A_name  => $A_chr,   # 6
                    ori     => 1,        # 7 - why always 1 ?
                };
            }
        }

    }
    # last component was a non-match
    else {

        # start a new align block
        $support->log("start a new align block\n");
        push @{ $match->{$R_chr} },
        {
            A_start => $A_start, # 0
            A_end   => $A_end,   # 1
            A_count => $j,       # 2
            R_start => $R_start, # 3
            R_end   => $R_end,   # 4
            R_count => $i,       # 5
            A_name  => $A_chr,   # 6
            ori     => 1,        # 7 - why always 1 ?
        };
    }
}

sub parse_projections {
    my ( $R, $A ) = @_;
    my (      $R_ctg_start, $R_ctg_end, $R_chr_start, $R_chr_end, $R_ctg_strand ) =
        @$R{qw(R_ctg_start   R_ctg_end   R_chr_start   R_chr_end   R_ctg_strand)};
    my (      $A_ctg_start, $A_ctg_end, $A_chr_start, $A_chr_end, $A_ctg_strand ) =
        @$A{qw(A_ctg_start   A_ctg_end   A_chr_start   A_chr_end   A_ctg_strand)};
    my $blocks = [];
    my $ori    = $R_ctg_strand * $A_ctg_strand;

    # check that contig slices overlap
    if ( $A_ctg_end < $R_ctg_start || $R_ctg_end < $A_ctg_start ) {
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
        {
            R_start => $R_chr_s,
            R_end   => $R_chr_e,
            A_start => $A_chr_s,
            A_end   => $A_chr_e,
            tag     => $start_offset < 0 ? -1 : 1,
            ori     => 0,
        };
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
    push @$blocks,
    {
        R_start => $R_chr_s,
        R_end   => $R_chr_e,
        A_start => $A_chr_s,
        A_end   => $A_chr_e,
        tag     => 0,
        ori     => $ori
    };

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
        {
            R_start => $R_chr_s,
            R_end   => $R_chr_e,
            A_start => $A_chr_s,
            A_end   => $A_chr_e,
            tag     => $end_offset > 0 ? -1 : 1,
            ori     => 0,
        };
    }
    return $blocks;
}

# $A_exists may not be necessary, but for fidelity with previous version...
#
sub store_blocks {
    my ( $blocks, $R_chr, $A_name, $A_exists, $match, $nomatch, $multiple ) = @_;

    foreach my $block (@$blocks) {
        my (          $R_s,   $R_e, $A_s,   $A_e, $tag, $ori ) =
            @$block{qw(R_start R_end A_start A_end tag ori)};
        if ($ori) {
            $support->log("start a new align block\n");
            push @{ $match->{$R_chr} },
            {
                A_start => $A_s,
                A_end   => $A_e,
                A_count => $A_exists,
                R_start => $R_s,
                R_end   => $R_e,
                R_count => 1,
                A_name  => $A_name,
                ori     => $ori,
            };
        }
        elsif($multiple && $tag == -1) {
            $support->log("start a new overlap non-align block\n");
            push @{ $nomatch->{$R_chr} },
            {
                A_start => $A_s,
                A_end   => $A_e,
                A_count => 1,
                R_start => $R_s,
                R_end   => $R_e,
                R_count => 1,
                A_name  => $A_name
            };
        }
    }
}

