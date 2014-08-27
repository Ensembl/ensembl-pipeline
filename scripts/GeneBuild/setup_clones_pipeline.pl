#!/usr/env perl

use strict;
use warnings;

use Getopt::Long;

use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw/throw warning/;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use FindBin;

# Connection to the target DB
my $host   = '';
my $port   = '3306';
my $user   = '';
my $pass   = '';
my $dbname = '';

my $xmlfile;
my $fasta_file;
my $verbose = 0;
my $alignment_ln = 'bacend';
my $linking_ln = 'clone_linking';
my $create_analysis = 0;
my $prepare_fasta = 0;
my $upload_info = 0;
my $create_cfg = 0;
my $input_ids = 0;
my $cfg_dir;
my $clipping = 0;
my $exonerate_exe = '/software/ensembl/genebuild/usrlocalensemblbin/exonerate-0.9.0';
my $chunk_dir;
my $submit_toplevel;
my $include_non_reference = 0;
my $fastasplit_exe = '/software/ensembl/bin/fastasplit_random';


&GetOptions (
        'dbhost|host=s' => \$host,
        'dbport|port=s' => \$port,
        'dbuser|user=s' => \$user,
        'dbpass|pass=s' => \$pass,
        'dbname=s'      => \$dbname,
        'xmlfile=s'     => \$xmlfile,
        'fastafile=s'  => \$fasta_file,
        'config_dir=s'  => \$cfg_dir,
        'exonerate=s' => \$exonerate_exe,
        'exonerate_ln=s' => \$alignment_ln,
        'clone_ln=s' => \$linking_ln,
        'chunkdir=s' => \$chunk_dir,
        'submit_toplevel=s' => \$submit_toplevel,
        'include_non_reference!' => \$include_non_reference,
        'create_analysis!' => \$create_analysis,
        'prepare_fasta!' => \$prepare_fasta,
        'upload_info!' => \$upload_info,
        'create_cfg!' => \$create_cfg,
        'input_ids!' => \$input_ids,
        'clip!' => \$clipping,
        'verbose+'      => \$verbose,
        );

throw("You need database connection parameters to run the setup!") unless ($host and $port and $user and $dbname);

if (!($create_analysis | $prepare_fasta | $upload_info | $create_cfg | $input_ids)) {
    $create_analysis = 1;
    $prepare_fasta = 1;
    $upload_info = 1;
    $create_cfg = 1;
    $input_ids = 1;
}

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
        -host => $host,
        -port => $port,
        -user => $user,
        -pass => $pass,
        -dbname => $dbname,
        );

$| = 1 if ($verbose);
sanity_check($db, $fastasplit_exe, $chunk_dir);

my $clones_info = prepare_fasta_file($xmlfile, $fasta_file, $clipping, $prepare_fasta);
create_analyses($db, $alignment_ln, $linking_ln) if ($create_analysis);
upload_clones_information($db, $clones_info, $alignment_ln) if ($upload_info);
add_input_ids($db, $alignment_ln, $linking_ln) if ($input_ids);
#create_configs($cfg_dir) if ($create_cfg);


sub prepare_fasta_file {
    my ($xmlfile, $fasta_file, $clipping, $writing_fasta) = @_;

    my %clone_infos;
    my $single_entry = 0;
    my $insert_size;
    my $insert_stdev;
    my $trace_direction;
    my $trace_id;
    my $trace_name;
    my $clone_id;
    my $clone_lib;
    my $accession;
    my $clip_quality_left;
    my $clip_quality_right;
    my $clip_vector_left;
    my $clip_vector_right;
    my %seen_skipped;
    print STDOUT "Reading XML..." if ($verbose);
    open(my $fh, $xmlfile) || die("Could not open $xmlfile");
    while (my $line = <$fh>) {
        if ($line =~ /^\s+<trace>/){
            $single_entry = 1;
            $insert_size = undef;
            $insert_stdev = undef;
            $trace_direction = undef;
            $trace_id = undef;
            $accession = undef;
            $clip_quality_left = undef;
            $clip_quality_right = undef;
            $clip_vector_left = undef;
            $clip_vector_right = undef;
        }

        elsif ($single_entry == 1) {
            if ($line =~ /^\s+<ti>(\w+)<\/ti>/i) {
                $trace_id = $1;
            }
            elsif ($line =~ /^\s+<trace_name>(\w+)<\/trace_name>/i) {
                $trace_name = $1;
            }
            elsif ($line =~ /^\s+<clone_id>(\w+[-]\w+)<\/clone_id>/i) {
                $clone_id = $1;
            }
            elsif ($line =~ /^\s+<insert_size>(\d+)<\/insert_size>/i) {
                $insert_size = $1;
            }
            elsif ($line =~ /^\s+<insert_stdev>(\d+)<\/insert_stdev>/i) {
                $insert_stdev = $1;
            }
            elsif ($line =~ /^\s+<library_id>(\S*)<\/library_id>/i) {
                $clone_lib = $1;
            }
            elsif ($line =~ /^\s+<trace_end>([FR])\w*<\/trace_end>/i) {
                $trace_direction = $1;
            }
            elsif ($line =~ /^\s+<accession>(\w+)<\/accession>/i) {
                $accession = $1;
            }
            elsif ($line =~ /^\s+<clip_quality_left>(\d+)<\/clip_quality_left>/i) {
                $clip_quality_left = $1;
            }
            elsif ($line =~ /^\s+<clip_quality_right>(\d+)<\/clip_quality_right>/i) {
                $clip_quality_right = $1;
            }
            elsif ($line =~ /^\s+<clip_vector_left>(\d+)<\/clip_vector_left>/i) {
                $clip_vector_left = $1;
            }
            elsif ($line =~ /^\s+<clip_vector_right>(\d+)<\/clip_vector_right>/i) {
                $clip_vector_right = $1;
            }
            elsif ($line =~ /^\s+<\/trace>/) {
                if ($accession) {
                    print "Clone with accession $accession, skipping\n" if ($verbose > 1);
                    $seen_skipped{$accession} = 1;
                    next;
                }
                next if (!$trace_id or !$insert_size or !$insert_stdev or !$trace_direction or !$trace_name);
# >918936606:CH243-100A1:F:CH243:184000:36800:1098268172037
                $clone_infos{$trace_name} = [$trace_id, $clone_id, $trace_direction, $clone_lib, $insert_size, $insert_stdev, $trace_name, 0];
                push(@{$clone_infos{$trace_name}}, $clip_quality_left, $clip_quality_right, $clip_vector_left, $clip_vector_right) if ($clipping);
                $single_entry = 0;
            }
        }
    }
    close($fh) || die("Could not close $xmlfile");
    print STDOUT "Done\n" if ($verbose);
    if ($verbose > 1) {
        foreach my $key (keys %seen_skipped) {
            print STDOUT "$key\n";
        }
    }
    print STDOUT "Processing data..." if ($verbose);
    my $seq_count = 0;
    my $readfasta = Bio::SeqIO->new(-format => 'fasta', -file => $fasta_file);
    my $writefasta = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$fasta_file.'_formated');
    my $min_length = 1000000;
    while (my $seq = $readfasta->next_seq) {
        if (exists $clone_infos{$seq->id}) {
            $seq_count++;
            push(@{$clone_infos{$seq->id}}, $seq->length);
            if ($clipping) {
                my $clip_left = $clone_infos{$seq->id}->[8] > $clone_infos{$seq->id}->[10] ? $clone_infos{$seq->id}->[8] : $clone_infos{$seq->id}->[10];
                my $clip_right = $clone_infos{$seq->id}->[9] < $clone_infos{$seq->id}->[11] ? $clone_infos{$seq->id}->[9] : $clone_infos{$seq->id}->[11];
                my $sequence = substr($seq->seq, $clip_left-1, $clip_right-$clip_left+1);
                $seq->seq($sequence);
                $min_length = length($sequence) if (length($sequence) < $min_length);
            }
            if ($writing_fasta) {
                $seq->id($clone_infos{$seq->id}->[0]);
                $writefasta->write_seq($seq);
                print STDOUT 'Writing ', $seq->id, "\n" if ($verbose == 2);
            }
        }
        else {
            print STDOUT 'Skipping ', $seq->id, "\n" if ($verbose == 2);
        }
    }
    print STDOUT "Done\n" if ($verbose);
    warning("Your shortest sequence is $min_length bp long!");
    if ($writing_fasta) {
        my $nchunk = int($seq_count/100);
        warning("Getting $nchunk files") if ($verbose);
        system($fastasplit_exe.' '.$fasta_file."_formated $nchunk $chunk_dir");
        throw("$fastasplit_exe failed: $?") if ($?);
    }
    return \%clone_infos;
}

sub create_analyses {
    my ($db, $alignment_ln, $linking_ln) = @_;

    print STDOUT 'Creating analyses...' if ($verbose);
    my $analysis_adaptor = $db->get_AnalysisAdaptor;
    my $rule_adaptor = $db->get_RuleAdaptor;
    my $submit_exonerate_analysis = Bio::EnsEMBL::Pipeline::Analysis->new(
            -logic_name => 'Submit_'.$alignment_ln,
            -module => 'Dummy',
            -input_id_type => 'CLONECHUNK',
            );
    $analysis_adaptor->store($submit_exonerate_analysis);
    my $exonerate_analysis = Bio::EnsEMBL::Pipeline::Analysis->new(
            -logic_name => $alignment_ln,
            -module => 'ExonerateAlignFeature',
            -program => 'exonerate',
            -program_file => $exonerate_exe,
            -input_id_type => 'CLONECHUNK',
            );
    $analysis_adaptor->store($exonerate_analysis);
    my $exonerate_rule = Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $exonerate_analysis );
    $exonerate_rule->add_condition($submit_exonerate_analysis->logic_name);
    $rule_adaptor->store($exonerate_rule);
    my $accumulator = Bio::EnsEMBL::Pipeline::Analysis->new(
            -logic_name => $alignment_ln.'_wait',
            -module => 'Dummy',
            -input_id_type => 'ACCUMULATOR',
            );
    $analysis_adaptor->store($accumulator);
    my $accumulator_rule = Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $accumulator );
    $accumulator_rule->add_condition($exonerate_analysis->logic_name);
    $rule_adaptor->store($accumulator_rule);
    my $submit_clone_analysis = Bio::EnsEMBL::Pipeline::Analysis->new(
            -logic_name => 'Submit_'.$linking_ln,
            -module => 'CloneEndsLinking',
            -input_id_type => 'CLONE_LIBRARY',
            );
    $analysis_adaptor->store($submit_clone_analysis);
    my $clone_link_analysis = Bio::EnsEMBL::Pipeline::Analysis->new(
            -logic_name => $linking_ln,
            -module => 'CloneEndsLinking',
            -input_id_type => 'CLONE_LIBRARY',
            );
    $analysis_adaptor->store($clone_link_analysis);
    my $clone_link_rule = Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $clone_link_analysis );
    $clone_link_rule->add_condition($accumulator->logic_name);
    $clone_link_rule->add_condition($submit_clone_analysis->logic_name);
    $rule_adaptor->store($accumulator_rule);
    print STDOUT "Done\n" if ($verbose);
}

sub upload_clones_information {
    my ($db, $clone_infos, $alignment_ln) = @_;

    print STDOUT 'Uploading clones information...' if ($verbose);
    my $count = 0;
    my $analysis_adaptor = $db->get_AnalysisAdaptor;
    my $analysis = $analysis_adaptor->fetch_by_logic_name($alignment_ln);
    throw("Could not find analysis $alignment_ln") unless ($analysis);
    my $sth = $db->dbc->prepare('INSERT INTO clones (trace_id, clone_id, direction, library, insert_size, insert_stdev, trace_name, length, clip_qleft, clip_qright, clip_vleft, clip_vright, analysis_id) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
    foreach my $trace_name (keys %$clone_infos) {
        $sth->bind_param(1, $clone_infos->{$trace_name}->[0]);
        $sth->bind_param(2, $clone_infos->{$trace_name}->[1]);
        $sth->bind_param(3, $clone_infos->{$trace_name}->[2]);
        $sth->bind_param(4, $clone_infos->{$trace_name}->[3]);
        $sth->bind_param(5, $clone_infos->{$trace_name}->[4]);
        $sth->bind_param(6, $clone_infos->{$trace_name}->[5]);
        $sth->bind_param(7, $clone_infos->{$trace_name}->[6]);
        $sth->bind_param(8, $clone_infos->{$trace_name}->[7]);
        if (exists $clone_infos->{$trace_name}->[8]) {
            $sth->bind_param(9, $clone_infos->{$trace_name}->[8]);
            $sth->bind_param(10, $clone_infos->{$trace_name}->[9]);
            $sth->bind_param(11, $clone_infos->{$trace_name}->[10]);
            $sth->bind_param(12, $clone_infos->{$trace_name}->[11]);
        }
        else {
            $sth->bind_param(9, 1);
            $sth->bind_param(10, $clone_infos->{$trace_name}->[7]);
            $sth->bind_param(11, 1);
            $sth->bind_param(12, $clone_infos->{$trace_name}->[7]);

        }
        $sth->bind_param(13, $analysis->dbID);
        $sth->execute();
        $count++;
    }
    print STDOUT "Done\n" if ($verbose);
}

sub add_input_ids {
    my ($db, $alignment_ln, $linking_ln) = @_;

    my $analysis_adaptor = $db->get_AnalysisAdaptor;
    print STDOUT 'Adding input ids...' if ($verbose);
    my $submit_exonerate_analysis = $analysis_adaptor->fetch_by_logic_name('Submit_'.$alignment_ln);
    throw("Could not find analysis Submit_$alignment_ln in your database") unless ($submit_exonerate_analysis);
    my $iif_chunks = Bio::EnsEMBL::Pipeline::Utils::InputIDFactory->new(
            -db => $db,
            -file => 1,
            -logic_name => $submit_exonerate_analysis->logic_name,
            -input_id_type => $submit_exonerate_analysis->input_id_type,
            -dir => $chunk_dir,
            -verbose => $verbose,
            -include_non_reference => 0,
            );
    $iif_chunks->generate_input_ids;
    $iif_chunks->store_input_ids;
    my $submit_clone_analysis = $analysis_adaptor->fetch_by_logic_name('Submit_'.$linking_ln);
    throw("Could not find analysis Submit_$linking_ln in your database") unless ($submit_clone_analysis);
    my $sic_adaptor = $db->get_StateInfoContainer;
    my $sth_clones = $analysis_adaptor->dbc->prepare('SELECT DISTINCT(clone_id) FROM clones');
    $sth_clones->execute();
    foreach (@{$sth_clones->fetchall_arrayref}) {
        $sic_adaptor->store_input_id_analysis($_->[0], $submit_clone_analysis, 'dummy');
    }
    print STDOUT "Done\n" if ($verbose);
}

sub sanity_check {
    my ($db, $fastasplit_exe, $chunk_dir) = @_;

    my $error_msg = '';
    my $dbhandle = $db->dbc->db_handle();
    my $pipeline_tables = 0;
    my $clones_table = 0;
    foreach my $table ($dbhandle->tables(undef,undef,undef,'TABLE')) {
        $pipeline_tables++ if ($table =~ /`job`/ or
                $table =~ /`job_status`/ or
                $table =~ /`rule_goal`/ or
                $table =~ /`rule_conditions`/ or
                $table =~ /`input_id_type_analysis`/ or
                $table =~ /`input_id_analysis`/
                );
        $clones_table = 1 if ($table =~ /`clones`/);
    }

    my ($pipeline_sql) = $FindBin::Bin =~ /(.*ensembl-pipeline)/;
    $pipeline_sql .= '/sql';
    if ($pipeline_tables != 6) {
        if (-e $pipeline_sql) {
            my $sql_query;
            open(SQL, $pipeline_sql.'/table.sql') || throw("Cannot open table.sql in $pipeline_sql");
            while (<SQL>) {
                next if (/^\s*#/);
                chomp;
                $sql_query .= $_. ' ';
            }
            close(SQL) || throw("Could not close table.sql in $pipeline_sql");
            $sql_query =~ s/\s+/ /g;
            foreach my $query (split(';', $sql_query)) {
                next unless ($query =~ /\w/);
                my $sth_tables = $db->dbc->prepare($query);
                $sth_tables->execute();
            }
        }
        else {
            $error_msg .= "Can't find $pipeline_sql";
        }
    }
    if (!$clones_table) {
        if (-e $pipeline_sql) {
            my $sql_query;
            open(SQL, $pipeline_sql.'/clones.sql') || throw("Cannot open clones.sql in $pipeline_sql");
            while (<SQL>) {
                next if (/^\s*#/);
                chomp;
                $sql_query .= $_. ' ';
            }
            close(SQL) || throw("Could not close clones.sql in $pipeline_sql");
            $sql_query =~ s/\s+/ /g;
            foreach my $query (split(';', $sql_query)) {
                next unless ($query =~ /\w/);
                next if ($query =~ /^DROP/);
                my $sth_tables = $db->dbc->prepare($query);
                $sth_tables->execute();
            }
        }
        else {
            $error_msg .= "Can't find $pipeline_sql";
        }
    }
    if ($prepare_fasta) {
        $error_msg .= "Can't find $fastasplit_exe" unless (-x $fastasplit_exe);
        $error_msg .= "$chunk_dir does not exist!" unless (-d $chunk_dir);
    }
    throw($error_msg) if ($error_msg);
}
