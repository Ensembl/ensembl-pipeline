#!/usr/local/bin/perl


# InputIdManager.pl
#
# Script for managing input IDs in an EnsEMBL Pipeline database
#
# Creator: Simon Potter <scp@sanger.ac.uk>
# Date of creation: 24.11.2000
# Last modified : SCP 7.3.2001
#
# Copyright Simon Potter 2000
#
# You may distribute this code under the same terms as perl itself
#
# TODO: Some bits of this could be better written, maybe provide
# options to extract stats from InputIdAnalysis like which analyses
# need to be run etc...


=head1 NAME

InputIdManager.pl

=head1 SYNOPSIS

InputIdManager.pl [-options]

Options

    -insert    insert input IDs
    -delete    delete input IDs
    -dbname    DB to connect to
    -dbuser    DB username
    -dbpass    DB password
    -dbhost    DB host
    -class     inputid class (e.g. contig)
    -analysis  analysisid (e.g. Genscan)
    -file      file of list of inputids (can be - for stdin)
    -inputid   a single inputid to store
    -allanals  remove all analysis ids for given inputid
    -check     check format of inputids (e.g. XY123456 for a clone)
    -reject    check format of inputids and reject
    -help      show help
    -info      show document info with perldoc

=head1 DESCRIPTION

B<InputIdManager.pl> facilitates additions/deletions
of input IDs into an EnsEMBL pipeline InputIdAnalysis table

=head1 CONTACT

Simon Potter: scp@sanger.ac.uk

=head1 BUGS

Insert list of bugs here!

=cut


use strict;
use FileHandle;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Pipeline::DBSQL::Obj;

$Getopt::Long::autoabbrev = 0;


my ($file, $filename, $delete, $insert, $allanals);
my ($id, $anal, $analysis, $check, $reject, $help, $info);

# defaults

my $dbname = $ENV{'ENS_DBNAME'} || undef;
my $dbuser = $ENV{'ENS_USER'}   || 'ensadmin';
my $dbpass = $ENV{'ENS_PASS'}   || undef;
my $dbhost = $ENV{'ENS_HOST'}   || undef;
my $class  = 'contig';


unless (GetOptions(
    "insert"      => \$insert,
    "delete"      => \$delete,
    "allanals"    => \$allanals,
    "dbname=s"    => \$dbname,
    "dbuser=s"    => \$dbuser,
    "dbhost=s"    => \$dbhost,
    "dbpass=s"    => \$dbpass,
    "inputid=s"   => \$id,
    "class=s"     => \$class,
    "analysis=s"  => \$analysis,
    "file=s"      => \$filename,
    "check"       => \$check,
    "reject"      => \$reject,
    "help"        => \$help,
    "info"        => \$info
)) {
    usage();
    exit 2;
}

if ($info) {
    exec("perldoc $0");
}

if ($help) {
    usage();
    exit 0;
}

if (($delete && $insert) || ($filename && $id)) {
    usage();
    exit 2;
}

my $db = Bio::EnsEMBL::Pipeline::DBSQL::Obj->new(
    -user   => $dbuser,
    -dbname => $dbname,
    -host   => $dbhost,
    -pass   => $dbpass
);

# Oh bugger...
die "Unsuccessful DB connection to $dbname as user $dbuser" unless $db;

my $sic = $db->get_StateInfoContainer;
my $analA = $db->get_AnalysisAdaptor;

if ($filename) {
    $file = new FileHandle;
    if ($filename eq '-') {
	$file = \*STDIN;
    }
    else {
	open $file, "< $filename" or die "Can't find file $filename";
    }
}

# for checking the input id string conforms to spec.
# e.g. AB000381 for clone, AB000381.1.1.35863 for contig

my %IDformat = (
    'clone',  qw{^[A-Z]+\d+$},
    'contig', qw{^[A-Z]+\d+\.\d+\.\d+\.\d+$}
);

if ($analysis) {
    $anal = $analA->fetch_by_newest_logic_name($analysis)
     or die "Can't find analysis $analysis";
}

if ($file) {
    while (<$file>) {
	chomp;
	($id) = split;
	next unless $id =~ /\S/;
	print "$id\n";
	insert_delete();
    }
}
else {
    insert_delete();
}

close $file if ($file && $file ne '-');



sub insert_delete {
    if (defined $IDformat{$class} && $id !~ $IDformat{$class}) {
	if ($check) {
	    print "Warning: is $id really a valid $class ID?\n";
	}
	if ($reject) {
	    print "Error: $id not a valid $class ID\n";
	    exit 1;
	}
    }
    if ($insert) {
	print "Inserting $id into DB\n";
	my $res = $sic->store_inputId_class_analysis($id, $class, $anal);
	print STDERR "Warning: inputid $id not stored\n" unless $res == 1;
    } else {
	if ($allanals) {
	    print "Deleting $id from DB\n";
	    my $res = $sic->delete_inputId($id);
	    print STDERR "Deleted $res lines\n";
	}
	else {
	    print "Deleting $id from DB\n";
	    my $res = $sic->delete_inputId_analysis($id, $anal->dbID);
	    print STDERR "Deleted $res lines\n";
	}
    }
}



sub usage {
print <<EOF;

Usage: InputIdManager.pl [-insert | -delete] [-file ... | -inputid ...] [-options]
options are: -dbname    database to connect to
             -dbhost    host to connect to
             -dbuser    username
             -dbpass    password
             -file      file containing list of IDs
             -inputid   an inputid to do something with
             -analysis  analysis logic_name (e.g. SubmitContig)
             -class     class name (e.g. contig)
             -check     check format of inputids (e.g. XY123456 for a clone)
             -reject    check format of inputids and reject
             -allanals  remove all analysis ids for given inputid
             -help      show help
             -info      show document info with perldoc

EOF
}
