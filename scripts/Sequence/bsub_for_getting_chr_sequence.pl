#!/usr/local/ensembl/bin/perl -w

use strict;
use IO::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $fh = new IO::File;

# connect to the database
my $dbhost;
my $dbuser = 'ensro';
my $dbname;
my $outdir;
my $coord_system;
my $coord_system_version;
my $dbpass;
my $dbport = 3306;
&GetOptions(
            'dbname:s'       => \$dbname,
            'dbhost:s'       => \$dbhost,
            'outdir:s'       => \$outdir,
            'dbuser:s' => \$dbuser,
            'dbport:s' => \$dbport,
            'dbpass:s' => \$dbpass,
            'coord_system:s' => \$coord_system,
            'coord_system_version:s' => \$coord_system_version,
	    );

unless( $dbhost &&  $dbname && $outdir && $coord_system){
    print STDERR "Captain Kirk! We don't have the database information, exiting the operation, beep!\n";
    print STDERR "Specify the directory (full path) where to send the chromosome files\n";
    print STDERR "Usage: $0 -dbname $dbname -dbhost $dbhost -outdir $outdir ".
      " -coord_system $coord_system -coord_system_version ".
        "$coord_system_version\n";
    exit(0);
}

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                           -host  => $dbhost,
                                           -user  => $dbuser,
                                           -dbname => $dbname,
                                           -port => $dbport,
                                           -pass => $dbpass,
                                          );


my $command = '/ecs2/work1/lec/code/briggsae/ensembl-pipeline/'.
  'scripts/Sequence/get_sequence_for_chr.pl';


my @chr_names = &get_chr_names($db, $coord_system, $coord_system_version);

foreach my $chr (@chr_names){
    my $command = "$command -seq_region_name $chr -outfile $outdir/$chr.fa ".
      "-mask -dust -softmask -dbname $dbname -dbhost $dbhost ".
        "-dbuser $dbuser -dbport $dbport -coord_system $coord_system";
    $fh->open("| bsub -q normal -m 'ecs2_hosts' -o $outdir/jobs/lsf-out-$chr");
    $fh->print($command);
    $fh->close;
}


sub get_chr_names{
    my $db = shift;
    my $cs = shift;
    my $cs_version = shift;

    my @chrs = @{$db->get_SliceAdaptor->fetch_all($cs, $cs_version)};
    my @chr_names;
    foreach my $chr(@chrs){
      push @chr_names, $chr->seq_region_name;
    }
    return @chr_names;
}
