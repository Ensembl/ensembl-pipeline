#!/usr/local/ensembl/bin/perl -w

use strict;
use IO::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


my $dbhost;
my $dbname;
my $outdir;

&GetOptions(
	    'dbhost:s'        => \$dbhost,
	    'dbname:s'        => \$dbname,
	    'outdir:s'        => \$outdir,
	    );

unless ( $dbhost && $dbname && $outdir ){
    print  "Usage: $0 -dbname -dbhost -outdir\n";
    exit(0);
}


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => 'ensro',
					    '-dbname' => $dbname,
					    );

print STDERR "connected to $dbname : $dbhost\n";

my @chrs = @{$db->get_ChromosomeAdaptor->fetch_all};

print STDERR scalar(@chrs)." chromosomes found\n";

my $fh = new IO::File;
foreach my $chr ( @chrs ){
    my $chr_name = $chr->chr_name;
    my $command = 
      "calculate_slices.pl -chr_name $chr_name -outfile $outdir/$chr_name.list -dbname $dbname -dbhost $dbhost";
    print STDERR "command $command\n";
    $fh->open("| bsub -q acari -C0 -o $outdir/lsf-out-$chr_name");
    $fh->print($command);
    $fh->close;
}
