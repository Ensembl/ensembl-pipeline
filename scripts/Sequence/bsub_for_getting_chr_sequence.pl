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

&GetOptions(
	    'dbname:s'       => \$dbname,
	    'dbhost:s'       => \$dbhost,
	    'outdir:s'       => \$outdir,
	    );

unless( $dbhost &&  $dbname && $outdir ){
    print STDERR "Captain Kirk! We don't have the database information, exiting the operation, beep!\n";
    print STDERR "Specify the directory (full path) where to send the chromosome files\n";
    print STDERR "Usage: $0 -dbname dbname -dbhost dbhost -outdir outdir\n";
    exit(0);
}

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
					   -host  => $dbhost,
					   -user  => $dbuser,
					   -dbname=> $dbname
					   );



my @chr_names = &get_chr_names($db);

foreach my $chr (@chr_names){
    my $command = "get_sequence_for_chr.pl -chr $chr -outfile $outdir/$chr.fa -mask -dust -softmask -dbname $dbname -dbhost $dbhost";
    $fh->open("| bsub -q acari -m ecs2_hosts -o $outdir/jobs/lsf-out-$chr");
    $fh->print($command);
    $fh->close;
}


sub get_chr_names{
    my $db = shift;
    my $q = qq(   SELECT   chr.name
		  FROM     chromosome chr
		  );
    
    my $sth = $db->prepare($q) || $db->throw("can't prepare: $q"); 
    my $res = $sth->execute    || $db->throw("can't execute: $q");
    
    my @chrs;
    while( my ($chr_name) = $sth->fetchrow_array ) {
	push (@chrs,$chr_name);
    }
    return @chrs;
}
