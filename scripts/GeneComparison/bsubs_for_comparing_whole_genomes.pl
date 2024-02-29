#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use warnings ;
use strict;
use IO::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $fh = new IO::File;

# connect to the database
my $dbhost;
my $dbuser = 'ensro';
my $dbname;
my $outdir;

GetOptions(
	    'dbname|db|D:s'       => \$dbname,
	    'host|dbhost|h:s'       => \$dbhost,
	    'outdir:s'       => \$outdir,
	    );

unless( $dbhost &&  $dbname && $outdir ){
  print STDERR "Usage: $0 -dbname dbname -dbhost dbhost -outdir outdir\n";
  print STDERR "Need database information and the directory (full path) where to send the chromosome files\n";
  print STDERR "It asumes 'compare_Exons.pl' is in your path\n";
  exit(0);
}

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
					   -host  => $dbhost,
					   -user  => $dbuser,
					   -dbname=> $dbname
					   );



my @chr_names = &get_chr_names($db);

foreach my $chr (@chr_names){
    my $command = "compare_Exons.pl -chr $chr";
    $fh->open("| bsub -q acari -m ecs2_hosts -o $outdir/lsf-out-$chr -e $outdir/$chr-comparison");
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
