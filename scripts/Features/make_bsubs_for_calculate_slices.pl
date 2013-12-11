#!/usr/bin/env perl


# Copyright [1999-2013] Genome Research Ltd. and the EMBL-European Bioinformatics Institute
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


my $dbhost;
my $dbname;
my $outdir;

GetOptions(
	    'dbhost|host|h:s'        => \$dbhost,
	    'dbname|db|D:s'        => \$dbname,
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
