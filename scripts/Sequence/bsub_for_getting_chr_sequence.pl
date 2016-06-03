#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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
my $coord_system;
my $coord_system_version;
my $dbpass;
my $dbport = 3306;
GetOptions(
            'dbname|db|D:s'       => \$dbname,
            'dbhost|host|h:s'       => \$dbhost,
            'outdir:s'       => \$outdir,
            'dbuser|user|u:s' => \$dbuser,
            'dbport|port|P:s' => \$dbport,
            'dbpass|pass|p:s' => \$dbpass,
            'coord_system|cs_name:s' => \$coord_system,
            'coord_system_version|cs_version:s' => \$coord_system_version,
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
