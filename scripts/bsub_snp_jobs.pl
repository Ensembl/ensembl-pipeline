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

use Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor; 
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

#my $dbtype     = 'rdb';
#my $host_snp   = 'ensrv3';
my $host_snp   = 'ecs1a';
my $port       = '';
#my $dbname_snp = 'homo_sapiens_snp_110';
my $dbname_snp = 'snp120';
my $dbuser     = 'ensadmin';
my $dbpass     = '***';
my $module_snp = 'Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor';

#&GetOptions ( 
#               'dbtype:s' => \$dbtype,
#               'host:s'   => \$host,
#               'port:n'   => \$port,
#               'dbname:s' => \$dbname, 
#               'dbuser:s' => \$dbuser,
#               'dbpass:s' => \$dbpass,
#               'module:s' => \$module,
#               );

my $locator_snp = "$module_snp/host=$host_snp;port=$port;dbname=$dbname_snp;user=ensro";

my $debug=0;

my $snpdb =  Bio::EnsEMBL::DBLoader->new($locator_snp); 

my $start_refnum = 1;
my $length_refnum = 10000;
#my $hitcount = $snpdb->get_Hitcount;
my $max_refsnpid = $snpdb->get_max_refsnpid;
#my $hitcount = 50000;
print "the total count in Hit is $hitcount\n";

my $count;
for ($start_refnum =1; $start_refnum <=$max_refsnpid; $start_refnum+=$length_refnum) {
  my $end_refnum = $start_refnum+$length_refnum-1;
  print "job_num is ", ++$count, "start_refnum is $start_refnum and end_refnum is $end_refnum\n";
  `bsub -q ecslong -o /work4/yuan/OUT1/out_$start_refnum -e /work4/yuan/ERROR1/error_$start_refnum /nfs/acari/yuan/ensembl/test/snp_cross_match_new.pl $start_refnum $end_refnum`;
}























































































