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
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Pipeline::Config::FirstEF qw (
						FEF_WRITEDBNAME
						FEF_WRITEDBHOST
						FEF_WRITEDBUSER
						FEF_WRITEDBPASS
						FEF_REFDBNAME 
						FEF_REFDBHOST
						FEF_REFDBUSER
					       );

my $input_id;
my $write = 0;

GetOptions( 
	     'input_id:s'  => \$input_id,
             'write'       => \$write
	     );




my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $FEF_WRITEDBNAME,
					     -host   => $FEF_WRITEDBHOST,
					     -user   => $FEF_WRITEDBUSER,
					     -pass   => $FEF_WRITEDBPASS);

my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $FEF_REFDBNAME, 
						-host   => $FEF_REFDBHOST,
						-user   => $FEF_REFDBUSER);

$db->dnadb($dnadb);

my $analysis_adaptor = $db->get_AnalysisAdaptor;
my $analysis = $analysis_adaptor->fetch_by_logic_name('firstef');

my $firstef = Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF->new(-input_id  => $input_id,
							       -db        => $db,
							       -analysis  => $analysis );

$firstef->fetch_input;
$firstef->run;

my @output = $firstef->output;

print "Found " . scalar @output . " features\n";

if ($write) {
  $firstef->write_output;
}
