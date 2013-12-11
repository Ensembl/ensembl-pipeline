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
use Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateGenomic;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$| = 1;

my $dnadb=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'homo_sapiens_core_110',-host=>'ensrv3',-user=>'ensro');

my $mousedb=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'mouse',-host=>'ensrv3',-user=>'ensadmin');

my $newdb=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'exonerate_genomic',-host=>'ecs1d',-user=>'ensadmin', -dnadb=>$dnadb);

my $estgen = Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateGenomic->new(-querydb => $dnadb, -targetdb => $mousedb, -dbobj => $newdb);
$estgen->fetch_input();
$estgen->run;



