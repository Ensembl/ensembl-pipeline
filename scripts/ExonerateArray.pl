#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;

my $host   = 'ecs4';
my $port       = '3351';
my $user   = 'ensro';
my $dbname     = 'homo_sapiens_core_20_34b';

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host =>$host,
                                             -user  =>$user,
                                             -port => $port,
                                             -dbname=>$dbname);

#my $query_file = "/acari/scratch4/ensembl/yuan/array/HG-U95Av2_probe_fasta";
my $query_file = "/acari/scratch4/ensembl/yuan/array/test.fa";
my $target_file = "/data/blastdb/Ensembl/Human/NCBI34/genome/masked_dusted/*";
my $exonerate = "/usr/local/ensembl/bin/exonerate-0.7.1";


my $exonarray = Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate2Array->new(
									 -db          => $db,
									 -query_file  => $query_file ,
									 -target_file => $target_file,
									 -exonerate   => $exonerate,
									);

$exonarray->fetch_input;
$exonarray->run;
$exonarray->write_output;
