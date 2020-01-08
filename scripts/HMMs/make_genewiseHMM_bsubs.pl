#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
my $species = 'mouse';
#my $species = 'rat';
my $dbhost = 'ecs2f';
my $dbuser = 'ensro';
my $dbname = "genewisedb_".$species;
my $outdir;
my $threshold = 50;


GetOptions(
	    'dbname|db|D:s'       => \$dbname,
	    'dbhost|host|h:s'       => \$dbhost,
	    'outdir:s'       => \$outdir,
	    );

my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
					   -host  => $dbhost,
					   -user  => $dbuser,
					   -dbname=> $dbname
					  );


my $input_ids = &get_input_ids($db,$threshold);
my %ids = %$input_ids;

foreach my $hmm ( keys %ids ){
  foreach my $dna_id ( keys %{$ids{$hmm}} ){
    
    my $command = "bsub -q acari -C0 -f  \"/ecs2/work1/eae/GeneWiseHMM/HMMs/$hmm.hmm > /tmp/$hmm.hmm\" -o /ecs2/work1/eae/GeneWiseHMM/genewiseHMM_".$species."_jobs_dir/$hmm-$dna_id.out -e /ecs2/work1/eae/GeneWiseHMM/genewiseHMM_".$species."_jobs_dir/$hmm-$dna_id.err -E \"/nfs/acari/eae/ensembl/ensembl-pipeline/scripts/HMMs/run_genewiseHMM.pl -check \" /nfs/acari/eae/ensembl/ensembl-pipeline/scripts/HMMs/run_genewiseHMM.pl -input_id $dna_id  -hmm /tmp/$hmm.hmm";
    print $command."\n";
  }    
}

sub get_input_ids{
  my $db = shift;
  my $threshold = shift;
  my $q = "SELECT protein_id, dna_id 
           FROM   genewisedb
           WHERE  bit_score > $threshold
             AND  bit_score < 10000";
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  my %ids;
  while( my ($protein_id,$dna_id) = $sth->fetchrow_array) {
    $ids{$protein_id}{$dna_id} = 1;
  }
  return \%ids;
}
