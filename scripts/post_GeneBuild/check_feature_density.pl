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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);

my $dbhost;
my $dbuser = 'ensro';
my $dbname;
my $dbpass = undef;

my $feature_name;

GetOptions(
	    'feature_name:s' => \$feature_name,
	    'host|dbhost|h:s'        => \$dbhost,
	    'dbname|db|D:s'        => \$dbname,
	   );

unless ( $feature_name ){
  print STDERR "script to check the density of features\n";
  
  print STDERR "Usage: $0 [-dbname -dbhost]\n";
  exit(0);
}


#my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
#					       '-host'   => $dnadbhost,
#					       '-user'   => $dnadbuser,
#					       '-dbname' => $dnadbname,
#					       '-pass'   => $dnadbpass,
#					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
#					    '-dnadb'  => $dnadb,
					   );


print STDERR "connected to $dbname : $dbhost\n";




my @chrs = &get_chrnames($db);

foreach my $chr (@chrs){
    print "Chr: $chr\n";
    my ($contigs,$position,$id) = &get_contigs($db,$chr);
    
    foreach my $contig (@$contigs){
	
      my ($count_dna_features,$count_protein_features) = &count_features($db,$contig);
      print $position->{$contig}."\t".$id->{$contig}."\t".$count_dna_features."\t".$count_protein_features."\n";
    }
}

sub get_chrnames{
    my $db   = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my @chrnames;

  my $q = qq( SELECT distinct(chromosome_id) 
	      FROM   assembly
	      );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr) = $sth->fetchrow_array) {
      push (@chrnames, $chr);
  }
  return @chrnames;
}


sub get_contigs{
  my $db   = shift;
  my $chr = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my %position;
  my %id;
  my @contigs;

  my $q = qq( SELECT s.chr_start,s.contig_id, c.contig_id 
	      FROM assembly s, contig c
	      WHERE s.contig_id = c.contig_id 
	      AND   s.chromosome_id = $chr
	      ORDER BY s.chr_start
	      );
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr_start, $raw_id, $id) = $sth->fetchrow_array) {
      push (@contigs, $raw_id);
      $position{$raw_id} = $chr_start;
      $id{$raw_id} = $id;
  }
  return (\@contigs,\%position,\%id);
}

sub count_features{
    my $db   = shift;
    my $contig = shift;
    
    if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
	die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
    }
    
    my %position;
    my %id;
    my @contigs;
    
    
    my $q = qq( SELECT count(dna_align_feature_id)
		FROM   dna_align_feature
		WHERE  contig=$contig
		);
    
    my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
    my $res = $sth->execute    || $db->throw("can't execute: $q");
    
    my $dna_align_feature_count = $sth->fetchrow_array;

    my $q = qq( SELECT count(protein_align_feature_id)
		FROM   protein_align_feature
		WHERE  contig=$contig
		);
    
    my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
    my $res = $sth->execute    || $db->throw("can't execute: $q");
    
    my $protein_align_feature_count = $sth->fetchrow_array;
    
    return ( $dna_align_feature_count,  $protein_align_feature_count );
}

