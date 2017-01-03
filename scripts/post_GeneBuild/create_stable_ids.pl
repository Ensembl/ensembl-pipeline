#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw); 
use Bio::SeqIO;
use Getopt::Long;
#use Bio::EnsEMBL::Pipeline::ESTConf;

my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $label; 
my $sql; 
my $dbport = '3306' ;
GetOptions(
	    'dbname|db|D:s'    => \$dbname,
	    'host|dbhost|h:s'    => \$dbhost,
            'port|dbport|P:i'    => \$dbport,
	    'label:s'     => \$label,
            'sql!'        => \$sql, 

	   );

unless ( $dbname && $dbhost && $label ){
  print STDERR "script to fake some stable ids for genes, transcripts, exons and translations\n";
  print STDERR "Usage: $0 -dbname -dbhost -label ( ENSEST, ENS, ENSRNO, ...)\n";
  exit(0);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
                                            '-port'   => $dbport,
					   );

print STDERR "connected to databse  $dbname : $dbhost\n";

my $gene_file        = "gene_stable_id";
my $transcript_file  = "transcript_stable_id";
my $translation_file = "translation_stable_id";
my $exon_file        = "exon_stable_id";

print "Results will be written to $gene_file $transcript_file etc\n";  

open (GENE,">$gene_file")               or die("unable to open file $gene_file");
open (TRANSCRIPT,">$transcript_file")   or die("unable to open file $transcript_file");
open (TRANSLATION,">$translation_file") or die("unable to open file $translation_file");
open (EXON,">$exon_file")               or die("unable to open file $exon_file");

my $count_gene = 0;
my $count_transcript = 0;
my $count_exon = 0;

my %seen_gene;
my %seen_exon;
my %seen_transcript;

my $gene_adaptor = $db->get_GeneAdaptor;
foreach my $gene_id( @{$gene_adaptor->list_dbIDs} ) { 

    my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id); 
  eval {

    ####################
    # create gene id 
    ####################
    $count_gene++;
    my $gene_id = $gene->dbID;
    my $gene_stable_id = &make_stable_id($label."G",$gene);
    #  927491 | ENSestG00000000001 |       1 | 0000-00-00 00:00:00 | 0000-00-00 00:00:00 |  
    if ( $sql ) { 
      make_sql($gene_id, $gene_stable_id, "gene") ;  
    } else { 
      print GENE "$gene_id\t$gene_stable_id\t1\t0\t0\n";
    }
    
    #######################
    # create transcript id
    #######################
    foreach my $trans ( @{$gene->get_all_Transcripts} ) {
      $count_transcript++;
      my $transcript_id = $trans->dbID;
      my $transcript_stable_id = &make_stable_id($label."T",$trans);
      #  927491 | ENSestT00000000001 |       1 
     if ( $sql ) { 
      make_sql($transcript_id, $transcript_stable_id, "transcript") ;  
     } else  { 
      print TRANSCRIPT "$transcript_id\t$transcript_stable_id\t1\n";
     } 
      #######################
      # create translation id
      ####################### 
      if ( defined  $trans->translation ) { 
        my $translation           = $trans->translation;
        my $translation_id        = $translation->dbID;
        my $translation_stable_id = &make_stable_id($label."P",$translation);
        #  927491 | ENSestP00000000001 |       1 
        if ( $sql ) { 
         make_sql($translation_id, $translation_stable_id, "translation") ;   
        } else { 
          print TRANSLATION "$translation_id\t$translation_stable_id\t1\n"; 
        }
      }
    }

    ####################
    # create exon ids 
    ####################
    my @exons = @{$gene->get_all_Exons}; 
    if  ( scalar(@exons) == 0 ) {  
        throw( "gene " . $gene->dbID . " " . $gene->biotype . " has 0 exons \n") ;  
    } 
    foreach my $exon ( @exons ){
      $count_exon++;
      my $exon_id = $exon->dbID;
      my $exon_stable_id = &make_stable_id($label."E",$exon);
      #  927491 | ENSestE00000000001 |       1 
        if ( $sql ) { 
          make_sql($exon_id, $exon_stable_id, "exon") ;   
        } else {  
          print EXON "$exon_id\t$exon_stable_id\t1\t0\t0\n";
        }
    }
  }; 
   warn $@ if $@;
}

close EXON;
close TRANSCRIPT;
close TRANSLATION;
close GENE;


print STDERR "Files gene_stable_id, transcript_stable_id, translation_stable_id , exon_stable_id written\n";  


sub make_sql{ 
  my ($id, $stable_id,$table) = @_ ;  

  print "insert into $table\_stable_id values ( $id, \"$stable_id\",1,now(),now());\n";;
}
 
sub make_stable_id{
  my ($label,$object) = @_;
  # $label is ENSestG, etc...

  my $number = $object->dbID;
  my $string = "$number";
  while ( length($string) < 11 ){
    $string = "0".$string;
  }
  $string = $label.$string;
  return $string;
}
