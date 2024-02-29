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


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)

    -dbport      For RDBs, what port to connect to (port= in locator)

    -dbname      For RDBs, what name to connect to (dbname= in locator)

    -dbuser      For RDBs, what username to connect as (dbuser= in locator)

    -dbpass      For RDBs, what password to use (dbpass= in locator)

    -dnadbname, dnadbhost, dnadbuser, dnadbpass, dnadbport  - Options for using a reference DNA database

    -file        File containing list of gene ids to be flagged as pseudogenes
    
    -analysis    Logic name for pseudogene analysis

=cut

use warnings ;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);


my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;
my $dnadbhost;
my $dnadbuser;
my $dnadbpass;
my $dnadbport;
my $dnadbname;
my $file;
my $logic_name;
my $analysis;

GetOptions( 
	    'dbhost|host|h=s'      => \$dbhost,
	    'dbname|db|D=s'      => \$dbname,
	    'dbuser|user|u=s'      => \$dbuser,
	    'dbpass|pass|p=s'      => \$dbpass,
	    'dbport|port|P=s'      => \$dbport,
	    'dnadbhost=s'   => \$dnadbhost,
	    'dnadbname=s'   => \$dnadbname,
	    'dnadbuser=s'   => \$dnadbuser,
	    'dnadbpass=s'   => \$dnadbpass,
	    'dnadbport=s'   => \$dnadbport,
	    'file=s'        => \$file,
	    'analysis=s'    => \$logic_name,
	   ) or die("couldn't get options");


&check_options();


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor (
					     -host    => $dbhost,
					     -user    => $dbuser,
					     -dbname  => $dbname,
					     -pass    => $dbpass,
					     -port    => $dbport,
					    );


if($dnadbhost){
  my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor (
						  -host    => $dnadbhost,
						  -user    => $dnadbuser,
						  -dbname  => $dnadbname,
						  -pass    => $dnadbpass,
						  -port    => $dnadbport,
						 );
  $db->dnadb($dnadb);
}

$analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
if (!defined $analysis){
  die("Can't fetch analysis for logic name [$logic_name] from [$dbname] on [$dbhost]\n");
}
  
my $ga = $db->get_GeneAdaptor;
  
open(FH, $file) or die("couldn't open ".$file);

my %genes_updated;
while(<FH>){
  chomp;
  my @values = split;
  my $gene_id = $values[0];
  if(!$genes_updated{$gene_id}){
    my $gene = $ga->fetch_by_dbID($gene_id);
    $gene->type($logic_name);
    $gene->analysis($analysis);

    # need to delete translation and set transcript->translation id to be 0
    print STDERR "updating $gene_id\n";
    $ga->update($gene);
    $genes_updated{$gene_id} = 1;
    foreach my $transcript(@{$gene->get_all_Transcripts}){
      $db->get_TranslationAdaptor->remove($transcript->translation);
      my $sth = $db->prepare("UPDATE transcript SET translation_id=0 WHERE transcript_id = ?");
      $sth->execute($transcript->dbID);
    }
  }
}

close FH or die ("Couldn't close [$file]\n");

### SUBROUTINES ###

sub check_options{
  &usage unless defined ($dbname && $dbhost && $dbuser);
  if(defined $dnadbname){
    &usage unless defined ($dnadbhost && $dnadbuser);
  }
  &usage unless defined $logic_name;
  &usage unless defined $file;

}

sub usage{
  print "Usage: update_pseudogenes.pl -dbname -dbuser -dbhost -dbpass [-dnadbname -dnadbuser -dnadbhost -dnadbpass] -file -analysis\n";
  exit(0);
}
