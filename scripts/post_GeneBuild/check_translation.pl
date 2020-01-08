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


=head1 NAME

  dump_translations.pl

=head1 SYNOPSIS
 
  dump_translations.pl

=head1 DESCRIPTION

dump_translations.pl dumps out the translations of all the genes in a database specified in GeneConf.pm
It\'s a stripped down version of gene2flat.

=head1 OPTIONS

=cut

use warnings ;
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case);

my $dbhost;
my $dbuser = 'ensro';
my $dbname;
my $dbpass = undef;

my $dnadbhost;
my $dnadbuser = 'ensro';
my $dnadbname;
my $dnadbpass = undef;

my $genetype;


my $tstable_id;
my $t_id;
my $gstable_id;
my $g_id;

GetOptions(
	    'tstable_id:s' => \$tstable_id,
	    't_id:s' => \$t_id,
	    'g_id:s' => \$g_id,
	    'gstable_id:s' => \$gstable_id,
	    'host|dbhost|h:s'        => \$dbhost,
	    'dbname|db|D:s'        => \$dbname,
	    'dnadbhost:s'     => \$dnadbhost,
	    'dnadbname:s'     => \$dnadbname,
	    'genetype:s'      => \$genetype,
	    );

unless ( $dbhost && $dbname && $dnadbhost && $dnadbname && ( $t_id || $tstable_id || $g_id || $gstable_id) ){
  print STDERR "script to check the translation from the transcripts or genes in a database\n";
  
  print STDERR "Usage: $0 [-dbname -dbhost -dnadbname -dnadbhost -genetype] -t_id -tstable_id -g_id -gstable_id\n";
  exit(0);
}




my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dnadbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dnadbpass,
					      );


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					   );



print STDERR "connected to $dbname : $dbhost\n";

my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*STDOUT ) ;

if ( $t_id){

  # first method
  my $tadaptor = $db->get_TranscriptAdaptor;
  my $trans    = $tadaptor->fetch_by_dbID($t_id);
  my $tseq     = $trans->translate();
  $tseq->desc("Transcript dbID: $t_id, transcript from TranscriptAdaptor");
  $seqio->write_seq($tseq);
  
  
}



