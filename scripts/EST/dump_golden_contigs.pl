#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

  dump_golden_contigs.pl

=head1 SYNOPSIS
 
 dump_golden_contigs.pl

=head1 DESCRIPTION

  dumps out masked golden contig sequences

=head1 OPTIONS

  -dbname
  -dbhost
  -path
  -out

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_REFDBNAME
					EST_REFDBHOST
					EST_REFDBUSER
				       );

$| = 1;

my $dbname     = '';
my $dbuser     = ''; # always use read only
my $dbhost     = '';
my $path       = '';
my $outfile    = 'out.fa';

# try ESTConf for configuration options
$dbname  = $EST_REFDBNAME;
$dbuser  = $EST_REFDBUSER;
$dbhost  = $EST_REFDBHOST;

# otherwise get them from the command line
GetOptions( 
	    'dbname|db|D:s'     => \$dbname,
	    'host|dbhost|h:s'     => \$dbhost,
	    'outfile:s'    => \$outfile,
	    'path|cs_version:s'       => \$path,
	   );

# usage
if(!defined $dbname    ||
   !defined $dbhost    ||
   !defined $path    
  ){
  print  "USAGE: dump_golden_contigs.pl -dbname dbname -host host -path path\n optional:\n -outfile to specify a filename other than out.fa\n";
  exit(1);
}

open (OUT, ">$outfile") or die "Can't open outfile [$outfile]\n";

# global stuff
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host   => $dbhost,
					    -dbname => $dbname,
					    -user   => $dbuser,
					   );

my $query = "select c.name from assembly a, contig c, meta m where a.contig_id=c.contig_id and a.type = m.meta_value and m.meta_key='assembly.default'";

my $sth = $db->prepare($query);
my $res = $sth->execute;

my $counter = 0;
while (my $id= $sth->fetchrow_array){
  $counter++;

  my $contig = $db->get_RawContigAdaptor->fetch_by_name($id);
  print OUT ">" . $contig->id . "\n" . $contig->get_repeatmasked_seq->seq . "\n";
  
  print STDERR "processed $counter contigs\n";
}

close OUT or die "Can't close outfile [$outfile]\n";
