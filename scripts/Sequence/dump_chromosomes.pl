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


=head1 NAME

  dump_chromosomes.pl

=head1 SYNOPSIS
 
=head1 DESCRIPTION

=head1 OPTIONS

  -repeatmask (off by default)
  -dbname
  -dbport
  -dbhost
  -dbuser
  -dbpass
  -outdir (the directory where the chromosome sequences will 
	   be written, one file per chromosome).
  -chr (specify a chromosome name to dump just one chromosome)

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

my $dbname;
my $dbport = 3306;
my $dbhost;
my $dbuser;
my $dbpass;
my $repeatmask;
my $dust;
my $outdir;
my $chr;

GetOptions('dbname|db|D:s'   => \$dbname,
	    'dbport|port|P:s'   => \$dbport,
	    'dbhost|host|h:s'   => \$dbhost,
	    'dbuser|user|u:s'   => \$dbuser,
	    'dbpass|pass|p:s'   => \$dbpass,
	    'repeatmask' => \$repeatmask,
	    'dust'       => \$dust,
	    'outdir:s'   => \$outdir,
	    'chr:s'      => \$chr);

# Whack out some help if lost...
if(!defined $dbname || !defined $dbhost || !defined $dbuser || !defined $outdir){

  print "Usage: dump_chromosomes.pl -dbname [name] -dbhost [host] -dbuser [user]\n" .
    "-dbport [port] -dbpass [pass] -outdir [directory for chr files] -repeatmask\n" .
      "-dust -chr [chr name]\n";

  exit(0);
}


# Connect to db and stuff

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $dbname,
					     -user   => $dbuser,
					     -host   => $dbhost,
					     -port   => $dbport,
					     -pass   => $dbpass);

my $slice_adaptor = $db->get_SliceAdaptor;
my $chromosome_adaptor = $db->get_ChromosomeAdaptor;

# Figure out what kind of masking is wanted

my @mask_logicnames = ();
push (@mask_logicnames, 'RepeatMask') if $repeatmask;
push (@mask_logicnames, 'Dust')       if $dust;

# Loop through each chromosome in our database and dump each sequence 
# into its own eponymous file.

foreach my $chromosome (@{$chromosome_adaptor->fetch_all}){

  print "Dumping chromosome " . $chromosome->chr_name;

  my $outfile = $outdir . '/' . $chromosome->chr_name . '.fa';

  if (-e $outfile){
    print " ... output file exists already - skipping\n";
    next
  }

  if ($chr && $chr ne $chromosome->chr_name){
    print " ... skipping\n";
    next
  }

  print "\n";

  # Open a Bio::SeqIO output stream

  my $seqio = Bio::SeqIO->new(-file   => ">$outfile",
			      -format => 'Fasta');


  my $slice = $slice_adaptor->fetch_by_chr_name($chromosome->chr_name);

  unless ($repeatmask){
    $seqio->write_seq($slice);
    next
  }

  my $rm_seq = $slice->get_repeatmasked_seq(\@mask_logicnames,1);
  $rm_seq->display_id($slice->display_id . " @mask_logicnames");


  $seqio->write_seq($rm_seq);

}
