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


=head1 NAME

run_blat

=head1 SYNOPSIS
 

=head1 DESCRIPTION


=head1 OPTIONS

    -input_id  The input id for the RunnableDB
    -runnable  The name of the runnable module we want to run
    -analysis  The number of the analysisprocess we want to run

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::SeqIO;

use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::Exonerate qw (
							      EST_REFDBNAME
							      EST_REFDBUSER
							      EST_REFDBHOST
							      EST_GENOMIC
							      EST_EXONERATE_OPTIONS
							     );


use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $runnable;
my $input_id;
my $write  = 0;
my $check  = 0;
my $params;
my $pepfile;
my $analysis;
my $query_seq;

# can override db options on command line
GetOptions( 
	    'runnable:s'    => \$runnable,
	    'analysis:s'    => \$analysis,
	    'write'         => \$write,
	    'check'         => \$check,
	    'query_seq:s'   => \$query_seq,
	   );

$| = 1;

die "No runnable entered" unless defined ($runnable);
(my $file = $runnable) =~ s/::/\//g;
require "$file.pm";

if ($check) {
   exit(0);
}

print STDERR "REFDB: $EST_REFDBHOST : $EST_REFDBUSER : $EST_REFDBNAME\n";

############################################################
# this is just a refdb where we get dna/assembly/chromosomes from
my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       -host             => $EST_REFDBHOST,
					       -user             => $EST_REFDBUSER,
					       -dbname           => $EST_REFDBNAME,
					       );

############################################################
# I should have the correct analysis in the refdb analysis table
my $analysis_obj = $refdb->get_AnalysisAdaptor->fetch_by_logic_name($analysis);


die "Need ana analysis object" unless( defined $analysis_obj );
die "No input entered" unless( defined ($query_seq));

############################################################
# convert the multiple fasta file with query est/rna sequences into an array of seq features
my @sequences;
my $seqio= Bio::SeqIO->new(
			   -format => "Fasta",
			   -file   => "$query_seq",
			  );

############################################################
# read the query and create Bio::Seqs
while( my $seq = $seqio->next_seq() ){
  if (defined($seq) && !($seq->seq eq '') && !($seq->display_id eq '') ){
    push( @sequences, $seq );
  }
  else{
    print STDERR "problems getting sequence ".$seq->display_id."\n".$seq->seq."\n";
  }
}
#print STDERR "got ".scalar(@sequences)." sequence objects\n";


############################################################
# create runnableDB
my $runobj = "$runnable"->new(-db         => $refdb,
			      -input_id   => \@sequences,
			      -rna_seqs   => \@sequences,
			      -analysis   => $analysis_obj,
			      -database   => $EST_GENOMIC,
			      -query_type => 'dna',
			      -target_type=> 'dna',
			      -options    => $EST_EXONERATE_OPTIONS,
			     );

$runobj->fetch_input;
$runobj->run;
 

my @out = $runobj->output;

if ($write) {
  $runobj->write_output;
}

############################################################
# do this if the file is local
#if ( $query_seq =~ /\/tmp\/\S+/ ){
#    unlink( $query_seq );
#}












