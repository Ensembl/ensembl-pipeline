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


use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Connection to the target DB
my $corehost   = 'ensembldb.ensembl.org';
my $coreport   = 5306;
my $coreuser   = 'anonymous';
my $rnaseqhost = 'ensembldb.ensembl.org';
my $rnaseqport = '5306';
my $rnasequser = 'anonymous';
my $coredbname;
my $rnaseqdbname;
my $id_file;
my $outfile;
my $exon;
my $verbose;
my $help;

&GetOptions (
        'corehost=s'     => \$corehost,
        'coreport=s'     => \$coreport,
        'coreuser=s'     => \$coreuser,
        'coredbname=s'   => \$coredbname,
        'rnaseqhost=s'   => \$rnaseqhost,
        'rnaseqport=s'   => \$rnaseqport,
        'rnasequser=s'   => \$rnasequser,
        'rnaseqdbname=s' => \$rnaseqdbname,
        'id_file=s'      => \$id_file,
        'outfile=s'      => \$outfile,
        'exon_info!'     => \$exon,
        'verbose!'       => \$verbose,
        'help!'          => \$help,
        );

&usage if (!$coredbname or !$rnaseqdbname or !$id_file or $help);

my @ensembl_list_ids;
my @xref_list_ids;

$outfile = $id_file.'.out' if (!$outfile);
# We want a file with ad id on each line
open(RF, $id_file) || die ("Could not open $id_file\n");
while (my $line = <RF>) {
    if ($line =~ /\s*(ENS\w+)/) {
        push(@ensembl_list_ids, $1);
    }
    elsif ($line =~ /\s*([a-zA-Z0-9._]+)/)  {
        push(@xref_list_ids, $1);
    }
    else {
        print "This id doesn't seem to be good! $line";
    }
}
close(RF);

# Connection to the core and the rnaseq databases
my $coredb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host => $corehost,
        -port => $coreport,
        -user => $coreuser,
        -dbname => $coredbname);
my $rnaseqdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host => $rnaseqhost,
        -port => $rnaseqport,
        -user => $rnasequser,
        -dbname => $rnaseqdbname);

# We need to fetch all the genes in the core database
my $corega = $coredb->get_GeneAdaptor();
my $rnaseqsa = $rnaseqdb->get_SliceAdaptor();
my $rnaseqdafa = $rnaseqdb->get_DnaAlignFeatureAdaptor();
my @genes;
foreach my $id (@ensembl_list_ids) {
    push(@genes, $corega->fetch_by_stable_id($id));
}
foreach my $id (@xref_list_ids) {
    push(@genes, @{$corega->fetch_all_by_external_name($id)});
}
print 'Number of genes to look for: ', scalar(@genes), "\n" if $verbose;

open(WF, '>'.$outfile) || die("Could not open $outfile to write the results!\n");
print WF "EnsEMBL_id\tTissue\tChromosome\tStart\tEnd\tStrand\tLength\tTotal_number_of_reads";
print WF "\tNumber_of_reads_intron1\tNumber_of_reads_intron2..." if ($exon);
print WF "\n";
# For each gene we will fetch the region (slice) where the gene is, but from
# the rnaseq database. After we will look at all the models that overlap the
# region and process them
foreach my $gene (@genes) {
    my $gene_slice = $gene->feature_Slice;
    my $rnaseqslice = $rnaseqsa->fetch_by_name($gene_slice->name);
    print WF $gene->display_id, "\t", $gene->analysis->logic_name, "\t", $gene->seq_region_name, "\t", $gene->start, "\t", $gene->end, "\t", $gene->strand, "\t", $gene->length, "\n";
    # We load the exons to speed up a bit
    my $transcripts = $rnaseqslice->get_all_Transcripts(1);
    foreach my $transcript (@$transcripts) {
        next unless ($gene->strand == $transcript->strand);
        print WF $transcript->display_id, "\t", $transcript->analysis->logic_name, "\t", $transcript->seq_region_name, "\t", $transcript->seq_region_start, "\t", $transcript->seq_region_end, "\t", $transcript->strand, "\t", $transcript->length, "\t";
        my $tslice = $transcript->feature_Slice();
        # We fetch the dna features that overlap the gene
        my $features = $rnaseqdafa->fetch_all_by_Slice($tslice);
        my @introns;
        my $start = 0;
        my $end = 0;
        # We want to know the boundaries of the introns
        # The exon are sorted, 5' end first
        foreach my $exon (@{$transcript->get_all_Exons()}) {
            if ($gene->strand == -1) {
                if($end) {
                    push(@introns, [$exon->display_id, $exon->seq_region_end+1, $end-1]);
                }
                $end = $exon->seq_region_start;
            }
            else {
                if($end) {
                    push(@introns, [$exon->display_id, $end+1, $exon->seq_region_start-1]);
                }
                $end = $exon->seq_region_end;
            }
        }
        my $index = 1;
        my $treads_count = 0;
        my $ereads_count;
        foreach my $intron (@introns) {
            # We loop inside the feature array to get the introns that matches this intron boundaries
            my $read_count = 0;
            foreach my $feature (@{$features}) {
                # We have to call seq_region_start/seq_region_end because we are on a sub slice
                # start/end will give the position starting from the begining of the gene
                next unless ( $intron->[1] == $feature->seq_region_start and $intron->[2] == $feature->seq_region_end);
                # We print all the features that matches the boundaries
                # We count all the reads that span over a splice site
                $read_count += $feature->score;
            }
            $ereads_count .= "\t".$read_count;
            $treads_count += $read_count;
#                    print WF "\tIntron $index  ", $intron->[1], ':', $intron->[2], "\tNumber of spanning reads: ", $read_count, "\n  Exon ", ++$index, '  ';
        }
        $treads_count .= $ereads_count if ($exon);
        print WF $treads_count, "\n";
    }
}
close(WF);

sub usage {
    print <<EOF
perl $0 -coredbname <dbname> [-corehost <dbhost>] [-coreport <dbport>] [-coreuser <dbuser>] -rnaseqdbname <dbname> [-rnaseqhost <dbhost>] [-rnaseqport <dbport>] [-rnasequser <dbuser>] -id_file <path_to_your_file> [-outfile <path_to_your_output_file>] [-exon_info] [-verbose]

   -coredbname   Name of the core database
   -corehost     Name of the host for the core DB, default is ensembldb.ensembl.org
   -coreport     Port for the connection, default is 5603
   -coreuser     User name, default is anonymous
   -rnaseqdbname Name of the rnaseq database, contains the RNASeq models
   -rnaseqhost   Name of the host for the rnaseq DB, default is ensembldb.ensembl.org
   -rnaseqport   Port for the connection, default is 5603
   -rnasequser   User name, default is anonymous
   -id_file      File that contains the list of HGNC gene names or Ensembl gene ids (BRCA2, ENSG00000083642,...)
   -outfile      Output file, if not set, it will concatenate '.out' to the name of id_file
   -exon_info    Display the read spanning intron count for each intron, otherwise it's the sum for the entire transcript
   -verbose
EOF
;
}
