#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# Connection to the rnaseq database
my $rnaseqdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host =>  'ensembldb.ensembl.org',
        -port =>  5306,
        -user =>  'anonymous',
        -dbname =>  'homo_sapiens_rnaseq_68_37');

# We need the GeneAdaptor and the DnaFeatureAdaptor
my $rnaseqsa = $rnaseqdb->get_SliceAdaptor();
my $rnaseqaa = $rnaseqdb->get_AnalysisAdaptor();
my $rnaseqdafa = $rnaseqdb->get_DnaAlignFeatureAdaptor();
foreach my $slice (@{$rnaseqsa->fetch_all('toplevel')}) {
    foreach my $analysis (@{$rnaseqaa->fetch_all()}) {
        # The analysis has a name finishing by rnaseq
        my $logic_name = $analysis->logic_name;
        next if ($logic_name !~ /rnaseq/);
        print 'Looking at genes from analysis: ', $logic_name, "\n";
        # We load the gene completely as we will use exons and transcripts
        foreach my $gene (@{$slice->get_all_Genes_by_type(undef, $logic_name, 1)}) {
            print $gene->display_id, ':', $gene->start, ':', $gene->end, ':', $gene->strand, ' on chromosome ', $gene->seq_region_name, "\n";
            # Rnaseq models have only one transcript per gene
            foreach my $transcript (@{$gene->get_all_Transcripts()}) {
                print $transcript->display_id, ':', $transcript->start, ':', $transcript->end, '  length: ', $transcript->length, "\n";
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
                            push(@introns, [$exon->display_id, $exon->end+1, $end-1]);
                        }
                        $end = $exon->start;
                    }
                    else {
                        if($end) {
                            push(@introns, [$exon->display_id, $end+1, $exon->start-1]);
                        }
                        $end = $exon->end;
                    }
                }
                my $index = 1;
                print '  Exon ', $index, '  ';
                foreach my $intron (@introns) {
                    print $intron->[0], '  ', $intron->[1], ':', $intron->[2], "\n";
                    # We loop inside the feature array to get the introns that matches this intron boundaries
                    my $read_count = 0;
                    foreach my $feature (@{$features}) {
                        # We have to call seq_region_start/seq_region_end because we are on a sub slice
                        # start/end will give the position starting from the begining of the gene
                        next unless ( $intron->[1] == $feature->seq_region_start and $intron->[2] == $feature->seq_region_end);
                        # We print all the features that matches the boundaries
#                       print "\t", $feature->display_id, ':', $feature->seq_region_start, ':', $feature->seq_region_end, "\tNumber of spanning reads: ", $feature->score, "\n";
                        # We count all the reads that span over a splice site
                        $read_count += $feature->score;
                    }
                    print "\tIntron $index  ", $intron->[1], ':', $intron->[2], "\tNumber of spanning reads: ", $read_count, "\n  Exon ", ++$index, '  ';
                }
                print "\n";
            }
        }
    }
}
