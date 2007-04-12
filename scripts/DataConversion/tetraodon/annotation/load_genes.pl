#!/usr/local/bin/perl -w

### gtf2ensembl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;


my ($dbhost,
    $dbname,
    $dbuser,
    $dbpass,
    $dbport,
    $test,
    $do_stable_ids,
    );

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'test=s'     => \$test,
            'stableids' => \$do_stable_ids,
);


$dbport = 3306 if not $dbport;
$dbuser = "ensro" if not $dbuser;

die "You must supply a database name and host" if not $dbname or not $dbhost;

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);


if ($test) {
  $test = Bio::SeqIO->new(-format => 'fasta',
                          -file => ">$test");
}

my $gtf = &parse_gtf(\*ARGV);
my ($ana_hash, $slice_hash) = &prepare_for_writing_genes($gtf, $db);
&write_ens_genes($gtf, $db, $ana_hash, $slice_hash);


sub parse_gtf {
  my( $fh ) = @_;
  
  # GTF -> EnsEMBL phase convention map
  my %gtf_ens_phase = (
                       0  =>  0,
                       1  =>  2,
                       2  =>  1,
                       '.' => -1,
                       );
  
  # Parse the GTF file, and store its data in the hash $gtf
  my $gtf = {};
  
  while (defined(my $line = <$fh>)) {
    next if $line =~ /^$/;
    next if $line =~ /^#/;
    chomp($line);
    my( $chr_name,
        $gene_type,
        $feature,
        $start,
        $end,
        $score,
        $strand,
        $phase,
        $group,
        ) = split(/\s+/, $line, 9);

    $gtf->{'chromosome'}{$chr_name} = 1;
    
    # Put strand into EnsEMBL convention
    if ($strand eq '+') {
      $strand = 1;
    }
    elsif ($strand eq '-') {
      $strand = -1;
    }
    else {
      $strand = 0;
    }
    
    # Put the phase into the EnsEMBL convention
    $phase = $gtf_ens_phase{$phase};
    unless (defined $phase) {
      die "Illegal phase in: $line";
    }
    
    # Parse the group field
    # (Not technically correct, because can have 1 tag with
    # multiple values.)
    my $tag_val = {};
    foreach my $tv (grep $_, split /\s*;\s*/, $group) {
      if ($tv =~ /(\S+)\s+\"([^\"]+)\"$/) {
        $tag_val->{lc($1)} = $2;
      } elsif ($tv =~ /(\S+)\s+([^\"]+)$/) {
        $tag_val->{lc($1)} = $2;
      } elsif ($tv =~ /(\S+)\s*=\s*\"([^\"]+)\"$/) {
        $tag_val->{lc($1)} = $2;
      } elsif ($tv =~ /(\S+)\s*=\s*([^\"]+)$/) {
        $tag_val->{lc($1)} = $2;
      } else {
        die "Illegal group field: $group\n";
      }
    }        
    
    my ($gene, 
        $transcript,
        $translation,
        $exon
        ) = ($tag_val->{'gene_id'},
             $tag_val->{'transcript_id'},
             $tag_val->{'translation_id'},
             $tag_val->{'exon_id'}
                  );
    
    die "Group field [$group] should contain gene_id and transcript_id\n" 
        if not $gene or not $transcript;
                         
    if ($translation and not $gtf->{'gene'}{$gene_type}{$gene}{$transcript}) {
      $gtf->{'gene'}{$gene_type}{$gene}{$transcript} = { description => "",
                                                         translation => $translation,
                                                         chromosome => $chr_name };
    }

    $feature = lc $feature;
    if ($feature eq 'exon') {
      push(@{$gtf->{'exon'}{$transcript}}, [$start, $end, $strand, $exon]);
    }
    elsif ($feature eq 'cds') {
      push(@{$gtf->{'cds'}{$transcript}}, [$start, $end, $strand, $phase]);
    }
    else {
      warn "Unknown feature '$feature'\n";
    }
  }
  
  return $gtf;
}


sub prepare_for_writing_genes {
  my ($gtf, $db) = @_;

  my (%ana_hash, %slice_hash);

  # analyses
  my $ana_adap = $db->get_AnalysisAdaptor;
  foreach my $logic_name (keys %{$gtf->{'gene'}}) {
    if (my $ana = $ana_adap->fetch_by_logic_name($logic_name)) {
      $ana_hash{$logic_name} = $ana;
    } else {
      $ana_hash{$logic_name} = Bio::EnsEMBL::Analysis->new(
                                                           -logic_name      => $logic_name,
                                                           -gff_source      => $logic_name,
                                                           -gff_feature     => 'gene');
    }
  }

  # slices
  my $slice_adap = $db->get_SliceAdaptor;
  foreach my $chr_name (keys %{$gtf->{'chromosome'}}) {
    $slice_hash{$chr_name} = $slice_adap->fetch_by_region('chromosome', $chr_name);
  }


  return (\%ana_hash, \%slice_hash);
}



sub write_ens_genes {
  my( $gtf, $db, $ana_hash, $slice_hash) = @_;

  my $gene_adaptor = $db->get_GeneAdaptor;

  # Loop through each type of gene    
  my $gene_counter = 1;

  foreach my $type (keys %{$gtf->{'gene'}}) {
    my $ana = $ana_hash->{$type};
    my $gene_hash = $gtf->{'gene'}{$type};
        
    # Loop through each gene
    foreach my $gene_name (keys %$gene_hash) {
      my $tran_hash = $gene_hash->{$gene_name};
      
      my $gene = Bio::EnsEMBL::Gene->new;
      $gene->stable_id($gene_name) if $do_stable_ids;
      $gene->biotype('protein_coding');
      $gene->source("Genoscope");
      $gene->version(1);
      $gene->analysis($ana);


      # Loop through each transcript
      TRAN: foreach my $tran_name (keys %$tran_hash) {
        #printf STDERR "  %20s\t%-s\n", $tran_name, $desc;

        my $chromosome = $tran_hash->{$tran_name}->{'chromosome'};
        my $translation_id = $tran_hash->{$tran_name}->{'translation'};
        
        my $tsct = Bio::EnsEMBL::Transcript->new;
        $tsct->stable_id($tran_name) if $do_stable_ids;                
        $tsct->version(1);
        $tsct->analysis($ana);

        # Get the CDS and exon data for this transcript
        my @mrna = sort {$a->[0] <=> $b->[0]} @{$gtf->{'exon'}{$tran_name}};
        my( @cds );
        if (my $c = $gtf->{'cds'}{$tran_name}) {
          @cds = sort {$a->[0] <=> $b->[0]} @$c;
        }
        
        # Get strand from first exon, and check is
        # the same for all the rest.
        my $strand = $mrna[0][2];
        my @all_ex = (@mrna, @cds);
        for (my $i = 1; $i < @all_ex; $i++) {
          my $this_strand = $all_ex[$i][2];
          if ($strand != $this_strand) {
            print STDERR "ERROR: Multiple strands in '$tran_name'\n";
            next TRAN;
          }
        }
                
        # Flip arrays into translation order if on opposite strands
        if ($strand == -1) {
          @mrna = reverse @mrna;
          @cds  = reverse @cds;
        }
        
        # Make an exon for each entry in the mrna array
        my( @exons );
        foreach my $m (@mrna) {
          my $start = $m->[0];
          my $end   = $m->[1];
          
          my $ex = Bio::EnsEMBL::Exon->new;
          $ex->version(1);
          $ex->stable_id($m->[3]) if $do_stable_ids;
          $ex->start($start);
          $ex->end($end);
          $ex->strand($strand);
          $ex->phase(-1);
          $ex->end_phase(-1);
          $ex->slice($slice_hash->{$chromosome});
          push(@exons, $ex);
        }
                
        # Make Translation and set exon phases if there is CDS info
        if (@cds) {
          my $tsl = Bio::EnsEMBL::Translation->new;
          $tsl->stable_id($translation_id) if $do_stable_ids;
          $tsl->version(1);

          my $j = 0;  # Points to first CDS segement
          my $last_translating_exon = undef;
          for (my $i = 0; $i < @exons; $i++) {
            my $ex = $exons[$i];
            my $cd = $cds[$j];
            
            # Does this overlap the next CDS segment?
            if ($ex->start <= $cd->[1] and $ex->end >= $cd->[0]) {
              
              # Add the phase from the CDS
              if ($strand == -1) {                
                if ($ex->end == $cd->[1]) {
                  $ex->phase($cd->[3]);
                } else {
                  $ex->phase(-1);
                }

                if ($ex->start == $cd->[0]) {
                  my $this_phase = $ex->phase > 0 ? $ex->phase : 0;
                  $ex->end_phase(($this_phase + ($cd->[1] - $cd->[0] + 1))%3);
                } else {
                  $ex->end_phase(-1);
                }

              } else {
                if ($ex->start == $cd->[0]) {
                  $ex->phase($cd->[3]);
                } else {
                  $ex->phase(-1);
                }

                if ($ex->end == $cd->[1]) {
                  my $this_phase = $ex->phase > 0 ? $ex->phase : 0;
                  $ex->end_phase( ($this_phase + ($cd->[1] - $cd->[0] + 1))%3);
                } else {
                  $ex->end_phase(-1);
                }
              }
              #$ex->phase($cd->[3]);
              #$ex->end_phase(($ex->length + $ex->phase)%3);
              
              if (not $tsl->start_Exon) {
                my $pos = $strand == 1 ? $cd->[0] : $cd->[1];
                $tsl->start_Exon($ex);
                $tsl->start(exon_coord($ex, $pos));
              }                

              $last_translating_exon = $ex;
              
              # Move pointer unless this is the last CDS segement
              $j++ unless $j == $#cds;
            } 
          }
                                        
          if ($j != $#cds) {
            my $found = $j + 1;
            my $cds_count = @cds;
            print STDERR "ERROR: didn't match all CDS segements to exons for '$tran_name' (found $found out of $cds_count)\n";
            next TRAN;
          } else {
            my $cd = $cds[$j];
            $tsl->end_Exon($last_translating_exon);
            my $pos = $strand == 1 ? $cd->[1] : $cd->[0];
            $tsl->end(exon_coord($last_translating_exon, $pos));
          }
          
          $tsct->translation($tsl);
        }
        
        # Add each exon to the transcript                
        foreach my $ex (@exons) {
          $tsct->add_Exon($ex);
        }
        
        # Add this transcript to the Gene
        $gene->add_Transcript($tsct);
      }
      
      # Write the gene
      eval{
        if ($test) {
          &print_gene($test, $gene);
        } else {
          $db->get_GeneAdaptor->store($gene);
        }
      };
      if ($@) {
        print STDERR "ERROR: Failed to write gene: $@";
        exit(1);
      } else {
        print STDERR "Written gene " . $gene_counter++ . "\n";
      }
    }
  }
}


# Returns genomic position in exon coordinates
sub exon_coord {
  my( $exon, $coord ) = @_;
  
  if ($exon->strand == 1) {
    my $start = $exon->start;
    return $coord - $start + 1;
  } else {
    my $end = $exon->end;
    return $end - $coord + 1;
  }    
}

sub print_gene {
  my ($test, $gene) = @_;
  
  my $std = $gene->stable_id;
  $std = $gene if not $std;

  foreach my $t (@{$gene->get_all_Transcripts}) {
    my $tid = $t->stable_id;
    $tid = $t if not $tid;

    foreach my $e (@{$t->get_all_Exons}) {
      my $eid = $e->stable_id;
      $eid = $e if not $eid;

      printf("chr: %s, start: %d, end: %d, strand: %d, phase: %d, end_phase: %d, gname: %s, tname: %s, ename: %s\n", 
             $e->slice->seq_region_name, 
             $e->start, 
             $e->end, 
             $e->strand, 
             $e->phase, 
             $e->end_phase, 
             $std, 
             $tid, 
             $eid); 
    }
    
    my $translation = $t->translation;
    printf("  Translation: %s (%d), %s (%d)\n", 
           $translation->start_Exon->stable_id ? $translation->start_Exon->stable_id :$translation->start_Exon, 
           $translation->start,
           $translation->end_Exon->stable_id ? $translation->end_Exon->stable_id : $translation->end_Exon,
           $translation->end);
    
    my $translation_seq = $t->translate;
    
    $test->write_seq($translation_seq);
  }
}


__END__

