#!/usr/local/bin/perl -w

### gtf2ensembl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::PredictionExon;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;


my ($dbhost,
    $dbname,
    $dbuser,
    $dbpass,
    $dbport,
    $test,
    );

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'test=s'     => \$test,
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

    next if lc($feature) ne "cds";

    $gtf->{'chromosome'}{$chr_name} = 1;
    
    # Put strand into EnsEMBL convention
    if ($strand eq '+') {
      $strand = 1;
    } elsif ($strand eq '-') {
      $strand = -1;
    } else {
      $strand = 0;
    }
    
    # Put the phase into the EnsEMBL convention
    $phase = $gtf_ens_phase{$phase};
    unless (defined $phase) {
      die "Illegal phase in: $line";
    }
    
    my ($transcript) = $group =~ /transcript_id\s*\"([^\"]+)\"/;
    
    if (not $gtf->{'transcript'}{$gene_type}{$transcript}) {
      $gtf->{'transcript'}{$gene_type}{$transcript} = { chromosome => $chr_name };
    }

    push(@{$gtf->{'cds'}{$transcript}}, [$start, $end, $strand, $phase, $score]);
  }
  
  return $gtf;
}


sub prepare_for_writing_genes {
  my ($gtf, $db) = @_;

  my (%ana_hash, %slice_hash);

  # analyses
  my $ana_adap = $db->get_AnalysisAdaptor;
  foreach my $logic_name (keys %{$gtf->{'transcript'}}) {
    if (my $ana = $ana_adap->fetch_by_logic_name($logic_name)) {
      $ana_hash{$logic_name} = $ana;
    } else {
      $ana_hash{$logic_name} = Bio::EnsEMBL::Analysis->new(
                                                           -logic_name      => $logic_name,
                                                           -gff_source      => $logic_name,
                                                           -gff_feature     => 'prediction');
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

  my $pred_tran_adaptor = $db->get_PredictionTranscriptAdaptor;

  # Loop through each type of gene    
  my $counter = 1;
  foreach my $type (sort keys %{$gtf->{'transcript'}}) {
    my $ana = $ana_hash->{$type};
    my $tran_hash = $gtf->{'transcript'}{$type};
        
    # Loop through each gene
    foreach my $tran_name (sort keys %$tran_hash) {
      my $chromosome = $tran_hash->{$tran_name}->{'chromosome'};
      
      my $tsct = Bio::EnsEMBL::PredictionTranscript->new;
      $tsct->analysis($ana);
      $tsct->display_id($tran_name);

      my @cds = sort {$a->[0] <=> $b->[0]} @{$gtf->{'cds'}{$tran_name}};
        
      # Get strand from first exon, and check is
      # the same for all the rest.
      my $strand = $cds[0][2];
      if ($strand == -1) {
        @cds  = reverse @cds;
      }
      
      foreach my $m (@cds) {
        my $start = $m->[0];
        my $end   = $m->[1];
        my $phase = $m->[3];
        my $score = $m->[4];

        my $ex = Bio::EnsEMBL::PredictionExon->new;
        $ex->start($start);
        $ex->end($end);
        $ex->strand($strand);
        $ex->score($score);
        $ex->phase($phase);
        $ex->slice($slice_hash->{$chromosome});
        $tsct->add_Exon($ex);
      }
              
      # test the translation; if it does not translate, we need to guess the phase
      &guess_phases($tsct, $tran_name);
  
      # Write the gene
      eval{
        if ($test) {
          &print_transcript($test, $tsct, $tran_name);
        } else {
          $db->get_PredictionTranscriptAdaptor->store($tsct);
        }
      };
      if ($@) {
        print STDERR "ERROR: Failed to write gene: $@";
        exit(1);
      } else {
        print STDERR "Done Prediction transcript " . $counter++ . " ";
        print STDERR "(" .  $tsct->dbID . ")" if not $test; 
        print STDERR "\n";
      }
    }
  }
}


sub guess_phases {
  my $tran = shift;
  my $tran_name = shift;

  my @exons = @{$tran->get_all_Exons};

  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons;
  }
                                                                                                            
  if ($tran->translate->seq =~ /\*/) {
    my $original_phase = $exons[0]->phase;

    my %phases = (0 => 1, 
                  1 => 1, 
                  2 => 1);
    delete $phases{$original_phase};

    my $found_good = 0;
    foreach my $remaining_phase (keys %phases) {
      my $endphase = $remaining_phase;

      foreach my $exon (@exons) {
        $exon->phase($endphase);
        $endphase = ($exon->phase + $exon->length)%3;
      }
                                                                                                            
      if ($tran->translate->seq !~ /\*/) {
        $found_good = 1;
        last;
      }
    }

    if (not $found_good) {
      # reset the phase to 0 and be done with it
      print "Error: Could not find translation for " . $tran_name . "\n";
      my $endphase = $original_phase;
            
      foreach my $exon (@exons) {
        $exon->phase($endphase);
        $endphase = ($exon->phase + $exon->length) % 3;
      }
    }
    
  }
}


sub print_transcript {
  my ($test, $tran, $name) = @_;
  
  print "# Transcript:\n";
  foreach my $e (@{$tran->get_all_Exons}) {      
    printf("chr: %s, start: %d, end: %d, strand: %d, phase: %d, end_phase: %d\n", $e->slice->seq_region_name, $e->start, $e->end, $e->strand, $e->phase, $e->end_phase); 
  }
  
  my $peptide = Bio::PrimarySeq->new(-id         => $name,
                                     -seq        => $tran->translate->seq,
                                     -moltype    => 'protein' );  
  $test->write_seq($peptide);
}


__END__

