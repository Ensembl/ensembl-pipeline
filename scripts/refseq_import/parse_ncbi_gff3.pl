#!/usr/bin/env perl
#
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

 parse_ncbi_gff3.pl - an NCBI GFF3 parser

=head1 SYNOPSIS

 parse_ncbi_gff3.pl [arguments]


=head1 DESCRIPTION

 Code to parse the NCBI RefSeq GFF3 files, create Ensembl
 objects and write to a mysql Ensembl schema database.
 The file uses RefSeq sequence accessions and therefore we
 need them in the seq_region_synonym table for this script
 to work. You can add them using ./load_refseq_synonyms.pl

Example:

     wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Anolis_carolinensis/GFF/ref_AnoCar2.0_top_level.gff3.gz

     gunzip ref_AnoCar2.0_top_level.gff3.gz 

     perl $BASE/scripts/parse_ncbi_gff3.pl --host genebuildn --dbname anolis_carolinensis_otherfeatures_nn_n \
          --user user --pass pass --port 3306--dnadbname anolis_carolinensis_core_nn_n \
          --dnadbhost genebuildn --dnadbuser user --dnadbpass pass -dnadbport 3306 \ 
          --infile ref_AnoCar2.0_top_level.gff3

=head1 AUTHOR

 Daniel Barrell <dbarrell@ebi.ac.uk>, Ensembl genebuild team

=head1 CONTACT

 Please post comments/questions to the Ensembl development list
 <http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
use 5.10.1; 
use feature qw(say);
use Getopt::Long;
use POSIX qw(strftime);
use File::stat;
use File::Basename;

use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use HTML::Entities qw(decode_entities);

use Data::Dumper;

my ($host, $user, $pass, $dbname, $infile, $dnahost, $dnaport, $dnadbname, $dnapass);
my $port = $dnaport = 3306;
my $dnauser = 'ensro';
my $cs = 'toplevel';
my $anal_name = 'refseq_import';
my ($verbose, $test_parse, $ignore_mt) = (0) x 3;
my (%genes, %transcripts, %exons, %CDSs, %xrefs, %missing_sequences, %sequences, %tRNAs,
  %curated_genomic, %microRNAs, %rRNAs, %ig_genes, %ribosomal_slippage ) = () x 13;

&GetOptions (
  'host:s'      => \$host,            'dnahost:s'   => \$dnahost, 
  'port:i'      => \$port,            'dnadbname:s' => \$dnadbname, 
  'user:s'      => \$user,            'dnauser:s'   => \$dnauser,
  'pass:s'      => \$pass,            'dnapass:s'   => \$dnapass,  
  'dbname:s'    => \$dbname,          'dnaport:s'   => \$dnaport,   
  'infile:s'    => \$infile,          'verbose!'    => \$verbose,
  'test_parse!' => \$test_parse,      'ignore_mt!'  => \$ignore_mt,
  'coord_system|cs:s' => \$cs,        'analysis:s'  => \$anal_name,
);
&usage unless defined ($dbname && $host && $user && $infile);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  '-host'     => $host,
  '-port'     => $port,
  '-user'     => $user,
  '-pass'     => $pass,
  '-dbname'   => $dbname) or die('Failed to connect to Ensembl database');

if (!$dnadbname) {
  throw ("You must provide host and dbname details to connect to DNA DB!");
}
my $dna_db =  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host     => $dnahost,
  -user     => $dnauser,
  -port     => $dnaport,
  -dbname   => $dnadbname,
  -pass     => $dnapass,
  -species  => 'dna_'.$db->species());
throw ("Your DNA_DB contains no DNA. You must provide a DNA_DB to connect to") unless (check_if_DB_contains_DNA($dna_db));
$db->dnadb($dna_db);

my $sa = $db->get_SliceAdaptor;
my $ga = $db->get_GeneAdaptor;
my $aa = $db->get_AnalysisAdaptor;
my $analysis = $aa->fetch_by_logic_name($anal_name);
my $timestamp = strftime("%F %X", localtime(stat($infile)->mtime));
my $infile_name = basename($infile);
if ($analysis) {
#    if ($timestamp eq $analysis->db_version) {
#        print STDOUT "\nYour $anal_name is already at the latest version, $timestamp. NO need to load it again\nIf you really want to do it: touch $infile\n and start again.\n\n";
#        exit(0);
#    }
#    else {
  if ($infile_name ne $analysis->db_file) {
    warning('Your old analysis had '.$analysis->db_file." as db_file, it will be update to $infile_name\n");
    $analysis->db_file($infile_name);
  }

  warning('Old db_version: '.$analysis->db_version."\nNew db_version: $timestamp\n");
  $analysis->db_version($timestamp);
  $aa->update($analysis) unless ($test_parse);
#    }
}
else {
    $analysis = Bio::EnsEMBL::Analysis->new(
        -logic_name => $anal_name,
        -db => 'RefSeq',
        -db_version => $timestamp,
        -db_file => $infile_name,
    );
}
# Counts of gene, transcript, exons, non-coding, curated genomic, rRNA, CDS, genes 
my ($gcnt, $tcnt, $ecnt, $ncnt, $cgcnt, $micnt, $rrcnt, $cdscnt, $igcnt, $genes_stored ) = (0) x 10;
my $version = 1;
my $status  = 'PUTATIVE';
my $time    = time();

# Now we know that the file exists, find the mitochondrial accession.
# We don't want to import annotation from the MT. 
my $MT_acc = `grep genome=mitochondrion $infile | cut -f1`; chomp $MT_acc;
say("MT accession: " . $MT_acc) if $MT_acc ne '';
say("!! Ignoring MT sequence !!") if $ignore_mt;

# For species with PAR: fetch regions so that we can skip RefSeq annotations on them.
my @par_regions;
my $par_srid;
my $sql = "SELECT seq_region_id,      "
         ."       seq_region_start,   "
         ."       seq_region_end      "
         ."  FROM assembly_exception  "
         ." WHERE exc_type = 'PAR'    "
         ." ORDER BY seq_region_start ";

my $sth = $db->dbc->prepare($sql);
$sth->execute(); my %row; $sth->bind_columns( \( @row{ @{$sth->{NAME_lc} } } ));
while ($sth->fetch()) {
  push (@par_regions, [$row{seq_region_start}, $row{seq_region_end}]);
  $par_srid = $row{seq_region_id}; # a sensible asumption that there will only be one 
  say("Par region exists: ",$par_srid,":",$row{seq_region_start},"-",$row{seq_region_end});
}


# ================================================================================  
# PARSE : A note about parsing, we try and store as little as possible until we know
#         we need to do something otherwise we are going to end up getting very
#         memory hungry. It is also important that we do not pass around variables
#         between line types: always look up the parent object using the parent_id
#         key if you need to know something about it or you need to inherit some 
#         information from it.
my $fh;
unless (open($fh,'<',$infile)) {
  die("Could not find " . $infile);
}
my @lines = <$fh>;
close($fh);                           

LINE: while (my $line = shift( @lines )) {   
  next LINE if $line =~ /^#/;
  my ($seqname,
      $source,
      $type,
      $start,
      $end,
      $score,
      $strand,
      $phase,
      $attributes) = split(/\t/, $line);
  next LINE if ($seqname eq $MT_acc && $ignore_mt);    
  next LINE unless $type =~ /^(gene|mRNA|transcript|exon|CDS|ncRNA|\w_gene_segment|primary_transcript|rRNA|tRNA)$/;
  next LINE if exists $missing_sequences{$seqname};
  
  # cols 1 and 2 can have HTML encoding which may affect us (indeed, NCBI now use %2C for commas in source)
  $seqname = decode_entities($seqname);
  $source  = decode_entities($source);
  
  # some data common to all lines let through the filter above:
  chomp($attributes); 
  my @atts = split (/;/, $attributes);
  my ($id) = grep s/^ID=(\w+\d+)$/$1/, @atts;

  # On the MT coding genes have CDS only and no exon. This means we can actually parse them as exons
  #$type='exon' if ( $seqname eq $MT_acc && $type eq 'CDS' );

  # SLICE block
  my $slice = $sa->fetch_by_region($cs, $seqname); 
  # It's important to get this right: Ensembl and RefSeq have different styles of handling slices
  # 1) Missing slices. Only report missing slices once. Missing slices are actually slices on
  # alternative assemblies
  if (! defined $slice && ! exists $missing_sequences{$seqname}) {
    $missing_sequences{$seqname}='';
    say("Slice not found " . $seqname);
    next LINE;
  } else {
    say("Slice " . $seqname . " found (" . $slice->name() . ")") unless exists $sequences{$seqname};
    $sequences{$seqname}='';
  }
  # 2) Ignore annotations on pseudo-autosomal regions (currently only X/Y for human)
  if (@par_regions && $slice->get_seq_region_id() == $par_srid) {
    foreach my $aref (@par_regions) {
      if ( ($start >= $$aref[0] && $start <= $$aref[1]) || ($end >= $$aref[0] && $end <= $$aref[1]) ) {
        say("In PAR region, skip...") if $type =~ /^(gene|tRNA|ncRNA|\w_gene_segment|rRNA)$/;
        next LINE;
      }
    }
  } # end SLICE block

  $source =~ s/\,.*$//;  # keep only one source, 'Curated Genomic,Gnomon' does not fit in the DB.

  $strand = -1 if $strand eq '-';
  $strand =  1 if $strand eq '+';
  
  # now deal with the 'type' specific lines, from which, we will make Ensembl objects


  # =============================================================================  
  # GENE : A gene line does not have a biotype so it may be set later. At this
  #        stage we can only really decide for sure if a gene is a pseudogene or
  #        not. When we get to the mRNA/transcript block we will learn more about
  #        this gene. 
  if ($type eq "gene") {

    my $biotype;                                   
    my ($gene_id) = grep s/Dbxref=GeneID:(\d+)/$1/, @atts;
    $gene_id =~ s/,.*$//;
    my ($pseudo)  = grep s/pseudo=(true)/$1/, @atts;

    my $description = "";
    ($description) = grep s/^description=(microRNA).*/$1/, @atts;

    # For the gene symbol first try the name attrib and failing that try the gene one. From what I've
    # seen they have both been the same
    my ($gene_symbol_src1) = grep s/^Name=(.+)/$1/, @atts;
    my ($gene_symbol_src2) = grep s/^gene=(.+)/$1/, @atts;
    if($gene_symbol_src1 ne $gene_symbol_src2) {
      warning("Parsed out two different values for the gene symbol for ".$id."\nVal1: ".
              $gene_symbol_src1."\nVal2: ".$gene_symbol_src2."\nThe script will choose Val1 over Val2 unless Val1 is undef");
    }
    my $gene_symbol;
    if($gene_symbol_src1) {
      $gene_symbol = $gene_symbol_src1;
    } else {
      $gene_symbol = $gene_symbol_src2;
    }

    if($gene_symbol) {
      my $xref = Bio::EnsEMBL::DBEntry->new(
                   -primary_id => $gene_id,
                   -display_id => $gene_symbol,
                   -dbname => 'EntrezGene',
                 );
      $xrefs{$id} = $xref;
    } else {
        warning("Could not find a gene symbol for GeneID: ".$id);
    }
    # make a new gene
    if (! exists($genes{$id})) {
      $gcnt++;
      say($type . " " . $id . " found") if ($verbose);

      if ($pseudo) {
        $biotype = 'pseudogene';

        # Pseudogenes rarely have a transcript line (see the TRANSCRIPT block for an
        # exception) therefore we use the Name as the stable_id for these.
        my ($name) = grep s/^Name=(.*)$/$1/, @atts;

        $transcripts{$id} = {
          id           => $id,
          parent_id    => $id,
          start        => $start,
          end          => $end,
          strand       => $strand,
          biotype      => $biotype,
          stable_id    => $name,
          slice        => $slice,
        };
        $tcnt++;

      } elsif (defined $description && $description eq "microRNA") {
        $biotype="miRNA";

      } else {
        $biotype='protein_coding';
        # It's possible that we find that this is misc_RNA when we come across a
        # transcript line for this gene, don't worry, we will update the gene 
        # biotype as and when this happens.
      }

      # we can make the Ensembl object for genes as we parse the file
      my $gene = Bio::EnsEMBL::Gene->new(
        -START        => $start,
        -END          => $end,
        -STRAND       => $strand,
        -STABLE_ID    => $gene_id,
        -SLICE        => $slice,
        -VERSION      => $version,
        -SOURCE       => $source,
        -STATUS       => $status,
        -ANALYSIS     => $analysis,
        -BIOTYPE      => $biotype,
        -CREATED_DATE => $time,
      );

      $genes{$id}=$gene;

      # set a flag if the source is 'Curated genomic'. These are a little
      # strange: they all seem to be pseudogenes but sometimes have exons and
      # sometimes don't. We are going to skip the single gene lines but try
      # and parse the normal looking pseudogenes.
      if ($source eq "Curated Genomic") {
        say("Found 'Curated Genomic'") if ($verbose);
        $genes{$id}{curated_genomic}='true';
      }

    } else {
      throw("Gene already exists " . $gene_id);
    }
  } # end GENE



  # ============================================================================  
  # TRANSCRIPT : There are two main types, mRNA for protein coding genes whilst
  #              transcript is used for misc_RNA. For some reason NCBI also has
  #              pseudogenes which claim they are misc_RNA in their exons. We 
  #              need to deal with that (example: GeneID:100556170). A note about
  #              transcripts: we cannot make real Ensembl objects until we have
  #              looked at every line in the file because we need the exons and
  #              we don't know their structure yet. Therefore we use our own data
  #              structure for now: %transcripts. Alternatively I could have used
  #              add_Exon(), I suppose, but there was a reason I didn't which I 
  #              haven't documented. Maybe that can change in the future?
  #              There are also variable/rearrangement genes which have multiple
  #              non-overlapping loci within the genome. We skip these for now.
  #              'primary_transcript' denotes trans-splicing microRNAs as far as
  #              I can see.
  elsif ($type =~ /^(mRNA|tRNA|transcript|primary_transcript|rRNA|ncRNA|\w_gene_segment)$/) {

    my $biotype;
    my ($parent_id) = grep s/^Parent=(.*)$/$1/, @atts;
    my ($transcript_id) = grep s/^transcript_id=(.*)$/$1/, @atts;

    # tRNAS/rRNAs on the MT do not have transcript_ids therefore use the parent
    # gene identifier as the transcript stable ID, they are 1:1 gene:transcript.
    # This is also the case for all tRNAscan-SE sources.
    if ($source eq 'tRNAscan-SE') {
      ($transcript_id) = grep s/^Dbxref=GeneID:(\d+).*$/$1/, @atts;
    }
    if ($seqname eq $MT_acc) {
      $transcript_id = $genes{$parent_id}{stable_id};
    }
    if($type =~ /^\w_gene_segment$/) {
      $transcript_id = $genes{$parent_id}{stable_id};
    }

    # Can't see how the code below makes any sense, if $parent_id is defined it doesn't
    # execute so if it does execute it means $parent_id is undef but then it checks to
    # see if an undef key exists in the hash? I'm commenting it out
#    unless (defined $parent_id) {
       # Anole has a tRNA *without* a parent! (maybe we should delete this block of
       # code altogether until tRNAs are implemented?)
#      if (exists $genes{$parent_id}) {
#        delete $genes{$parent_id};
#        $gcnt--;
#      }
#      next TRANSCRIPT; # Seems wrong
#    }


    # biotypes need checking
    if ($type eq "transcript") {
      # some pseudogenes have transcripts, this is NCBIs odd mix of pseudo and misc_RNA
      if ($genes{$parent_id}{biotype} eq "pseudogene") {
        # remove the fake transcript and replace with this one
        delete $transcripts{$parent_id};
        $tcnt--;
        $biotype ="pseudogene";
      } else {
        ($biotype) =  grep s/^gbkey=(.*)$/$1/, @atts;
        # update parent gene biotype since it will now be protein coding.
        $genes{$parent_id}{biotype}=$biotype;
      }

    } elsif ($type eq "mRNA") {
      $biotype = 'protein_coding';
    } elsif ($type eq "primary_transcript") {
      $biotype = $genes{$parent_id}{biotype};
      $microRNAs{$parent_id}=$id;
      $micnt++;
    } elsif ($type eq "rRNA") {
      $biotype = "rRNA";
      $genes{$parent_id}{biotype} = $biotype;
      $rRNAs{$parent_id}=$id;
      $rrcnt++;
    } elsif ($type eq "tRNA") {
      $biotype = "tRNA";
      $genes{$parent_id}{biotype} = $biotype;
      $tRNAs{$parent_id}=$id;
    } elsif ($type eq "ncRNA") {
      ($biotype) = grep s/^ncrna_class=(\w+)$/$1/, @atts;
      $genes{$parent_id}{biotype} = $biotype;
      # found a case where a primary transcript had the transcript ID but the ncRNA
      # does not, therefore use the gene_id instead...
      unless (defined $transcript_id) {
        ($transcript_id) = grep s/^Dbxref=GeneID:(\d+)/$1/, @atts;
        $transcript_id =~ s/,.*$//;
      }
      $ncnt++;
    } elsif($type =~ /^(\w)_gene_segment$/) {
       $biotype = 'IG_'.$1.'_gene';
       $genes{$parent_id}{biotype}=$biotype;
       $ig_genes{$parent_id}=$id;
       $transcript_id = $genes{$parent_id}{stable_id};
       $igcnt++;
    } else {
      throw("Not yet supported: " . $type)
    }

    # make a new transcript
    unless (exists($transcripts{$id})) {
      say($type . " " . $id . " found") if ($verbose);
      $transcripts{$id} = {
        id           => $id,
        parent_id    => $parent_id,
        start        => $start,
        end          => $end,
        strand       => $strand,
        biotype      => $biotype,
        stable_id    => $transcript_id,
        slice        => $slice,
      };
      $tcnt++;
    }
  } # end TRANSCRIPT

  # ============================================================================
  # EXON :  Exons on their own are generally only used for non_coding objects. 
  #         But it's not that easy: occasionally we get a protein coding gene
  #         which has gene, mRNA and exon lines but no CDS (e.g: gene8491).
  #         Having parsed the exons we can find translation location information
  #         in the CDS lines if they exist.
  #         Gene segment lines (IG genes) are treated as exons also. Since IG 
  #         genes have no transcript line we need to fake one up here.  
  elsif ($type =~ /^(exon)$/) {
    my @exon_atts = split (/;/, $attributes);
    my ($exon_id) = grep s/^ID=(.*)$/$1/, @exon_atts;
    my ($parent_id) = grep s/^Parent=(.*)$/$1/, @exon_atts;
    my $biotype;

#    if (my ($segment_type) = $type =~ /^(\w)_gene_segment$/) {
#      $biotype = 'IG_'.$segment_type.'_gene';
      # Fake up a transcript.
#      unless (exists $transcripts{$exon_id}) {
#        say("Faking up transcript for a ", $type);
        # Unfortunately, multiple transcripts in the gene will get the same stable
        # id but that's just how RefSeq model them, without an explicit mRNA line.
#        my ($transcript_id) = grep s/^Dbxref=GeneID:(\d+)/$1/, @exon_atts;
#        $transcript_id =~ s/,.*$//;
        
#        $transcripts{$id} = {
#          id           => $exon_id,
#          parent_id    => $parent_id,
#          start        => $start,
#          end          => $end,
#          strand       => $strand,
#          biotype      => $biotype,
#          stable_id    => $transcript_id,
#          slice        => $slice,
#        };
#        $tcnt++;
        # We can also set the gene biotype to the correct IG_?_gene biotype.
#        $genes{$parent_id}{biotype}=$biotype;
#     }                                                 
      # Finally, make sure we now set the parent_id to be the transcript_id
      # We are dealing with 'exon' lines in the file here - not transcripts!
#      $parent_id = $exon_id;
#    }

    $biotype = $transcripts{$parent_id}{biotype};
    if ($biotype eq "protein_coding") {
      $phase = 0; # If the gene has CDS we will update this later
    } else {
      $phase = -1; # non coding objects always have a phase of -1
    }

    # make a new exon
    $ecnt++;
    say($type . " " . $exon_id . " found") if ($verbose); 

    my $exon = Bio::EnsEMBL::Exon->new(
      -START        => $start,
      -END          => $end,
      -STRAND       => $strand,
      -PHASE        => $phase,
      -END_PHASE    => $phase,
      -STABLE_ID    => $exon_id,
      -BIOTYPE      => $biotype,
      -VERSION      => $version,
      -CREATED_DATE => $time,
      -SLICE        => $slice,
    );
    push(@{$exons{$parent_id}}, $exon);
  } # end EXON


  # ==========================================================================  
  # CDS : For protein coding genes this line contains the phase information
  #       We also need to store the start and end extents so that we can set
  #       the translation when we create transcripts later on.
  elsif ($type eq "CDS") {
    my @cds_atts = split (/;/, $attributes);
    my ($cds_id) = grep s/^ID=(.*)$/$1/, @cds_atts;
    my ($parent_id) = grep s/^Parent=(.*)$/$1/, @cds_atts;
    # and we will use this as the translation.stable_id...
    my ($protein_id) = grep s/protein_id=(.*)$/$1/, @cds_atts;

    # skip if exception=ribosomal slippage
    my ($cds_exception) = grep s/^exception=(.*)$/$1/, @cds_atts;
    # Note: some species, e.g. Anole, have an exception that can be longer than solely
    # ribosomal slippage. e.g: 'exception=ribosomal slippage%2C annotated by transcript or proteomic data'
    if (defined $cds_exception && $cds_exception =~ /^ribosomal\sslippage/x) {
      say("Found 'Ribosomal Slippage'") if ($verbose);
      # for now - we need to delete these transcripts, we cannot store
      # the multiple translations in Ensembl
      $ribosomal_slippage{$parent_id}=$parent_id;    
      next LINE;
    }

    # generate end_phase
    my $exon_length = ( $end - $start ) + 1; 
    my $end_phase = ( $exon_length + $phase ) % 3;

    # make a new CDS, we can use exon objects
    $cdscnt++;
    say($type . " " . $cds_id . " found") if ($verbose); 

    my $cds = Bio::EnsEMBL::Exon->new(
      -START        => $start,
      -END          => $end,
      -PHASE        => $phase,
      -END_PHASE    => $end_phase,
      -STABLE_ID    => $protein_id,
    );
    push(@{$CDSs{$parent_id}}, $cds);


    # For MT coding genes the file does not provide an exon or transcript object
    # Fortunately these are all 1:1 gene:CDS so we can fake both up here:
    if ($seqname eq $MT_acc) {
      say("Faking up transcript and exon for protein coding MT CDS: " . $protein_id);
      my $biotype = 'protein_coding';
      my $exon = Bio::EnsEMBL::Exon->new(
        -START        => $start,        
        -END          => $end,
        -STRAND       => $strand,
        -PHASE        => $phase,
        -END_PHASE    => $phase,
        -STABLE_ID    => $protein_id,   
        -BIOTYPE      => $biotype,
        -VERSION      => $version,
        -CREATED_DATE => $time,
        -SLICE        => $slice,
      );
      $ecnt++;
      push(@{$exons{$parent_id}}, $exon);  
                         
      $transcripts{$id} = {
        id           => $parent_id,
        parent_id    => $parent_id,
        start        => $start,
        end          => $end,
        strand       => $strand,
        biotype      => $biotype,
        stable_id    => $genes{$parent_id}{stable_id},
        slice        => $slice,
      };
      $tcnt++;
    }

  } # end CDS 


  
  # ============================================================================
  else {
    throw("I'll never get here!");
  }
} # end LINE

# ==============================================================================  
# Now deal with adding the exons, CDSs and translations to new transcript objects
# and then these transcripts to genes.
# This is where we finally make real Ensembl transcript and translation objects.

TRANSCRIPT: foreach my $k (keys %transcripts) {  
  # first check that we have exons. CAVEAT: In the case of 'Curated Genomic'
  # lines it seems that if it's a single exon gene it sometimes has an exon line,
  # but sometimes does not. All these lines in rat were pseudogenes. I am skipping
  # the gene lines without exons until I can find better documentation on them.
  my $parent_id = $transcripts{$k}{parent_id};
  if (exists $genes{$parent_id}{curated_genomic} && ! exists $exons{$k}) {
    say("Skipping curated_genomic transcript with no exons: " . $k) if $verbose;
    $curated_genomic{$parent_id}=$parent_id;
    next TRANSCRIPT;
  }

  # if this transcript has ribosomal slippage we want to delete it
  if (exists $ribosomal_slippage{$k}) {
    say("Deleting ribosomal slippage transcript: " . $k) if ($verbose);
    delete $transcripts{$k};
    $tcnt--;
    # if that was the only transcript for the gene then delete the gene also...
    my $found;
    RS: foreach my $kk (keys %transcripts) {
      if ($transcripts{$kk}{parent_id} eq $parent_id) {
        $found = 'true';
        last RS;
      }
    }
    unless ($found) {
      say("This ribosomal slippage transcript is the only one available for gene " .
        $parent_id . " => delete it.") if ($verbose);
      delete $genes{$parent_id};
      $gcnt--;
    }
    next TRANSCRIPT;
  }

  # Are there exons?
  if (exists $exons{$transcripts{$k}{id}}) {
    # Pseudogenes on the negative strand have their exons defined backwards,
    # in essence, we flip them here BUT we use the same logic for all types
    # of exons to ensure we have them in the correct rank therefore we avoid
    # making any assumptions of the file order.
    my @exons = @{$exons{$transcripts{$k}{id}}}; 
    if ($exons[0]->strand() == 1) {
      @exons = sort {$a->start() <=> $b->start()} @exons;
    } else {
      @exons = sort {$b->start() <=> $a->start()} @exons;
    } 
    # retrospective note: file order is probably now assumed for this parser :-/


    # TRANSLATION BLOCK - by far the most complicated block, therefore verbosely commented
    #                     plus I've left a lot of the debugging lines in since it is very
    #                     likely you'll have to add new code as more weird and wonderful
    #                     translations appear in the GFF3 files.

    my $translation;
    # Use this to look up the correct CDS array:
    my $transcript_id = $transcripts{$k}{id};
    my $transcript_biotype = $transcripts{$k}{biotype};

    # This code might not be needed in future if gene segment becomes a parent of CDS
#    if($transcript_biotype =~ /^IG_[CDJV]_gene$/) {
#      my $transcript_parent_id = $transcripts{$k}{parent_id};
#      say "Found IG biotype, checking for CDS associated with parent gene (".$transcript_parent_id.")";
#      if(exists $CDSs{$transcript_parent_id}) {
#        $CDSs{$transcript_id} = $CDSs{$transcript_parent_id};
#        say "Set CDS parent to: ".$transcript_parent_id;
#      } else {
#        say "No CDS present"
#      }
#    }

    # If we have CDSs we need to generate a translation
    if (exists $CDSs{$transcript_id}) {
      # These are the data we need to create a translation object
      my ($start_exon, $end_exon, $start_rank, $end_rank, $seq_start, $seq_end, $start_cds, $end_cds, $tr_stable_id);
      # sort the CDSs, in the same way as we did the exons
      my @CDS = @{$CDSs{$transcript_id}};
      if ($exons[0]->strand() == 1) {
        @CDS = sort {$a->start() <=> $b->start()} @CDS;
        $start_cds = ${$CDSs{$transcript_id}}[0];
        $end_cds   = ${$CDSs{$transcript_id}}[-1];
      } else {
        @CDS = sort {$b->start() <=> $a->start()} @CDS;
        $start_cds = ${$CDSs{$transcript_id}}[-1];
        $end_cds   = ${$CDSs{$transcript_id}}[0];
      }
      # found a case (human PEG10, GeneID:23089) where the number of CDS > number
      # of exons, i.e >1 CDS per exon. For example:
      #
      #     exon1           intron       exon2
      #
      #     ==========######-------------######|######=========
      #     
      #     5'UTR     {CDS1}             {CDS2}|{CDS3}    3'UTR
      #
      # The second case which breaks the parser is where a single exon has > 1 CDS
      # like this:
      #
      #     ================------------====######|#####=======
      #
      # We would need to store multiple translations here and at the moment we cannot do this.
      # These are cases of ribosomal slippage in the CDS: exception=ribosomal slippage
      # therefore we will deal with these in the CDS, not here 
      my $num_exons = scalar(@exons);
      my $num_CDS   = scalar(@CDS);
      if ($num_CDS > $num_exons) {
        warning("More CDSs than number of exons! Transcript id: ". $transcript_id);
        delete $genes{$parent_id};
        $gcnt--;
        next TRANSCRIPT;
      }

      # Loop through exons and find which is the start and end exon for the 
      # translation. Whilst we are doing this we also need to update the 
      # exon phase information. Single exon coding transcripts can be dealt
      # with separately. Single CDS in multi-exon transcripts complicate 
      # matters here so expect a lot of edge cases to test.
      if ($#exons > 0) {
        # multi exon
        say("MULTI EXON (" . $num_exons . ") with " . $num_CDS . " CDS") if ($verbose);
        EXON_EXTENT: for my $i (0 .. $#exons) {
          say("i = " . $i . " cds start: " . $CDS[0]->start() .
              " exon extent: " . $exons[$i]->start() . "-" .
              $exons[$i]->end()) if ($verbose);
          if ( $exons[$i]->end() >= $start_cds->start() && $exons[$i]->end() <= $start_cds->end() ) {
            say("IN search for START") if ($verbose);
            $start_exon = $exons[$i];
            $start_rank = $i;
            $exons[$i]->phase($start_cds->phase());
            $exons[$i]->end_phase($start_cds->end_phase());
            $tr_stable_id = $start_cds->stable_id(); # only need do this once
            # If this is a single CDS model we need to compute the end_exon here
            if ($num_CDS == 1) {
              say("SINGLE CDS in MULTI EXON transcript") if ($verbose);
              $end_exon = $start_exon;
              $end_rank = $start_rank;
              last EXON_EXTENT;
            }
          }
          elsif ( $exons[$i]->start() >= $end_cds->start() && $exons[$i]->start() <= $end_cds->end() ) {
            say("IN search for END") if ($verbose);
            $end_exon = $exons[$i];
            $end_rank = $i;
            $exons[$i]->phase($end_cds->phase());
            $exons[$i]->end_phase($end_cds->end_phase());
            $tr_stable_id = $start_cds->stable_id(); # only need do this once
            if ($num_CDS == 1) {
              say("SINGLE CDS in MULTI EXON transcript") if ($verbose);
              $start_exon = $end_exon;
              $start_rank = $end_rank;
              last EXON_EXTENT;
            }
          }
          elsif ($num_CDS == 1 &&
                 $CDS[0]->start() >= $exons[$i]->start() &&
                 $CDS[0]->start() < $exons[$i]->end()) {
            # found single cds that lies within one exon of a multi exon transcript
            say("SINGLE CDS in ONE EXON of MULTI EXON transcript") if ($verbose);
            $start_exon = $end_exon = $exons[$i];
            $start_rank = $end_rank = $i;
            $exons[$i]->phase($end_cds->phase());
            $exons[$i]->end_phase($end_cds->end_phase());
            $tr_stable_id = $start_cds->stable_id(); # only need do this once
            last EXON_EXTENT;
          }
          else {
            #TODO: see if we can move TRANSLATEABLE_EXON: loop logic below to this else block
            # it would mean thinking about if we were in a coding section or not.....
            # ... and as these different edge cases for multi-exon/single-cds keep cropping up
            # it becomes more difficult to do....
            say("Not a start or an end coding exon.") if ($verbose);
            next EXON_EXTENT;
          }
        } # EXON_EXTENT    
      } else {
        # single exon
        say("SINGLE EXON and CODING") if ($verbose);
        $start_exon = $exons[0];
        $start_rank = 0;
        $end_exon = $exons[0];
        $end_rank = 0;
        $exons[0]->phase($start_cds->phase());
        $exons[0]->end_phase($start_cds->end_phase());
        $tr_stable_id = $start_cds->stable_id();
      }

      # now we know the start and end of the coding exons, update the
      # translatable exons with the phase information from the CDSs.
      # Note that the CDS indices are still 0..-1 so they need their
      # own index:
      say("We have been looking at the coding region for ",$k) if ($verbose);
      unless (defined($start_rank) && defined($end_rank)) {
        warning("Found a case where we cannot compute the translation, you will need to write more code to deal with this case: ",$k);
        delete $genes{$parent_id};
        $gcnt--;
        next TRANSCRIPT;
        # Here's the entry I was getting problems with:
        # NC_000014.9	Curated Genomic	gene	106037902	106038345	.	-	.	ID=gene28742;Name=IGHV2-5;Dbxref=GeneID:28457,HGNC:5576,IMGT%2FGENE-DB:IGHV2-5;description=immunoglobulin heavy variable 2-5;gbkey=Gene;gene=IGHV2-5;gene_synonym=IGHV25,VH
        # NC_000014.9	Curated Genomic	CDS	106038300	106038345	.	-	0	ID=cds46977;Parent=gene28742;Note=2-5;Dbxref=GeneID:28457,HGNC:5576,IMGT%2FGENE-DB:IGHV2-5;exception=rearrangement required for product;gbkey=CDS;gene=IGHV2-5
        # NC_000014.9	Curated Genomic	CDS	106037902	106038213	.	-	2	ID=cds46977;Parent=gene28742;Note=2-5;Dbxref=GeneID:28457,HGNC:5576,IMGT%2FGENE-DB:IGHV2-5;exception=rearrangement required for product;gbkey=CDS;gene=IGHV2-5
        # NC_000014.9	Curated Genomic	V_gene_segment	106038300	106038345	.	-	.	ID=id689865;Parent=gene28742;Dbxref=GeneID:28457,HGNC:5576,IMGT%2FGENE-DB:IGHV2-5;gbkey=V_segment;gene=IGHV2-5;part=1%2F2
        # NC_000014.9	Curated Genomic	V_gene_segment	106037902	106038223	.	-	.	ID=id689865;Parent=gene28742;Dbxref=GeneID:28457,HGNC:5576,IMGT%2FGENE-DB:IGHV2-5;gbkey=V_segment;gene=IGHV2-5;part=2%2F2
        # This looks like:
        #                    3'            <----          5'
        # gene                =============================
        # CDS                 ======                 ======
        # V_gene_segment      =======                ======
        #                          ^
        # Emailed Terence Murphy at NCBI - it's a mistake. This also:
        #
        #                    3'            <----          5'
        # gene                =============================
        # CDS                 ======                 ======
        # V_gene_segment      =============================
      }
      my $CDS_index = 0;
      TRANSLATEABLE_EXON: for my $i ($start_rank .. $end_rank) {
        $exons[$i]->phase($CDS[$CDS_index]->phase());
        $exons[$i]->end_phase($CDS[$CDS_index]->end_phase());
        $CDS_index++;
      } # TRANSLATEABLE_EXON
      # adjust for strandedness         
      if ($transcripts{$k}{strand} eq '1') {
        say("computing for plus strand") if ($verbose);
        $seq_start = ( $start_cds->start() - $start_exon->start() ) + 1;
        $seq_end = ( $end_cds->end() - $end_exon->start() ) + 1;
      } else {
        # swap start and end exon        
        say("computing for minus strand") if ($verbose);
        $seq_start = ($end_exon->end() - $end_cds->end()) + 1;
        $seq_end = $start_exon->end() - $start_cds->start() + 1;
        my $tmp = $start_exon; $start_exon = $end_exon; $end_exon = $tmp;
      }
      if ($verbose) {
        say("SEQ_START: " . $seq_start);
        say("SEQ_END:   " . $seq_end);
        say("About to make translation");
      }
      $translation = Bio::EnsEMBL::Translation->new(
        -START_EXON => $start_exon,
        -END_EXON   => $end_exon,
        -SEQ_START  => $seq_start,
        -SEQ_END    => $seq_end,
        -STABLE_ID  => $tr_stable_id,
      );
      # some genes can have miscRNA transcripts AND CDSs. In these cases the gene
      # biotype will have been set to 'miscRNA'. Fix them here.
      unless ($genes{$parent_id}->biotype() eq 'protein_coding') {
        $genes{$parent_id}->biotype('protein_coding'); 
      }
    } # TRANSLATION BLOCK

    # inherit the parent gene's source
    my $source = $genes{$parent_id}{source};

    say("\n\n\nAbout to make transcript") if ($verbose);

    my $transcript = Bio::EnsEMBL::Transcript->new(
      -START        => $transcripts{$k}{start},
      -END          => $transcripts{$k}{end},
      -STRAND       => $transcripts{$k}{strand},
      -STABLE_ID    => $transcripts{$k}{stable_id},
      -SLICE        => $transcripts{$k}{slice},
      -BIOTYPE      => $transcripts{$k}{biotype},
      -VERSION      => $version,
      -SOURCE       => $source,
      -STATUS       => $status,      
      -ANALYSIS     => $analysis,
      -CREATED_DATE => $time,
      -EXONS        => \@exons,
    );

    say("Check if gene parent found") if ($verbose);
    if (exists $genes{$parent_id}) {
      my $stabid = $transcript->stable_id();
      say("Adding transcript " . $stabid . " to gene " . $parent_id) if ($verbose);
      $genes{$parent_id}->add_Transcript($transcript);
    } else {
      # this should never happen (edit: so why isn't it thrown?)
      say("Parent Gene not found for transcript: " . $transcripts{$k}{stable_id});
      next TRANSCRIPT;
      $tcnt--;
    }
    # finally add the translation
    if (defined $translation) {
      $transcript->translation($translation);
    } else {
      say("No translation for transcript " . $transcript->stable_id());
    }
    # and we are done.

  } else {
    # This should never be called
    warning("Transcript " . $transcripts{$k}{id} . " must have at least one exon. Deleting transcript.");
    delete $transcripts{$k};
    $tcnt--;
  }
} # end TRANSCRIPT



# ==============================================================================  
# store objects to the database
unless ($test_parse) {
  my $dbea = $db->get_DBEntryAdaptor;
  STORE: for my $k (keys %genes) {
    my $gene = $genes{$k};
    unless(ref($gene) eq "Bio::EnsEMBL::Gene") {
     warning("Issue with the following gene id: ".$k."\nNot a Gene object, skipping. Dump:\n".Dumper($gene));
     next;
    }

    say("gene start:end:strand is ". join( ":", map { $gene->$_ } qw(start end strand) )) if ($verbose);
    if ($gene->get_all_Transcripts()) {
      if($xrefs{$k}) {
         $gene->add_DBEntry($xrefs{$k});
         $gene->display_xref($xrefs{$k});
      }
      # code to handle HAP/PATCH regions, i.e. transform to toplevel
      $gene = $gene->transform('toplevel') unless $gene->slice()->is_toplevel();
      #fix_duplicate_transcript_stable_ids($gene);
      $ga->store($gene, undef, undef, 1);
      $genes_stored++;
    } else {
      warning("Failed to save gene on " . $gene->seqname() . " @ " . $gene->seq_region_start());
      $gcnt--;
    }
  } # end STORE
}

#TODO: better final report
say("------------------------------------------");
say($gcnt . " genes parsed");
say($tcnt . " transcripts parsed");
say($ecnt . " exons parsed");
say($ncnt . " non coding genes parsed");
say($cgcnt . " Curated Genomic genes parsed");
say($rrcnt . " rRNA genes parsed");
say(scalar(keys %tRNAs) . " tRNA genes parsed");
say("------------------------------------------");
say($genes_stored, " genes stored");

exit(0);


# Subs                      

sub check_if_DB_contains_DNA {
  my ($dba) = @_;
  my $sql_command = "select count(*) from dna";
  my $dnasth = $dba->dbc->prepare($sql_command);
  $dnasth->execute();
  my @dna_array = $dnasth->fetchrow_array;
  if ( $dna_array[0] > 0 ) {
    say("Your DNA_DB " . $dba->dbc->dbname . " contains DNA sequences.") if ( $verbose );
    return 1;
  } else {
    say("Your DNA_DB " . $dba->dbc->dbname . " does not contain DNA sequences.") if ( $verbose );
    return 0;
  }
}

sub fix_duplicate_transcript_stable_ids {
  # Sometimes the stable id for transcripts in a gene are identical. For example if they
  # used the gene stable id as their stable id (gene_segments). Sometimes it can be a mix
  # of unique and duplicated entries (miRNAs did this). This subroutine will take a gene
  # and version any non-unique transcript stable ids

  # NOTE, need to change from using db ids as they don't exist yet
  my $gene = @_;
  my @transcripts = @{$gene->get_all_Transcripts()};
  my $transcript_count = scalar(@transcripts);
  unless($transcript_count > 1) {
    next;
  }

  my %transcript_id_hash = ();
  my $dup_count = 0;
  foreach my $transcript (@transcripts) {
    my $transcript_db_id = $transcript->dbID();
    my $transcript_sid = $transcript->stable_id;
    unless(defined($transcript_id_hash{$transcript_sid})) {
      $transcript_id_hash{$transcript_sid} = [];
    }
    if($transcript_id_hash{$transcript_sid}) {
      $dup_count++;
    }
    push(@{$transcript_id_hash{$transcript_sid}},$transcript_db_id);
  }

  if($dup_count) {
    foreach my $transcript_sid (keys(%transcript_id_hash)) {
      my @sid_array = @{$transcript_id_hash{$transcript_sid}};
      unless(scalar(@sid_array) > 1) {
        next;
      }
      for(my $i=0; $i<scalar(@sid_array); $i++) {
        say "Warning: updating duplicate transcript stable id, set stable_id='".$transcript_sid.
            "_".($i+1)."' where transcript_id=".$sid_array[$i]." and gene_id=".$gene->dbID;
        
      }
    }
  }
}

sub usage {
say("

Example RefSeq gff3 file for parsing:

     ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38_top_level.gff3.gz


Example usage: 

     perl parse_ncbi_gff3.pl -dbhost host -dbuser user -dbpass *** -dbname dbname \
       -dbport 3306 -dnadbname dnadbname -dnadbhost dnadbhost -dnadbuser dnadbuser \
       -dnadbport 3306 -dnadbpass dnadbpass \
       -infile ref_GRCh38_top_level.gff3 -test_parse


Script options:

    -dbname       Database name

    -host       Database host

    -port       Database port

    -user       Database user

    -pass       Database password

    -dnadbname    DNA Database name

    -dnahost    DNA Database host

    -dnauser    DNA Database user
    
    -dnaport    DNA Database port
    
    -dnapass    DNA Database pass

    -infile       RefSeq gff3 file


Optional arguments:

    -coord_system_name    Coordinate system to use

    -verbose              Increase verbosity of output messages

    -test_parse           Parse the file but do not store anything in the database

    -ignore_mt            Do not parse lines belonging to the mitochondrial sequence    
    
");
exit;    
}
