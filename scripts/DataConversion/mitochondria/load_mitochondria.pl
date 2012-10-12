#!/usr/local/ensembl/bin/perl -w

=head1 Synopsis

load_mitochondria.pl 

=head1 Description

Parses mitochondrial genes out of genbank format NC_ file  and writes them to the database specified.
Can also load the sequence of the chromosome, the corresponding assembly entries and the appropriate attribute types if needed.

=head1 Config

All configuration is done through MitConf.pm

=cut

use strict;
use MitConf;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::SeqEdit;
use Bio::EnsEMBL::Pipeline::Tools::GeneBankParser;


use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace throw warning);

use constant {
    INFO_TYPE         => 'DIRECT',
    INFO_TEXT         => 'Imported from GeneBank',
    REFSEQ_PEP        => 'RefSeq_peptide',
    REFSEQ_XPEP       => 'RefSeq_peptide_predicted',
    STATUS            => 'KNOWN',
    CHECK_CODON_TABLE => 2
};

my %EXTERNAL_DB = ( 'GeneID'               => 'EntrezGene',
                    'UniProtKB/Swiss-Prot' => 'Uniprot/SWISSPROT',
                  );
my $help;
my @genes;;
my $translation; # JUST INCASE
my $count;
#use scaffold or supercontig:
my $scaffold = 'scaffold';
my %opt;
my $path = undef;

 # options submitted with commandline override MitConf.pm 

GetOptions( \%opt,
            '-h|help',
            'dbhost=s',
            'dbuser=s',
            'dbpass=s',
            'dbport=i',
            'dbname=s',
            'contig=s',         # MIT_CONTIG_SEQNAME
            'chromosome=s',     # MIT_CHROMOSOME_SEQNAME
            'scaffold=s',    # MIT_SCAFFOLD_SEQNAME
            'clone=s',          # MIT_CLONE_SEQNAME
            'toplevel=s',       # MIT_TOPLEVEL
            'gene_type=s',
            'trna_type=s',
            'rrna_type=s',
            'codon_table=i',
            'name=s',
            'genbank_file=s',
            'noxref!',
            'path=s', );
            # or &usage();

if ( $opt{path}){
  print "you specify path: ", $opt{path},"\n";
  $path =  $opt{path};
}

if (    $opt{dbhost} && $opt{dbuser}
     && $opt{dbname} && $opt{dbpass}
     && $opt{dbport} ) {
  $MIT_DBHOST = $opt{dbhost};
  $MIT_DBUSER = $opt{dbuser};
  $MIT_DBPASS = $opt{dbpass};
  $MIT_DBPORT = $opt{dbport};
  $MIT_DBNAME = $opt{dbname};
}

$MIT_GENBANK_FILE = $opt{genbank_file} if $opt{genbank_file};
$MIT_LOGIC_NAME   = $opt{logic_name}   if $opt{logic_name};
$MIT_NAME         = $opt{name}         if $opt{name};
$MIT_TOPLEVEL     = $opt{toplevel}     if $opt{toplevel};
$MIT_CODON_TABLE  = $opt{codon_table}  if $opt{codon_table};
$MIT_GENE_TYPE    = $opt{gene_type}    if $opt{gene_type};
$MIT_TRNA_TYPE    = $opt{trna_type}    if $opt{trna_type};
$MIT_RRNA_TYPE    = $opt{rrna_type}    if $opt{rrna_type};

unless (    $MIT_DBHOST && $MIT_DBUSER
         && $MIT_DBNAME && $MIT_GENBANK_FILE
         && !$help ) {
  warning( "Can't run without MitConf.pm values: "
      . "MIT_DBHOST $MIT_DBHOST "
      . "MIT_DBUSER $MIT_DBUSER "
      . "MIT_DBNAME $MIT_DBNAME "
      . "MIT_DBPASS $MIT_DBPASS "
      . "MIT_GENBANK_FILE $MIT_GENBANK_FILE "
      . "MIT_LOGIC_NAME $MIT_LOGIC_NAME "
      . "MIT_NAME $MIT_NAME "
      . "MIT_TOPLEVEL $MIT_TOPLEVEL "
      . "MIT_CODON_TABLE  $MIT_CODON_TABLE "
      . "MIT_GENE_TYPE $MIT_GENE_TYPE "
      . "MIT_TRNA_TYPE $MIT_TRNA_TYPE "
      . "MIT_RRNA_TYPE $MIT_RRNA_TYPE " );
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}
############################
#PARSE GENBANK FILE FIRST
open(my $fh, $MIT_GENBANK_FILE) || die("Could not open GeneBank file $MIT_GENBANK_FILE\n");
my $genebank = Bio::EnsEMBL::Pipeline::Tools::GeneBankParser->new($fh);
my $genebank_hash = $genebank->next_entry;
close($fh);
my $features = $genebank_hash->{-features};

####################
#Open db for writing

my $output_db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( '-host'     => $MIT_DBHOST,
                                      '-user'     => $MIT_DBUSER,
                                      '-pass'     => $MIT_DBPASS,
                                      '-dbname'   => $MIT_DBNAME,
                                      '-port'     => $MIT_DBPORT,
                                      '-no_cache' => 1, );

########################
# Check taxon is correct

my $meta_container = $output_db->get_MetaContainer();
my $db_taxon_id = $meta_container->get_taxonomy_id(); 
my $gb_taxon_id; 
foreach my $feature (@$features) {
    if (exists $feature->{source}) {
        $gb_taxon_id = $feature->{source}->{db_xref}->{taxon};
        last;
    }
}

if ($db_taxon_id != $gb_taxon_id) {
    throw("Your taxon ids differ: $db_taxon_id and $gb_taxon_id\n");
}

###########################################
# write chromosome if it is not already there

my $dbe_adaptor = $output_db->get_DBEntryAdaptor;
my $slice_adaptor = $output_db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('toplevel', $MIT_NAME);
my $slices_ref =  &get_chromosomes($genebank_hash, $output_db);
my %slices = %{$slices_ref};

if ($slice){
  print "Found chromosome ".$slice->name."\n" ;
}  else {
  print "There is no chromosome in $MIT_DBNAME called $MIT_NAME\nHave found the following slices to load:\n";
  foreach my $cs (keys %slices){
    print "$cs    \t".$slices{$cs}->name."\n";
  }
  print "Do you want to load the chromosome and "
      . "the associated assembly entries from the genbank file?(Y/N) ";
  my $answer = <>;
  chomp $answer;
  if ($answer eq "y" or $answer eq "Y"){
    &load_chromosomes(\%slices, $output_db, $genebank_hash->{-seq} );
    $slice = $slice_adaptor->fetch_by_region('toplevel',$MIT_NAME);
  }
  else {
    print "Ok gene load aborted.\n";
    exit 0;
  }
}

#########################################
# Check that there is an entry in seq_region_attrib 
# for the MT chromosome. It needs to use the codon table 2
my $has_correct_codon_table = check_codon_table($output_db,$slice);
if ($has_correct_codon_table) {
  print "Using codon table 2 for translations. (This is correct.)\n";
} else {
  throw(   "Cannot find seq_region_attrib entry for MT chromosome. "
         . "Need to specify value=2 for codon_table." );
}


#########################################
#Check that there are not genes in the MT chromosome before uploading the new ones.

my @mt_genes = @{$slice->get_all_Genes}; 

if (@mt_genes && scalar(@mt_genes) > 0){
  print "There are already ",scalar(@mt_genes)," genes in $MIT_NAME chromosome\n";

  print "Do you want to remove them?(Y/N)\n";
  my $g_answer = <>;
  chomp $g_answer;
  if ($g_answer eq "y" or $g_answer eq "Y"){
    my $gene_adaptor = $output_db->get_GeneAdaptor;
    foreach my $mt_gene(@mt_genes){  
      print "Removing gene: ",$mt_gene->dbID,"\n";
      $gene_adaptor->remove($mt_gene);
    }
    print "Genes removed succesfully, moving to new genes load\n";
  }else{
    print "You choose not to remove the genes\n"
        . "Do you want to keep loading the MT genes? "
        . "(This may create duplicated entries)(Y/N)?\n";
    my $load_answer = <>;
    chomp $load_answer;
    if ($load_answer eq "y" or $load_answer eq "Y"){
      print "Loading genes without removing existing ones, Thanks\n";
    }else{
      print "You choose to abort gene loading. Program exiting\n";
      exit 0;
    }
  }
}
  
#########################################
#Fetch analysis object

my $logic_name;

if ($MIT_LOGIC_NAME){
  $logic_name = $MIT_LOGIC_NAME;
} else {
  $logic_name = 'MT_genbank_import';
  print "Cannot find MIT_LOGIC_NAME - using standard logic name ensembl\n" if $MIT_DEBUG;
}
if (!$MIT_DB_FILE) {
  warning("MIT_DB_FILE not defined");
}
if (!$MIT_DB_VERSION) {
  warning("DB_VERSION not defined");
}

my $ensembl_analysis = $output_db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
if(!defined $ensembl_analysis){
  #croak "analysis $logic_name not found\n";
  $ensembl_analysis =
    Bio::EnsEMBL::Analysis->new( -logic_name => 'MT_genbank_import', -db_version => $MIT_DB_VERSION );
  if (defined $MIT_DB_FILE) {
    $ensembl_analysis->db_file($MIT_DB_FILE);
  }
  print "You have no ".$logic_name." defined creating new object\n";
}


 #########################################################
 # Create gene objects from array of hashes
 # contining parsed coords
 # index 0 contains the accession 
 # index 1 contains source data (taxon and sequence file)
 # subsequent indecies contain all the other annotations 
 # 

my $length = scalar(@$features);
my $entry_count = 1;
for (my $index = 1; $index < $length; $index++) {
    my $feature = $features->[$index];
    my ($key) = keys %$feature;
    next unless ($key eq 'CDS' or $key =~ /[tr]RNA/);
    printf "ENTRY %d\n========\n", $entry_count;
    if (@{$feature->{$key}->{exons}} > 1) {
        if ($feature->{$key}->{note} =~ /frameshift/i) {
            warning('There is a frameshift in '.$feature->{$key}->{product}.", the gene has more than one exon!\n");
        }
        else {
            throw($feature->{$key}->{product}." has more than one exon!\n");
        }
    }
    if (exists $feature->{$key}->{note} and $feature->{$key}->{note} =~ /tRNAscan-SE/) {
        warning('Skipping '.$feature->{$key}->{product}.' : '.$feature->{$key}->{note}."\n");
        next;
    }
    #############
    # TRANSCRIPTS

    my $transcript = new Bio::EnsEMBL::Transcript;
    my $exon_number = 0;
    my $start_exon;
    my $end_exon;
    foreach my $exon_hash (@{$feature->{$key}->{exons}}) {
        my $exon = new Bio::EnsEMBL::Exon;
        $exon->start($exon_hash->{start}); 
        $exon->end($exon_hash->{end});
        $exon->strand($exon_hash->{strand});
        $exon->slice($slice);
        $exon->phase(0);
        $exon->end_phase(($exon->end - $exon->start + 1)%3);
        $transcript->add_Exon($exon);
        $start_exon = $exon if ($exon_number == 0);
        $end_exon = $exon;
        ++$exon_number;
    }
    $transcript->start_Exon($start_exon);
    $transcript->end_Exon  ( $end_exon );
    my $type;
    my $gene = new Bio::EnsEMBL::Gene;
    if ($key =~ /(\w)RNA/) {
        if ($1 eq 't') {
            $type = $MIT_TRNA_TYPE;
        }
        elsif ($1 eq 'r') {
            $type = $MIT_RRNA_TYPE;
        }
        else {
            $type = 'UNKNOWN';
            warning('Unknow type for '.$key."\n");
        }
    }
    elsif ($key eq 'CDS') {
        #############
        # TRANSLATION

        throw('Translation table is '.$feature->{$key}->{transl_table}.' instead of '.$MIT_CODON_TABLE."\n") unless ($MIT_CODON_TABLE eq $feature->{$key}->{transl_table});
        my $translation = new  Bio::EnsEMBL::Translation(
                             -START_EXON => $start_exon,
                             -END_EXON   => $end_exon,
                             -SEQ_START  => 1,
                             -SEQ_END    => $transcript->length,
                             );

        my %h_dbentry;
        foreach my $dbentry (@{$dbe_adaptor->fetch_all_by_name($feature->{$key}->{gene})}) {
            next unless $dbentry->info_type eq 'UNMAPPED';
            push(@{$h_dbentry{$dbentry->primary_id}}, $dbentry);
        }
        if (!exists $opt{noxref}) {
            foreach my $db_xref_key (keys %{$feature->{$key}->{db_xref}}) {
                next if ($db_xref_key eq 'GI');
                if (exists $h_dbentry{$feature->{$key}->{db_xref}->{$db_xref_key}}) {
                    # add dbentry to transcript only if they are of unmapped type
                    warning('Xref already exists for '.$feature->{$key}->{db_xref}->{$db_xref_key}."\t".$feature->{$key}->{gene}."\n");
                    foreach my $dbentry (@{$h_dbentry{$feature->{$key}->{db_xref}->{$db_xref_key}}}) {
                        $dbentry->info_type(INFO_TYPE);
                        $dbentry->info_text(INFO_TEXT);
                        $gene->add_DBEntry($dbentry);
                    }
                }
                else {
                    if (exists $EXTERNAL_DB{$db_xref_key}) {
                        my $db_entry = Bio::EnsEMBL::DBEntry->new(
                          -adaptor     => $dbe_adaptor,
                          -primary_id  => $feature->{$key}->{db_xref}->{$db_xref_key},
                          -display_id  => $feature->{$key}->{gene},
                          -description => $feature->{$key}->{product},
                          -dbname      => $EXTERNAL_DB{$db_xref_key},
                          -info_type   => INFO_TYPE,
                          -info_text   => INFO_TEXT,
                        );
                        $gene->add_DBEntry($db_entry);
                    }
                    else {
                        warning('Not using xref '.$feature->{$key}->{db_xref}->{$db_xref_key}. ' from '.$db_xref_key);
                    }
                }
            }
            if (exists $feature->{$key}->{protein_id}) {
                my $dbname = $feature->{$key}->{protein_id} =~ /^NP/ ? REFSEQ_PEP : REFSEQ_XPEP;
                my $dbxentry = Bio::EnsEMBL::DBEntry->new(
                  -adaptor     => $dbe_adaptor,
                  -primary_id  => $feature->{$key}->{protein_id},
                  -display_id  => $feature->{$key}->{gene},
                  -description => $feature->{$key}->{product},
                  -dbname      => $dbname,
                  -info_type   => INFO_TYPE,
                  -info_text   => INFO_TEXT,
                );
                $translation->add_DBEntry($dbxentry);
            }
        }
        if (exists $feature->{$key}->{note} and $feature->{$key}->{note} =~ /frameshift/) {
            warning('There is a frameshift in: '.$feature->{$key}->{gene}."\nNote: ".$feature->{$key}->{note});
        }
        if (exists $feature->{$key}->{transl_except}) {
            if ($feature->{$key}->{transl_except} =~ /pos:(\d+),aa:(\w\+)/) {
                my $alt_seq;
                my $pos = $1;
                if ($2 eq 'TERM') {
                    $alt_seq = 'AA' if ($feature->{$key}->{transl_except} =~ /TAA/);
                    my $seq_edit = Bio::EnsEMBL::SeqEdit->new(
                          -CODE    => '_rna_edit', 
                          -NAME    => 'rna_edit', 
                          -DESC    => 'RNA edit', 
                          -START   => $pos,
                          -END     => $pos, 
                          -ALT_SEQ => $alt_seq
                          );
                    $translation->add_Attributes($seq_edit->get_Attribute());
                }
            }
        }
        if (exists $feature->{$key}->{fragment}) {
            warning('The gene '.$feature->{$key}->{gene}." is fragmented, a methionine will be added!\n");
        }
        $transcript->translation($translation);
        if ($transcript->translate()->seq() =~ /\*/) {
          throw("Stop codon found in translation ".$transcript->translate()->seq());
        }
        if ($transcript->translate()->seq() =~ /^[^M]/) {
          warning("Adding SeqEdit for non-methionine start codon in translation ".$transcript->translate()->seq());
          my $seqedit = Bio::EnsEMBL::SeqEdit->new(
                -CODE    => 'amino_acid_sub', 
                -NAME    => 'Amino acid substitution', 
                -DESC    => 'Some translations have been manually curated for amino acid substitiutions. For example a stop codon may be changed to an amino acid in order to prevent premature truncation, or one amino acid can be substituted for another.', 
                -START   => 1,
                -END     => 1, 
                -ALT_SEQ => 'M' 
                );
          $transcript->translation()->add_Attributes($seqedit->get_Attribute());
        }
        $type = $MIT_GENE_TYPE;
    }
    eval {
      $gene->biotype($type);
      $gene->analysis($ensembl_analysis);
      $gene->status(STATUS);
      $gene->description($feature->{$key}->{product});
      $transcript->biotype($type);
      $transcript->analysis($ensembl_analysis);
      $gene->add_Transcript($transcript);
      $count++;
    };
    if ($@){
      print "Error: $@\n";
      exit;
    } 
    printf "\t%-12s %5d\n\t%-12s %5d\n\t%-12s %s\n\t%-12s %s\n***********************************************\n\n", 'Start:',  $gene->start, 'End:', $gene->end, 'Description:', $gene->description, 'Biotype:', $gene->biotype;
    ++$entry_count;
    push @genes,$gene; 
  } 


print "Have ".scalar(@genes)." gene objects\n";   

print " LOADING : Do you want to load them into the db ? (Y/N) ";
  my $answer = <>;
  chomp $answer;
  if ($answer eq "y" or $answer eq "Y"){
    &load_db($output_db,\@genes);
  }
else{
  exit 0;
}


######
#TEST

exit 0 unless $MIT_DEBUG;

print "\n\n################################################################
# Testing gene load\n\n" if $MIT_DEBUG;
my $new_genes_adaptor = $output_db->get_GeneAdaptor;
my @new_genes = @{$new_genes_adaptor->fetch_all_by_Slice($slice)};
if (!scalar(@new_genes)) {
  throw("No genes loaded");
}  

foreach my $new_gene (@new_genes){
  print ref($new_gene)."\t with dbID ".$new_gene->dbID."\n" ;
  foreach my $new_transcript (@{$new_gene->get_all_Transcripts()}) {
    print "Transcript:\n".$new_transcript->seq->seq."\n"; 
    if ($new_transcript->translation) {
      print "Translation:\n".$new_transcript->translate()->seq()."\n" if $MIT_DEBUG;
    }
  }
}

exit 0;

sub load_db(){
  my ($output_db,$genes)=@_;
  print "storing genes\n" ; 
  #############
  # STORE GENES

  foreach my $gene(@{$genes}){
    print "Loading gene ",$gene,"\t" if $MIT_DEBUG;
    my $dbid = $output_db->get_GeneAdaptor->store($gene);
    print "dbID = $dbid\n" if $MIT_DEBUG;
    my $stored_gene = $output_db->get_GeneAdaptor->fetch_by_dbID($dbid);
  }
  return 0;
}

################################
# Get the sequence if requested

sub get_chromosomes {
  my ($genbank_ref,$output_db,) = @_;
  my %slices;

  my %assembly;
  $assembly{$MIT_TOPLEVEL} = $MIT_NAME;
  if ($MIT_SCAFFOLD_SEQNAME) {
      $assembly{$scaffold} = $MIT_SCAFFOLD_SEQNAME;
  }
  else {
      $assembly{$scaffold} = $genebank_hash->{-accession};
  }
  if ($MIT_CLONE_SEQNAME) {
      $assembly{'clone'} = $MIT_CLONE_SEQNAME;
  }
  else {
      if ( $genebank_hash->{_comment} =~ /reference\s+sequence\s+was\s+derived\s+from\s+(\w+)./ ) {
          $assembly{'clone'} = "$1";
      }
  }
  if ($MIT_CONTIG_SEQNAME) {
      $assembly{'contig'} = $MIT_CONTIG_SEQNAME;
  }
  else {
      if ( $genebank_hash->{_comment} =~ /reference\s+sequence\s+was\s+derived\s+from\s+(\w+)./ ) {
          $assembly{'contig'} = $1.'.1'.$genebank_hash->{_length};
      }
  }

  my $csa = $output_db->get_CoordSystemAdaptor();
  my $sa  = $output_db->get_SliceAdaptor();
  # Get all coord systems in the database:
  # Make a slice for each coord system

  foreach my $cs (@{$csa->fetch_all()}) {
    my $name = $cs->name;
    $name =  'top_level' if ($cs->name eq $MIT_TOPLEVEL);
    $name =  'seq_level' if ($cs->is_sequence_level);
    if ($assembly{$cs->name}){
      $slices{$name}  = Bio::EnsEMBL::Slice->new
      (
        -coord_system      => $cs,
        -start             => 1,
        -end               => $genbank_ref->{_length},
        -strand            => 1,
        -seq_region_name   => $assembly{$cs->name},
        -seq_region_length => $genbank_ref->{_length},
      );
   
    }
  }

  # Die before storing anything unless you have sequences that are top level and seq level
  # Unless you only have one coord system in which case you set it to both seq and top level
  die "Cannot find seq_level coord system" unless $slices{'seq_level'};
  die "Cannot find top_level coord system $MIT_TOPLEVEL"
    unless (    scalar( keys %slices ) > 1 && $slices{'top_level'}
             or scalar( keys %slices ) == 1 );

return \%slices;

}

sub load_chromosomes {
  my ( $slices, $output_db, $seq_ref ) = @_;
  my $sa = $output_db->get_SliceAdaptor();
  my $aa = $output_db->get_AttributeAdaptor();
  # Store slices, add the mit codon usage atribute
  # add the sequence if the coord system is contig
  # add the top level attribute if the coord system is chromosome
  # Make mitochondrial codon usage attribute

  push my @codonusage,
    Bio::EnsEMBL::Attribute->new( -CODE        => 'codon_table',
                                  -NAME        => 'Codon Table',
                                  -DESCRIPTION => 'Alternate codon table',
                                  -VALUE       => $MIT_CODON_TABLE, );
  # Make top level seq attribute
  push my @toplevel,
    Bio::EnsEMBL::Attribute->new(
                    -CODE        => 'toplevel',
                    -NAME        => 'Top Level',
                    -DESCRIPTION => 'Top Level Non-Redundant Sequence Region',
                    -VALUE       => 1 );

  foreach my $cs (  sort keys %slices ) {
    print "Slice " . $slices{$cs}->name . "\n" if $MIT_DEBUG;
    if ( $cs eq 'seq_level' ) {
      $sa->store( $slices{$cs}, \$seq_ref );
      $aa->store_on_Slice( $slices{$cs}, \@codonusage );
      $slices{$cs}->adaptor($sa);
      print "Storing seqlevel \n" if $MIT_DEBUG;
      # If only have 1 coord systen it needs to be both seq_level
      # and top level
      if ( scalar( keys %slices ) == 1 ) {
        $aa->store_on_Slice( $slices{$cs}, \@toplevel );
      }
      next;
    }
    print "Storing slice \n" if $MIT_DEBUG;
    $sa->store( $slices{$cs} );
    $aa->store_on_Slice( $slices{$cs}, \@codonusage );
    $slices{$cs}->adaptor($sa);
    if ( $cs eq 'top_level' ) {
      $aa->store_on_Slice( $slices{$cs}, \@toplevel );
    }
  }

  # if you only have 1 coordsystem dont need an assembly
  return 0 if ( scalar( keys %slices ) == 1 );

  # load the assembly
  # Load a chromosome - scaffold entry in the asembly, if these
  # coord stestems exist

  if ( $slices{'top_level'} && $slices{$scaffold} ) {
    print "Making assembly for chromosome vs scaffold\n" if $MIT_DEBUG;
    &load_assembly( $slices{'top_level'}->get_seq_region_id,
                    $slices{'top_level'}->start,
                    $slices{'top_level'}->end,
                    $slices{$scaffold}->get_seq_region_id,
                    $slices{$scaffold}->start,
                    $slices{$scaffold}->end,
                    1,
                    $output_db );
  }

  # Load assemby tables for each other coord system vs seq_level

  foreach my $cs ( keys %slices ) {
    print "Slice " . $slices{$cs}->name . "\n" if $MIT_DEBUG;
    next if ( $cs eq 'seq_level' );
    print "Making assembly for $cs vs seq level\n" if $MIT_DEBUG;
    &load_assembly( $slices{$cs}->get_seq_region_id,
                    $slices{$cs}->start,
                    $slices{$cs}->end,
                    $slices{'seq_level'}->get_seq_region_id,
                    $slices{'seq_level'}->start,
                    $slices{'seq_level'}->end,
                    1,
                    $output_db );
  }
  return 0;
} ## end sub load_chromosomes

##################################################################
# Do the sql statement to load the values into the assembly table

sub load_assembly {
  my ( $chr_id, $chr_start, $chr_end,
       $contig, $contig_start, $contig_end,
       $contig_ori, $db ) = @_;

  if ( !$contig ) {
    die "contig id must be defined for this to work\n";
  }
  my $sql = "insert into assembly(asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori) values(?, ?, ?, ?, ?, ?, ?)";
  my $sth = $db->dbc->prepare($sql);
  $sth->execute( $chr_id, $chr_start, $chr_end, $contig, $contig_start,
                 $contig_end, $contig_ori );
  return 0;
}

sub check_codon_table {
  my ($out_db,$mt_slice) = @_;
  my $found;
  my $mt_slice_dbID = $slice_adaptor->get_seq_region_id($mt_slice);

  my $sql = "SELECT sra.value ".
            "FROM seq_region_attrib sra, attrib_type att ".
            "WHERE sra.attrib_type_id = att.attrib_type_id ".
            "AND att.code = 'codon_table' ".
            "AND sra.seq_region_id = ". $mt_slice_dbID;
  my $sth = $out_db->dbc->prepare($sql) or die "sql error!";
  $sth->execute();
  my $val = $sth->fetchrow();
  $sth->finish;


  if ($val == CHECK_CODON_TABLE) {
    $found = 1;
  }
  return $found;
}
##################################################
# Example of genbank entry


#     CDS             join(13818..13986,16435..16470,18954..18991,20508..20984,
#                     21995..22246,23612..23746,25318..25342,26229..26701)
#                     /gene="COX1"
#                     /locus_tag="Q0045"
#                     /EC_number="1.9.3.1"
#                    /note="Subunit I of cytochrome c oxidase, which is the
#                     terminal member of the mitochondrial inner membrane
#                     electron transport chain; one of three
#                     mitochondrially-encoded subunits;
#                     go_component: respiratory chain complex IV (sensu
#                     Eukaryota) [goid GO:0005751] [evidence IPI] [pmid
#                     1331058];
#                     go_function: cytochrome-c oxidase activity [goid
#                     GO:0004129] [evidence IDA] [pmid 1331058];
#                     go_process: aerobic respiration [goid GO:0009060]
#                     [evidence IMP] [pmid 9724417]"
#                     /codon_start=1
#                     /evidence=experimental
#                     /transl_table=3
#                     /product="Cox1p"
#                     /protein_id="NP_009305.1"
#                     /db_xref="GI:6226519"
#                     /db_xref="SGD:S000007260"
#                     /db_xref="GeneID:854598"
#                     /translation="MVQRWLYSTNAKDIAVLYFMLAIFSGMAGTAMSLIIRLELAAPG
#                     SQYLHGNSQLFNVLVVGHAVLMIFFLVMPALIGGFGNYLLPLMIGATDTAFPRINNIA
#                     FWVLPMGLVCLVTSTLVESGAGTGWTVYPPLSSIQAHSGPSVDLAIFALHLTSISSLL
#                     GAINFIVTTLNMRTNGMTMHKLPLFVWSIFITAFLLLLSLPVLSAGITMLLLDRNFNT
#                     SFFEVSGGGDPILYEHLFWFFGHPEVYILIIPGFGIISHVVSTYSKKPVFGEISMVYA
#                     MASIGLLGFLVWSHHMYIVGLDADTRAYFTSATMIIAIPTGIKIFSWLATIHGGSIRL
#                     ATPMLYAIAFLFLFTMGGLTGVALANASLDVAFHDTYYVVGHFHYVLSMGAIFSLFAG
#                     YYYWSPQILGLNYNEKLAQIQFWLIFIGANVIFFPMHFLGINGMPRRIPDYPDAFAGW
#                     NYVASIGSFIATLSLFLFIYILYDQLVNGLNNKVNNKSVIYNKAPDFVESNTIFNLNT
#                     VKSSSIEFLLTSPPAVHSFNTPAVQS"
