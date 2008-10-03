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
use Carp;
#use Bio::EnsEMBL::Pipeline::Config::MitConf qw(%MitConf);
use MitConf;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis;
use Bio::SeqIO;

use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace throw);


my $help;
my @genes;;
my $translation; # JUST INCASE
my $count;
#use scaffold or supercontig:
my $supercontig = 'supercontig';
my %opt;


 # options submitted with commandline override MitConf.pm 

GetOptions(
           \%opt ,
           '-h|help'    , 
           'dbhost=s' , 
           'dbuser=s' , 
           'dbpass=s' , 
           'dbport=i' , 
           'dbname=s' , 
           'contig=s',      # MIT_CONTIG_SEQNAME 
           'chromosome=s',  # MIT_CHROMOSOME_SEQNAME 
           'supercontig=s', # MIT_SUPERCONTIG_SEQNAME
           'clone=s',       # MIT_CLONE_SEQNAME 
           'toplevel=s',    # MIT_TOPLEVEL 
           'gene_type=s', 
           'trna_type=s',
           'rrna_type=s', 
           'codon_table=i', 
           'name=s' , 
           'genbank_file=s' , 
           ) ; # or &usage();

if ($opt{dbhost} && $opt{dbuser} && $opt{dbname} && $opt{dbpass} && $opt{dbport} ) {  
  $MIT_DBHOST  = $opt{dbhost} ; 
  $MIT_DBUSER = $opt{dbuser} ;  
  $MIT_DBPASS = $opt{dbpass} ; 
  $MIT_DBPORT = $opt{dbport} ; 
  $MIT_DBNAME = $opt{dbname} ; 
}

$MIT_GENBANK_FILE = $opt{genbank_file} if $opt{genbank_file} ; 
$MIT_LOGIC_NAME =  $opt{logic_name} if $opt{logic_name} ; 
$MIT_NAME =  $opt{name} if $opt{name} ; 
$MIT_TOPLEVEL =  $opt{toplevel} if $opt{toplevel} ; 
$MIT_CODON_TABLE =  $opt{codon_table} if $opt{codon_table} ; 
$MIT_GENE_TYPE =  $opt{gene_type} if $opt{gene_type} ; 
$MIT_TRNA_TYPE =  $opt{trna_type} if $opt{trna_type} ; 
$MIT_RRNA_TYPE =  $opt{rrna_type} if $opt{rrna_type} ; 


unless ($MIT_DBHOST && $MIT_DBUSER && $MIT_DBNAME && $MIT_GENBANK_FILE && !$help){
  warn("Can't run without MitConf.pm values:
MIT_DBHOST $MIT_DBHOST 
MIT_DBUSER $MIT_DBUSER 
MIT_DBNAME $MIT_DBNAME
MIT_DBPASS $MIT_DBPASS
MIT_GENBANK_FILE $MIT_GENBANK_FILE
MIT_LOGIC_NAME $MIT_LOGIC_NAME
MIT_NAME $MIT_NAME
MIT_TOPLEVEL $MIT_TOPLEVEL
MIT_CODON_TABLE  $MIT_CODON_TABLE
MIT_GENE_TYPE $MIT_GENE_TYPE
MIT_TRNA_TYPE $MIT_TRNA_TYPE
MIT_RRNA_TYPE $MIT_RRNA_TYPE
");
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}
############################
#PARSE GENBANK FILE FIRST

my ($genbank_ref,$assembly_ref)  = &_parse_coordinates; 

# filter the genbank-entries in the file because the tRNA's are 
# redundant ( ( some tRNA features are annotated by tRNAscan-SE 
# which we don't want to load ) 

my ($s,$e) ;  
my %nonred_entries ; 

my %why1 = %{ ${$genbank_ref}[0] } ;
my %why2 = %{ ${$genbank_ref}[1] } ;


my $save_entry ; 
GENBANK_ENTRY: for my $g (@$genbank_ref) {   

   my %entry = %$g ;  

   # 
   #  filter / find entries :
   #  kill all entries which have at leat 2 times the same start/end and an additional note 
   #  "anotated by tRNA-scan-SE" 
    
   #  get start/end and note out of the entry   
   my $start="" ; 
   my $end = "" ; 
   my $note ="";

   foreach my $k (sort keys %entry) {  

     if ($k=~/start/) {   
       my @starts = @{$entry{$k}} ;    
       throw ( " not a single-exon-gene ") if scalar(@starts) > 1  ;  
       $start = join ("-", @starts)  ;   

     } elsif ($k=~/end/) {  
       my @starts = @{$entry{$k}} ;   
       throw ( " not a single-exon-gene ") if scalar(@starts) > 1  ;   
       $end = join ("-", @starts)  ;  

     } elsif ($k=~m/note/){  
       $note = $entry{$k} ;  
     } elsif ($k=~m/organism/){      
        $save_entry = $g ;
        next GENBANK_ENTRY  ; 
     } 
   }   
   
     if (!exists $nonred_entries{"$start-$end"} && $start && $end ){
       if ( $note!~m/tRNAscan-SE/) {  
         $nonred_entries{"$start-$end"}=\%entry;   
       }else {  
         print " feature with coords $start-$end will be skipped because $note\n" ; 
       } 
     }
}

#my @genbank = @{$genbank_ref};  
my @genbank ; 
push @genbank , $save_entry ; 
push @genbank , values %nonred_entries ;
 
my %assembly = %{$assembly_ref};

##########################
# Get chromosome sequence 

my $seq_file = Bio::SeqIO->new(
			       -file   => $MIT_GENBANK_FILE,
			       -format => 'genbank',
			      );
my $chromosome_seq = $seq_file->next_seq;

####################
#Open db for writing

my $output_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $MIT_DBHOST,
						   '-user'   => $MIT_DBUSER,
						   '-pass'   => $MIT_DBPASS,
						   '-dbname' => $MIT_DBNAME,
						   '-port'   => $MIT_DBPORT,	
						  );

########################
# Check taxon is correct

my $meta_container = $output_db->get_MetaContainer();
my $db_taxon_id = $meta_container->get_taxonomy_id(); 
my $gb_taxon_id ; 

my %genome ; 

foreach my $entry (@genbank ){ 
  if ($$entry{source}){  
     $gb_taxon_id = $$entry{'db_xref'};   
     #  print "\n\n\n$gb_taxon_id \n\n\n" ; 
      %genome = %{$entry} ; 
  }
}

#$gb_taxon_id =~ s/^taxon://;
#if ($db_taxon_id  == $gb_taxon_id){
#  print "Genbank taxon id $gb_taxon_id matches the taxon id in the meta table\n";
#}else {
#  croak "taxonomy id in the genbank file $gb_taxon_id is not equal to the id found in the meta table $db_taxon_id \n";
#}
#
###########################################
# write chromosome if it is not already there

my $dbe_adaptor = $output_db->get_DBEntryAdaptor;
my $slice_adaptor = $output_db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('toplevel',$MIT_NAME);
my $slices_ref =  &get_chromosomes($assembly_ref,$genbank_ref, $output_db);
my %slices = %{$slices_ref};

if ($slice){
  print "Found chromosome ".$slice->name."\n" ;
}  else {
  print "There is no chromosome in $MIT_DBNAME called $MIT_NAME\nHave found the following slices to load:\n";
  foreach my $cs (keys %slices){
    print "$cs    \t".$slices{$cs}->name."\n";
  }
  print "Do you want to load the chromosome and the associated assembly entries
from the genbank file?(Y/N) "; 
  my $answer = <>;
  chomp $answer;
  if ($answer eq "y" or $answer eq "Y"){
    &load_chromosomes(\%slices, $output_db, $chromosome_seq );
    $slice = $slice_adaptor->fetch_by_region('toplevel',$MIT_NAME);
  }
  else {
    print "Ok gene load aborted.\n";
    exit 0;
  }
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
    print "You choose not to remove the genes\n";
    print "Do you want to keep loading the MT genes?(This may create duplicated entries)(Y/N)?\n";
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

my $logic_name = 'ncRNA';
my $ncRNA_analysis = $output_db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
if(!defined $ncRNA_analysis){
#croak "analysis $logic_name not found\n";
  $ncRNA_analysis = Bio::EnsEMBL::Analysis->new(
                                                -logic_name => 'ncRNA',
                                               );
  print "You have no ".$logic_name." defined creating new object\n";
}

if ($MIT_LOGIC_NAME){
  $logic_name = $MIT_LOGIC_NAME;
} else {
  $logic_name = 'ensembl';
  print "Cannot find MIT_LOGIC_NAME - using standard logic name ensembl\n" if $MIT_DEBUG;
}
my $ensembl_analysis = $output_db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
if(!defined $ensembl_analysis){
#croak "analysis $logic_name not found\n";
  $ensembl_analysis = Bio::EnsEMBL::Analysis->new(
                                                  -logic_name => 'ensembl',
                                                 );
  print "You have no ".$logic_name." defined creating new object\n";
}


 #########################################################
 # Create gene objects from array of hashes
 # contining parsed coords
 # index 0 contains the accession 
 # index 1 contains source data (taxon and sequence file)
 # subsequent indecies contain all the other annotations 
 # 

print " have " . scalar(@genbank) . " ENTRIES \n" ;  

for (my $i=1; $i <= $#genbank; $i++){


   
   my %h = %{$genbank[$i]} ;  
   print "ENTRY $i\n";  
   foreach ( keys %h) { 
     print "$_\t$h{$_}\n"  ;   
     if ( /start/ || /end/) { 
      my @tmp = @{$h{$_}}; 
       print "STARTS :\n" ;  
      print join ("\t" , @tmp ) . "\n" ;  
     }
   } 
   print "\nXXXXXXXXXXXXXXXXXXXXXX\n" ; 


  my $desc = $genbank[$i]{'product'};
  my $type = $genbank[$i]{'type'};
  my $status = "KNOWN";

  my $acc;
  if ($genbank[$i]{'protein_id'}){
    $acc =  $genbank[$i]{'protein_id'};
  }
  else{
    $acc = $why1{'accession'}."  ".${$genbank[$i]{'start'}}[0]."  ".${$genbank[$i]{'end'}}[0]; 
  }

  my $strand =  $genbank[$i]{'strand'};
  if ($type eq 'tRNA'){
    $type = $MIT_TRNA_TYPE;
  }
  if ($type eq 'rRNA'){
    $type = $MIT_RRNA_TYPE;
  }
  print "\n***********************************************************\n" ; 
  # By default analysis is set to ncRNA
  my $analysis = $ncRNA_analysis;
  my $transcript = new Bio::EnsEMBL::Transcript;
  my $start_exon;
  my $end_exon;
  my $total_exons = scalar(@{$genbank[$i]{'start'}}-1);
  # exons are stored in an array within the hash incase of multi exon genes

  for (my $exon_number=0; $exon_number <= $total_exons ; $exon_number ++){
    my $start = $genbank[$i]{'start'}[$exon_number];
    my $end   = $genbank[$i]{ 'end' }[$exon_number];

    ###########################################################
    # Dodgy hack to solve problem of frameshifts in pseudogenes

    if ($exon_number > 0 && $start <= $genbank[$i]{ 'end' }[$exon_number-1]){
      $start = $genbank[$i]{ 'end' }[$exon_number-1]+1;
    }

    #######
    # EXONS

    my $exon = new Bio::EnsEMBL::Exon;
    $exon->start($start); 
    $exon->end($end);
    $exon->strand($strand);
    $exon->slice($slice);
    $exon->phase(0);
    $exon->end_phase(($exon->end - $exon->start + 1)%3);

    #############
    # TRANSCRIPTS

    $transcript->add_Exon($exon);
    $start_exon = $exon if($exon_number == 0);
    $end_exon = $exon if($exon_number == $total_exons);
  } 

  $transcript->start_Exon($start_exon);
  $transcript->end_Exon  ( $end_exon );

  #############
  # TRANSLATION

  if ($type eq 'CDS'){
    $translation = new  Bio::EnsEMBL::Translation(
						 -START_EXON => $start_exon,
						 -END_EXON   => $end_exon,
						 -SEQ_START  => 1,
						 -SEQ_END    => $transcript->length,
						 );
    $transcript->translation($translation);
    $analysis = $ensembl_analysis;
    $type = $MIT_GENE_TYPE;
   }
  $transcript->biotype($type);
  #############
  # Make  genes
   print "create new trans " . $transcript->seq_region_start."\n" ; 

  my $gene = new Bio::EnsEMBL::Gene;
  eval {
    $gene->biotype($type);
    $gene->analysis($analysis);
    $gene->status($status);
    $gene->description($desc);
    $transcript->biotype($type);
    $transcript->analysis($analysis);
    $gene->add_Transcript($transcript);
    $count++;
  };
  if ($@){
    print "Error: $@\n";
    exit;
  } 
  print " pushing " . $gene->seq_region_start . "\n" ; 
  push @genes,$gene; 
}


print "Have ".scalar(@genes)." gene objects\n";   

@genes = sort {$a->seq_region_start <=> $b->seq_region_start } @genes ;  

foreach my $gene (@genes ) { 
  print  $gene->seq_region_start . "\n" ;  
}  

print " TESTING : Do you want to load them into the db ? (Y/N) ";
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
my @new_genes;
eval{
  @new_genes = @{$new_genes_adaptor->fetch_all_by_Slice($slice)};
};
  if($@){
    print "error in fetching from slice ".$slice->name."\n$@\n";
  }

foreach my $new_gene (@new_genes){
  print $new_gene."\t".$new_gene->dbID."\n" ;
  my $new_transcript = @{$new_gene->get_all_Transcripts()}[0];
  eval{
    if ($new_transcript->translate()->seq()){
      my $codon = $new_transcript->seq->seq;
      print "Transcript:\n$codon\n" if $MIT_DEBUG;
      print "Translation:\n".$new_transcript->translate()->seq()."\n" if $MIT_DEBUG;
    }
  };
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
  }
  return 0;
}



#######################################
# PARSE COORDINATES OUT OF GENBANK FILE

sub _parse_coordinates(){

  print "Parsing coords\n" if $MIT_DEBUG;

  my $desc;
  my $acc;
  my $start;
  my $end;
  my $strand;
  my $type;
  my $go = 'stop';
  my $index=-1;
  my @entries;
  my %entry;
  my %assembly;
  my $comment;
  my @genbank_file;
  $assembly{$MIT_TOPLEVEL} = $MIT_NAME;

#################################################
# Read genbank file into array, join broken lines


  open (GENBANK,$MIT_GENBANK_FILE) or die "Cannot open genbank file $MIT_GENBANK_FILE\n";
  my $line;
  while (<GENBANK>){
    chomp;
    $_ =~ s/^\s{5}//;
    $line .= $_;
    next if ($line =~ /\,$/);
    push @genbank_file,$line;
    $line = "";
  }
  close GENBANK;

####################################
# Array of hashes holds genbank data

my $first_entry = 0  ; 

  # read all entries and than filter the predictions 

  for (@genbank_file){
    # ignore sequence data
    next if $_ =~ /^\d+/;  

    my @sources = split/\s+/ ; 

    # make comma delimited
    $_ =~ s/\s+/,/g;
    # comment line holds contig name
    if ($_ =~ /^COMMENT/){
      $go = "comment";
    }
    # ignore features section
    if ($_ =~ /^FEATURES/){
      $go = "stop";
    }
    if ($go eq "comment"){
      $comment .= $_;
      next;
    }
    # stop parsing unless line starts wih a comma:
    unless ($_ =~ /^\,/ ){
      $go = "stop";
    }
    # Types of entries to parse all others are ignored, could be extended if needed
    if ($_ =~ /^tRNA/ || $_ =~ /^rRNA/ || $_ =~ /^CDS/ || $_ =~ /^ACCESSION/ || $_ =~ /^source/ ){ 
      $index++;
      push(@entries,%entry);
      # splits the line into words and parses them into the hash
      $_ =~ s/^\,//;
      $_ =~ s/\,/\./g;
      $_ =~ s/\)//g;
      my @string = split(/\./,$_); 

        if ($_ =~ /ACCESSION/){
	  $entries[$index]{'accession'}=$string[1]; 
	} 

	else{
	  if  ($string[1]=~ 'complement'){
	    $entries[$index]{'strand'}='-1';
	  }
	  else{
	    $entries[$index]{'strand'}='1';
	  }
	  $string[0] =~ s/\s+//g;
	  $string[1] =~ s/\D+//g;

	  $entries[$index]{'type'}=$string[0];

	  ##################################################
	  # Pushes starts and stops into array to accomodate
	  # The rare occurances of multiexon genes
	  # assumes coords alternate between start and stop

	  my $key = 'start';
	  for(my $i =1 ; $i <= $#string ; $i++){
	    if ($string[$i] =~ /\d+/){
	      push @{$entries[$index]{$key}} , $string[$i];
	      if ($key eq 'end'){ 
		  $key = 'start';
		  next;
		}
	      if ($key eq 'start'){
		$key = 'end'; 
		next;
	      }
	    }
	  }
	}
      # got the first line of the entry, now parse subsequent lines
      $go = 'go';
    }

    if ($go eq 'go'){
      my @temp = split(/\"/,$_);
      if ($temp[1]){
	$temp[0] =~ s/\W+//g;
	$temp[1] =~ s/,+/ /g;
	$entries[$index]{$temp[0]}=$temp[1];
      }
      else {
	if($temp[0] && $temp[0] eq ',/pseudo'){
	  $entries[$index]{'type'} = 'pseudogene';
	}
      }
    }
  }

  for (my $array_index =0; $array_index <= $#entries ;$array_index ++){
    # Check for trans-splicing events
    if ( $entries[$array_index]{'note'} && $entries[$array_index]{'note'} eq 'trans-splicing'){
      print "Genbank file contains trans-splicing event that this script cannot parse\nexiting...\n";
      exit 0;
    }
  }
  
  #############################################################
  # Use config file vales if they are present, otherwise get
  # seq names from GFF file

  if ($MIT_SUPERCONTIG_SEQNAME){
    $assembly{$supercontig} = $MIT_SUPERCONTIG_SEQNAME;
  }
  else {
    $assembly{$supercontig} = $entries[0]{'accession'};
  }  
  if ($MIT_CLONE_SEQNAME){
    $assembly{'clone'} = $MIT_CLONE_SEQNAME;
  }
  else {
    if ($comment =~ /reference,sequence,was,derived,from,(\w+)./){
      $assembly{'clone'} = "$1";
    }
  }  
  if ($MIT_CONTIG_SEQNAME){
    $assembly{'contig'} = $MIT_CONTIG_SEQNAME
  }
  else { 
    if ($comment =~ /reference,sequence,was,derived,from,(\w+)./){
      $assembly{'contig'} = "$1.".@{$entries[1]{'start'}}[0].".".@{$entries[1]{'end'}}[0];  
    }
  }
  return \@entries,\%assembly;
}

################################
# Get the sequence if requested

sub get_chromosomes{
  my ($assembly_ref,$genbank_ref,$output_db,) = @_;
  my %assembly = %{$assembly_ref};
  my @genbank =@{$genbank_ref};
  my %slices;
  my $csa = $output_db->get_CoordSystemAdaptor();
  my $sa  = $output_db->get_SliceAdaptor();
  # Get all coord systems in the database:
  # Make a slice for each coord system

  foreach my $cs (@{$csa->fetch_all()}) {
    my $name = $cs->name;
    $name =  'top_level' if ($cs->name eq $MIT_TOPLEVEL);
    $name =  'seq_level' if ($cs->is_sequence_level);
    if($assembly{$cs->name}){
      $slices{$name}  = Bio::EnsEMBL::Slice->new
	(
	 -coord_system      => $cs,
	 -start             => $genbank[1]{'start'}[0],
	 -end               => $genbank[1]{'end'}[0],
	 -strand            => 1,
	 -seq_region_name   => $assembly{$cs->name},
	 -seq_region_length => $genbank[1]{'end'}[0]- $genbank[1]{'start'}[0]+1,
	 -adaptor           => $sa
	)
      }
  }

  # Die before storing anything unless you have sequences that are top level and seq level
  # Unless you only have one coord system in which case you set it to both seq and top level
  die "Cannot find seq_level coord system" unless $slices{'seq_level'};
  die "Cannot find top_level coord system $MIT_TOPLEVEL" unless 
    (scalar(keys %slices) > 1  && $slices{'top_level'} or scalar(keys %slices) == 1);
  
return \%slices;

}

sub load_chromosomes{
  my ($slices,$output_db,$seq_ref)=@_;
  my $sa  = $output_db->get_SliceAdaptor();
  my $aa  = $output_db->get_AttributeAdaptor();
  # Store slices, add the mit codon usage atribute 
  # add the sequence if the coord system is contig
  # add the top level attribute if the coord system is chromosome
  # Make mitochondrial codon usage attribute

  push my @codonusage , Bio::EnsEMBL::Attribute->new
    (-CODE => 'codon_table',
     -NAME => 'Codon Table',
     -DESCRIPTION => 'Alternate codon table',
     -VALUE       => $MIT_CODON_TABLE,
    );  
  # Make top level seq attribute
  push my @toplevel , Bio::EnsEMBL::Attribute->new
    (-CODE => 'toplevel',
     -NAME => 'Top Level',
     -DESCRIPTION => 'Top Level Non-Redundant Sequence Region',
     -VALUE => 1
    );

  foreach my $cs (keys %slices){
    print "Slice ".$slices{$cs}->name."\n"if $MIT_DEBUG;
    if ($cs eq 'seq_level'){
      $sa->store($slices{$cs},\$seq_ref->seq);
      $aa->store_on_Slice($slices{$cs}, \@codonusage);    
      $slices{$cs}->adaptor($sa);
      print "Storing seqlevel \n"if $MIT_DEBUG;
      # If only have 1 coord systen it needs to be both seq_level
      # and top level
     if (scalar(keys %slices) == 1){
	$aa->store_on_Slice($slices{$cs}, \@toplevel);
      }
      next ;
    }
    print "Storing slice \n"if $MIT_DEBUG;
    $sa->store($slices{$cs});
    $aa->store_on_Slice($slices{$cs}, \@codonusage);
    $slices{$cs}->adaptor($sa);
    if ($cs eq 'top_level'){
      $aa->store_on_Slice($slices{$cs}, \@toplevel);
    }
  }

  # if you only have 1 coordsystem dont need an assembly 
  return 0 if (scalar(keys %slices) == 1);

  # load the assembly
  # Load a chromosome - supercontig entry in the asembly, if these
  # coord stestems exist

  if ($slices{'top_level'} && $slices{$supercontig}){
    print "Making assembly for chromosome vs supercontig\n" if $MIT_DEBUG;
    &load_assembly
      (
       $slices{'top_level'}->get_seq_region_id,
       $slices{'top_level'}->start,
       $slices{'top_level'}->end,
       $slices{$supercontig}->get_seq_region_id,
       $slices{$supercontig}->start,
       $slices{$supercontig}->end,
       1,
       $output_db
      );
  }

  # Load assemby tables for each other coord system vs seq_level

  foreach my $cs (keys %slices){
    print "Slice ".$slices{$cs}->name."\n"if $MIT_DEBUG;
    next if ($cs eq 'seq_level');
    print "Making assembly for $cs vs seq level\n"if $MIT_DEBUG;
    &load_assembly 
      (
       $slices{$cs}->get_seq_region_id,
       $slices{$cs}->start,
       $slices{$cs}->end,
       $slices{'seq_level'}->get_seq_region_id,
       $slices{'seq_level'}->start,
       $slices{'seq_level'}->end,
       1,
       $output_db
      )
    }
  return 0;
}

##################################################################
# Do the sql statement to load the values into the assembly table

sub load_assembly{
  my ($chr_id, $chr_start, $chr_end, $contig, $contig_start, $contig_end, $contig_ori, $db) = @_;

  if(!$contig){
    #print STDERR "trying to insert into ".$chr_id." ".$chr_start." ".$chr_end."\n";
    die "contig id must be defined for this to work\n";
  }
  my $sql = "insert into assembly(asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori) values(?, ?, ?, ?, ?, ?, ?)";
  # print "$sql\n";
  my $sth = $db->dbc->prepare($sql);
  $sth->execute($chr_id, $chr_start, $chr_end, $contig, $contig_start, $contig_end, $contig_ori); 
return 0;
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
