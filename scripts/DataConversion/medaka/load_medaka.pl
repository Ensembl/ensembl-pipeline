#!/usr/local/ensembl/bin/perl -w

=head1 NAME

load_medaka.pl 

=head1 DESCRIPTION

Parses predicted medaka genes out of the given GFF file and writes them to the $output_db specified.
- Fill in the necessary variables in the config file, MedakaConf.pm. 
- It should not be necessary to change anything in load_medaka.pl unless the format of the gff file has changed, for example)


=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis;
use Bio::SeqIO;
use MedakaConf;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning);
my $help;
my %opt; 

# options submitted with commandline override MedakaConf.pm 
GetOptions(
           \%opt ,
           '-h|help'    , 
           'dbhost=s' , 
           'dbuser=s' , 
           'dbpass=s' , 
           'dbport=i' , 
           'dbname=s' ,            
           'gene_type=s', 
           'logic_name=s' , 
           'gff_file=s' ,
	   'coord_system=s' ,
	   'coord_system_version=s' , 
           ) ; 

if ($opt{dbhost} && $opt{dbuser} && $opt{dbname} && $opt{dbpass} && $opt{dbport} ) {  
  $MED_DBHOST = $opt{dbhost} ; 
  $MED_DBUSER = $opt{dbuser} ;  
  $MED_DBPASS = $opt{dbpass} ; 
  $MED_DBPORT = $opt{dbport} ; 
  $MED_DBNAME = $opt{dbname} ; 
}

$MED_GFF_FILE     = $opt{gff_file} if $opt{gff_file} ; 
$MED_LOGIC_NAME   = $opt{logic_name} if $opt{logic_name} ; 
$MED_GENE_TYPE    = $opt{gene_type} if $opt{gene_type} ; 
$MED_COORD_SYSTEM = $opt{coord_system} if $opt{coord_system} ; 
$MED_COORD_SYSTEM_VERSION = $opt{coord_system_version} if $opt{coord_system_version} ;

unless ($MED_DBHOST && $MED_DBUSER && $MED_DBNAME && $MED_GFF_FILE && !$help){
  warn("Can't run without MedakaConf.pm values:
MED_DBHOST $MED_DBHOST 
MED_DBUSER $MED_DBUSER 
MED_DBNAME $MED_DBNAME
MED_DBPASS $MED_DBPASS
MED_GFF_FILE $MED_GFF_FILE
MED_LOGIC_NAME $MED_LOGIC_NAME
MED_GENE_TYPE $MED_GENE_TYPE
MED_COORD_SYSTEM $MED_COORD_SYSTEM
MED_COORD_SYSTEM_VERSION $MED_COORD_SYSTEM_VERSION
MED_DNA_DBNAME $MED_DNA_DBNAME
MED_DNA_DBHOST $MED_DNA_DBHOST  
MED_DNA_DBUSER $MED_DNA_DBUSER
MED_DNA_DBPASS $MED_DNA_DBPASS
MED_DNA_DBPORT $MED_DNA_DBPORT
");
  $help = 1;
}

if ($help) {
    exec('perldoc', $0);
}

# open db for writing genes
my $output_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $MED_DBHOST,
						   '-user'   => $MED_DBUSER,
						   '-pass'   => $MED_DBPASS,
						   '-dbname' => $MED_DBNAME,
						   '-port'   => $MED_DBPORT,	
						  );


my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $MED_DNA_DBHOST,
						'-user'   => $MED_DNA_DBUSER,
						'-pass'   => $MED_DNA_DBPASS,
						'-dbname' => $MED_DNA_DBNAME,
						'-port'   => $MED_DNA_DBPORT,	
					       );


# end of set-up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# get all the scaffolds in the dna_db
print STDERR "Fetching all $MED_COORD_SYSTEM...\n";
my $dna_sa = $dna_db->get_SliceAdaptor();
my @scaffolds = @{$dna_sa->fetch_all($MED_COORD_SYSTEM,$MED_COORD_SYSTEM_VERSION)};
foreach my $scaff (@scaffolds){
  print STDERR "name = ".$scaff->name."\n" if $MED_DEBUG;
}



# Parse GFF file
print STDERR "Parsing $MED_GFF_FILE...\n";
parse_gff(\@scaffolds, $output_db, $dna_sa); 



## check that it's all working so far
#foreach my $g (sort keys %gene){
#  foreach my $exon (sort {$a <=> $b} keys %{$gene{$g}}){
#    print $g."\t".$exon."\t@{$gene{$g}{$exon}}\n" if $MED_DEBUG;
#  }
#}
#
#
#
## fix the frame / phase
#my $fixed_phases = fix_frames(\%gene);
#
#
#
#
#
## make gene objects
#my @gene_objs = @{make_gene_objs($fixed_phases, \@scaffolds, $output_db, $dna_sa)};
#print "Have ".scalar(@gene_objs)." gene objects.\n";   
#@gene_objs = sort {$a->seq_region_start <=> $b->seq_region_start } @gene_objs ;  
#foreach my $g (@gene_objs){
#  print "got gene_stable_id ".$g->stable_id."\n" if $MED_DEBUG;
#}
#
#
#
## load genes into db
#&load_genes($output_db, \@gene_objs);
print STDERR "Done!\n";



# start of subroutines
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Parse GFF file. Use a hash of hashes.
# OuterKey = gene_stable_id
# InnerKey = exon_number
sub parse_gff {
  my ($scaffolds, $output_db, $dna_sa) = @_;
  # example of GFF file.
  #	scaffold1       UTOLAPRE05100100001     initial-exon    129489  129606  .       +       0
  #	scaffold1       UTOLAPRE05100100001     internal-exon   129920  130027  .       +       2
  #	scaffold1       UTOLAPRE05100100001     internal-exon   130753  130839  .       +       2
  #	scaffold1       UTOLAPRE05100100001     final-exon      131859  132262  .       +       2
  #	scaffold6469    UTOLAPRE05100120178     single-exon     1604    2746    .       -       0

  # read in the file. 
  my %gff;
  open (GFF,$MED_GFF_FILE) or die "Cannot open gff file $MED_GFF_FILE\n";
  my $line;
  my $count;
  while (<GFF>){
    chomp;
    $line = $_;
    print STDERR "$line\n" if $MED_DEBUG;
    next if ($line =~ m/^\#/);
    my @fields = split/\s+/, $line;
    for (my $i=0; $i<scalar(@fields); $i++) {
      print STDERR $i."\t".$fields[$i]."\n" if $MED_DEBUG;
    }
    if ($fields[2] eq 'initial-exon' || $fields[2] eq 'single-exon') {
      print STDERR "Resetting count to 0, ".$fields[2]."\n" if $MED_DEBUG;
      $count = 0;      
    } else {
      $count++;
    }     
    #gene_id -> exon_number = (start, end, strand, frame, scaffold, feature)
    @{$gff{$fields[1]}{$count}} = ($fields[3], $fields[4], $fields[6], $fields[7], $fields[0], $fields[2]);     
    if ($fields[2] eq 'final-exon' || $fields[2] eq 'single-exon') {
      my $fixed_frames = fix_frames(\%gff);
      my @gene_obj = @{make_gene_objs($fixed_frames, $scaffolds, $output_db, $dna_sa)};
      load_genes($output_db, \@gene_obj);
      %gff = ();
    }
    $line = '';
  }    
  close GFF; 
}

sub load_genes(){
  my ($output_db,$genes) = @_;
#  print "Storing genes...\n" ; 
  foreach my $gene (@{$genes}){
    print STDERR "Loading gene ".$gene->stable_id."\t" ;#if $MED_DEBUG;
    my $dbid = $output_db->get_GeneAdaptor->store($gene);
    print "dbID = $dbid\n" if $MED_DEBUG;
  }
  return 0;
}

sub fix_frames {
  my ($gene_ref) = @_; 
  my %gene = %{$gene_ref};
  foreach my $g (keys %gene){
    foreach my $e (keys %{$gene{$g}}){
      if ($gene{$g}{$e}[3] == 3) {
        $gene{$g}{$e}[3] = 0;
      } elsif ($gene{$g}{$e}[3] == 1) {
        $gene{$g}{$e}[3] = 2;
      } elsif ($gene{$g}{$e}[3] == 2) {
        $gene{$g}{$e}[3] = 1;
      }
    }
  }
  return \%gene;
}


# create gene objects from hash %gene
sub make_gene_objs {
#  print STDERR "Creating Gene objects...\n";
  my ($gene_ref, $scaffolds, $output_db, $dna_sa) = @_;
  my %gene_hash = %{$gene_ref};  
  my $analysis = $output_db->get_AnalysisAdaptor->fetch_by_logic_name($MED_LOGIC_NAME);
  if(!defined $analysis){
    $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => 'gff_prediction');
    warning("You have no ".$MED_LOGIC_NAME." defined; creating new object.\n");
  }
  my $logic_name;
  my $gene;
  my $transcript;
  my $exon;
  my $translation;
  my $start_exon;
  my $end_exon;
  my $total_exons;
  my $scaffold_name;
  my $slice;
  my @genes;
      
  # loop through hash and make objects
  GENE: foreach my $gff_gene (keys %gene_hash) {
    $transcript = new Bio::EnsEMBL::Transcript;
    $total_exons = scalar(keys %{$gene_hash{$gff_gene}});
    EXON: foreach my $gff_exon (sort {$a <=> $b} keys %{$gene_hash{$gff_gene}}) {
      # make exon objects
      # gene_id -> exon_number = (start, end, strand, frame, scaffold, feature)
      $exon = new Bio::EnsEMBL::Exon;
      $exon->start(@{$gene_hash{$gff_gene}{$gff_exon}}[0]); 
      $exon->end(@{$gene_hash{$gff_gene}{$gff_exon}}[1]);
      if(@{$gene_hash{$gff_gene}{$gff_exon}}[2] eq '+'){
	$exon->strand(1);
      }else{
	$exon->strand(-1);
      }
      $exon->phase(@{$gene_hash{$gff_gene}{$gff_exon}}[3]);
      $exon->end_phase(($exon->end - $exon->start + 1)%3);
      $exon->analysis($analysis);
      foreach my $scaffold (@{$scaffolds}){
        if ($scaffold->name =~ m/^($MED_COORD_SYSTEM\:$MED_COORD_SYSTEM_VERSION\:$MED_PREFIX@{$gene_hash{$gff_gene}{$gff_exon}}[4]\:1\:.*)/){
	  $scaffold_name = $1;
	  #name = scaffold:MEDAKA1:HdrR_200510_scaffold980:1:63005:1 
	  print STDERR "Match! ".$scaffold->name." with ".$scaffold_name."\n" if $MED_DEBUG;
	}
      }
      $exon->slice($dna_sa->fetch_by_name($scaffold_name));   
      $transcript->add_Exon($exon); 
      if (@{$gene_hash{$gff_gene}{$gff_exon}}[5] eq 'initial-exon' || @{$gene_hash{$gff_gene}{$gff_exon}}[5] eq 'single-exon') {
        $start_exon = $exon;
      }
      if (@{$gene_hash{$gff_gene}{$gff_exon}}[5] eq 'final-exon' || @{$gene_hash{$gff_gene}{$gff_exon}}[5] eq 'single-exon') { 
        $end_exon = $exon;
      }
    } #EXON
    $transcript->start_Exon($start_exon);
    $transcript->end_Exon($end_exon);
    $transcript->analysis($analysis);    
    $translation = new  Bio::EnsEMBL::Translation(
						 -START_EXON => $start_exon,
						 -END_EXON   => $end_exon,
						 -SEQ_START  => 1,
						 -SEQ_END    => $transcript->length,
						 );  
    $transcript->translation($translation);   
    $transcript->biotype($MED_GENE_TYPE);
    # create genes
    $gene = new Bio::EnsEMBL::Gene;
    eval {
      $gene->biotype($MED_GENE_TYPE);
      $gene->analysis($analysis);
      $gene->stable_id($gff_gene);
      $gene->add_Transcript($transcript);
    };  
    if ($@){
      print "Error: $@\n";
      exit;
    }
    print STDERR " Made Gene obj ".$gene->stable_id." @ ".$gene->seq_region_start."\n" if $MED_DEBUG; 
    push @genes,$gene;     						   
  } #GENE
  return \@genes;
}                   
