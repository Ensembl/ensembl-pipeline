#!/usr/local/ensembl/bin/perl -w

# dump_seq_into_fastA.pl
# it reads a bit of sequence and dump it into a fasA file, to eb able to view it in Apollo

use strict;
use diagnostics;

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

# get a contig with a piece-of/entire  chromosome

my $input_id;
my $hmm = "/tmp/HMMs.pfam";

### rat ###
#my $dbhost    = 'ecs2f';
#my $dbuser    = 'ensro';
#my $dbname    = 'genewisedb_rat';
#my $dnadbname = 'rat_Nov02';
#my $dnadbhost = 'ecs2a';

### mouse ###
my $dbhost    = 'ecs2f';
my $dbuser    = 'ensro';
my $dbname    = 'genewisedb_mouse';
my $dnadbname = 'mus_musculus_core_11_3';
my $dnadbhost = 'ecs2d';


&GetOptions(
	    'input_id:s'     => \$input_id,
	    'hmm:s'          => \$hmm,
	   );

unless( $input_id && $hmm ){
  print STDERR "Usage: $0 -input_id -hmm\n";
  exit(0);
}

# connect to the database
my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       -host  => $dnadbhost,
					       -user  => $dbuser,
					       -dbname=> $dnadbname
					      );

my $db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
					     -host  => $dbhost,
					     -user  => 'ensadmin',
					     -pass  => 'ensembl',
					     -dbname=> $dbname,
					     -dnadb => $dnadb,
					    );

my ($chr,$start,$end);
if ($input_id =~ /(\S+)\.(\d+)-(\d+)/){
  $chr   = $1;
  $start = $2;
  $end   = $3;
}

my $slice   = 
  $db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end)->get_repeatmasked_seq(['RepeatMask']);;
my $genomic = "/tmp/genomic.$$";

open( GENOMIC, ">$genomic") ||  warn("Could not open $genomic $!");

my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			     '-fh'     => \*GENOMIC);

$seqout->write_seq($slice);

my $command = "genewisedb -pfam $hmm -dnas $genomic -alg 623 -aln 0";

open( RUN, "$command |");

#### output is of type ####
#
##High Score list
##Protein ID                 DNA Str  ID                        Bits Evalue
#--------------------------------------------------------------------------
#Protein SH2                 DNA [+] EM:HSVAVONCO             96.62
#Protein SH3                 DNA [+] EM:HSVAVONCO             78.11


#NAME  ABC
#NAME  OR
#NAME  spinfull_2
#NAME  p450


while(<RUN>){
  chomp;
  my ( $protein_string, $hmm_id, $dna, $strand, $dna_id, $bit_score, $e_value) = split ;
  print STDERR $_."\n";
  if ($protein_string){
    next unless $protein_string eq 'Protein';
    next if ($hmm_id eq 'ID' || 
	     $hmm_id eq 'info');
    if ( $strand eq '[+]' ){
      $strand = +1;
    }
    else{
      $strand = -1;
    }
    
    unless ( $dna_id){
      $dna_id = '$slice';
    }
    print STDERR "$hmm_id\t$strand\t$dna_id\t$bit_score\n";
    
    #&store($db,$hmm_id,$strand,$dna_id,$bit_score);
  }
}

sub store{
  my ($db,$hmm_id,$strand,$dna_id,$bit_score) = @_;
  
  my $sql = qq(
	       INSERT into genewisedb ( protein_id, strand, dna_id, bit_score )
	       VALUES ( ?, ?, ?, ? ) 
	      );

  
  my $sth =  $db->prepare($sql);
  
  $sth->execute( $hmm_id,$strand,$dna_id,$bit_score);
}
 
