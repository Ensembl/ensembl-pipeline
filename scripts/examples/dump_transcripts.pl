#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/examples/dump_transcripts.pl,v $
# $Revision: 1.2 $

# # # 
# You'll need bioperl and the ensembl core checkout in your PERL5LIB
#
# Example commandline
#   perl dump_transcripts.pl -species human -logicname adipose_rnaseq -dbtype rnaseq
#   perl dump_transcripts.pl -species 'Mus musculus' -logicname ensembl -dbtype core 
# # # 

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

# database
my $species;
my $host = 'ensembldb.ensembl.org';
my $user = 'anonymous';
my $port = 5306;

# genomic location
my $logic_name;

# database type to fetch data from
my $db_type;

# files to write output to
my $fasta_output_file;
my $info_output_file;

# usage
my $help = '';
$|=1;

if ( !GetOptions( 'logicname|l=s'  => \$logic_name,
                  'species|s=s'    => \$species,
                  'dbtype|d=s'     => \$db_type,
                  'fasta|f=s'      => \$fasta_output_file,
                  'info|i=s'       => \$info_output_file,
                  'help|h!'        => \$help )
     || !( defined($logic_name) && defined($species) && defined($db_type))
     || $help )
{
  print <<END_USAGE;

Usage:
  $0 --species=species --logicname=logic_name --dbtype=db_type --fasta=transcript_outfile.fa --info=transcript_outfile.txt
  $0 --help

    --species / -s      Name of species

    --logicname / -l    Logic_name for Gene models 

    --dbtype / -d       Database type, eg. core, otherfeatures, rnaseq

    --fasta / -f        Fasta file containing the transcript sequences

    --info / -i         Text file containing exon coordinates for each transcript

    --help / -h         To see this text

Example usage:

  $0 -s human -l brain_rnaseq -d rnaseq -f transcripts.fa -i transcripts.txt

or 

  $0 -s mouse -l ensembl -d core -f transcripts.fa -i transcripts.txt

END_USAGE

  exit(1);
} ## end if ( !GetOptions( 'logicname|c=s'...))

# # # usage

# connect to database:
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db( '-host' => $host, 
                                  '-port' => $port,
                                  '-user' => $user, );

# get adaptors
my $gene_adaptor = $registry->get_adaptor( $species, $db_type, 'Gene' );
my $slice_adaptor = $registry->get_adaptor( $species, $db_type, 'Slice' );
if (defined $gene_adaptor && defined $slice_adaptor) {
  print "\nConnected to $species database\n\n";
} else {
  print "\nERROR: not connected to $species database\n\n";
}

# open the output files
if (!defined $fasta_output_file) {
  $fasta_output_file = $gene_adaptor->db->dbc->dbname."_".$logic_name."_transcripts.fa";
}
open(FASTA, ">$fasta_output_file") or die ("Can't read $fasta_output_file $! \n");
if (!defined $info_output_file) {
  $info_output_file = $gene_adaptor->db->dbc->dbname."_".$logic_name."_transcripts.txt";
}
open(INFO, ">$info_output_file") or die ("Can't read $info_output_file $! \n");
 

# we need to go through the slices
my $slices = $slice_adaptor->fetch_all('toplevel',undef, 1); 

foreach my $slice (@$slices) {
  # fetch genes from rnaseq
  my $genes = $slice->get_all_Genes($logic_name);
  print STDERR "Fetched ".scalar(@$genes)." genes from slice ".$slice->name."\n";
  # Loop through genes and print gene and transcript info 
  foreach my $gene (@$genes) {
    print INFO "GENE ".$gene->stable_id." chr ".$gene->slice->seq_region_name." start ".$gene->start." end ".$gene->end." strand ".$gene->strand."\n";
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      print INFO "  TRANSCRIPT ".$transcript->stable_id." chr ".$transcript->slice->seq_region_name." start ".$transcript->start." end ".$transcript->end." strand ".$transcript->strand."\n";
      my $num_exons=0;
      foreach my $exon (@{$transcript->get_all_Exons}) {
        $num_exons++;
        print INFO "    EXON ".$num_exons." chr ".$exon->slice->seq_region_name." start ".$exon->start." end ".$exon->end." strand ".$exon->strand."\n";
      } # exon loop
      # print fasta sequence for transcript
      print FASTA ">".$transcript->stable_id."\n".$transcript->seq->seq."\n";
    } # transcript loop
  } # gene loop
} # slice loop
print STDERR "DONE!\n";
