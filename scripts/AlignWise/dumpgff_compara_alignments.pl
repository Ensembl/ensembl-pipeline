#!/usr/local/ensembl/bin/perl -w

############################################################
# script to get a multiple alignment from a set of species
# given a focus species and a region in the species
#
# The alignment is localised in the focus species, so that
# there are no gaps in the focus sequence.
#
# the multiple alignment is used to run alignwise

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::Seq;
use Bio::EnsEMBL::Pipeline::Runnable::AlignWise;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $verbose = 0;
my $input_id;
my $alignment_type = "WGA";
my $focus_species  = 'Homo sapiens';
my $compara_host   =  'ecs2f';
my $compara_name   = 'ensembl_compara_16_1';
my $compara_port   = '3306';
my $compara_conf   = '/nfs/acari/eae/ensembl/ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf';
my @target_species = ( 'Mus musculus', 'Rattus norvegicus');

GetOptions(
	   "compara_host=s"   => \$compara_host,
	   "compara_port=i"   => \$compara_port,
	   "compara_name=s"   => \$compara_name,
	   'input_id:s'       => \$input_id,
	   "alignment_type=s" => \$alignment_type,
	   "focus_species=s"  => \$focus_species,
	   "compara_conf=s"   => \$compara_conf,
	  );


############################################################
# get chr_name start and end from input_id
unless( $input_id ){
  &usage;
  exit(0);
}
$input_id =~/(\S+)\.(\d+)-(\d+)/;
my $chromosome = $1;
my $start = $2;
my $end = $3;

############################################################
# use compara database
my $compara = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(-host      => $compara_host,
							  -user      => 'ensro',
							  -dbname    => $compara_name,
							  -port      => $compara_port,
							  -conf_file => $compara_conf,
							 );

############################################################
# you will work with as many target species as you
# put in your compara_conf file from the specified ones in @target_species
# Check that they all actually have alignment information


############################################################
# get genome dbs from compara
my $all_genome_dbs = $compara->get_GenomeDBAdaptor()->fetch_all();

############################################################
# the focus genome db
my ( $focus_genomedb ) = grep { $_->name() eq $focus_species } @$all_genome_dbs;

unless( $focus_genomedb ){
  print STDERR "no genome db found for focus species: $focus_species\n";
  &usage();
  exit(0);
}

############################################################
# the other genome dbs
my ( @other_genomedb )   = grep { $_->name() ne $focus_species } @$all_genome_dbs;

############################################################
# check that the dbs are in compara_conf
my @wanted_other_genomes;
foreach my $other_genome (@other_genomedb) {
  my $name     = $other_genome->name;
  my $assembly = $other_genome->assembly;
  if ( grep { $_ eq $name } @target_species ){
    if (defined $compara->get_db_adaptor($name,$assembly)) {
      push @wanted_other_genomes, $other_genome;
    }
  }
}

my @target_genomes = @wanted_other_genomes;

############################################################
# focus_genomedb is a GenomeDB object
my $focus_slice =  $focus_genomedb->db_adaptor->get_SliceAdaptor
  ->fetch_by_chr_start_end( $chromosome, $start, $end );

############################################################
# get the dna_frags of chromosome type for the focus_genome
my $dnafrags = $compara->get_DnaFragAdaptor()->fetch_all_by_GenomeDB_region( 
									    $focus_genomedb,
									    "Chromosome",
									    $chromosome,
									    $start,
									    $end
									   );

my $gaa = $compara->get_GenomicAlignAdaptor;

my $count = 1;
my $alignments_found = 0;



############################################################
# go over each of the target species:
 TARGET:
for( my $i_species = 0; $i_species<=$#target_genomes; $i_species++ ) {
    my $qy_gdb = $target_genomes[$i_species];
    
    # go over each dna_frag of the focus species:
  FRAG:
    foreach my $df ( @$dnafrags ) {
	print STDERR "dnafrag: ".$df->start."-".$df->end."\n" if $verbose;
	
	############################################################
	#calculate coords relative to start of dnafrag
	my $df_start = $start - $df->start + 1;
	my $df_end   = $end   - $df->start + 1;
	
	############################################################
	# constrain coordinates so they are completely within the dna frag
	my $len = $df->end - $df->start + 1;
	$df_start = ($df_start < 1)  ? 1 : $df_start;
	$df_end   = ($df_end > $len) ? $len : $df_end;
	print STDERR "coordinates relative to start of dnafrag: $df_start - $df_end\n" if $verbose;
	
	############################################################
	# fetch all alignments (HSPs as stored in compara) 
	# in the region we are interested in
	my $genomic_aligns = $gaa->fetch_all_by_DnaFrag_GenomeDB($df,
								 $qy_gdb,
								 $df_start,
								 $df_end,
								 $alignment_type
								 );
	
	print STDERR "genomic aligns for ".$qy_gdb->name."\n";
	
      ALIGN:
	for my $align ( @$genomic_aligns ) {
	    print STDERR 
		"consensus_start: ".$align->consensus_start." ".
		    "consensus_end: ".$align->consensus_end." ".
			"query_start: ".$align->query_start." ".
			    "query_end: ".$align->query_end." ".
				"query_strand: ".$align->query_strand." ".
				    "alignment_type: ".$align->alignment_type." ".
					"score: ".$align->score." ".
					    "perc_id: ".$align->perc_id." ".
						"cigar_line: ".$align->cigar_line."\n" if $verbose;
	    
	    
	    # a dna_frag from the focus region
	    my $cfrag = $align->consensus_dnafrag();
	    
	    # a dna_frag from the target region
	    my $qfrag = $align->query_dnafrag();
	    
	    ############################################################
	    # get the aligned string such that gaps in the focus seq are removed
	    # and the corresponding insertions in the target seq are deleted as well
	    my $aligned_string = 
		$align->sequence_align_string( 
					       $cfrag->contig(), 
					       $qfrag->contig(),
					       "QUERY", 
					       "FIX_CONSENSUS" 
						   );
	    if( ! $aligned_string ) { 
		next; 
	    } 
	    else {
		$alignments_found = 1; 
	    }
	    ############################################################
	    # calculate the positions of the target aligned fragment in the
	    # focus sequence fragment local coordinate system
	    #
	    # cfrag->start             : start of the frag in chromosome coord. (usually  = 1 )
	    # target_slice->chr_start : start of the target slice in chr coord
	    # align->consensus_start   : position of this align-piece in the target dna_frag (chr) coord.
	    
	    my ( $p_start, $p_end );
	    $p_start = $align->consensus_start() + $cfrag->start() - $focus_slice->chr_start();
	    $p_end   = $align->consensus_end()   + $cfrag->start() - $focus_slice->chr_start();
	    
	    ############################################################
	    # print GFF
	    
	    my $homology_region_id = ($align->query_dnafrag->name).".".($align->query_start)."-".($align->query_end);

	    my @parts = split /\s+/, $qy_gdb->name;
	    my $name = join '_', @parts;
	    
	    my $strand = '+';
	    if ( $align->query_strand == -1 ){
		$strand = '-';
	    }

	    print STDOUT 
		$count."\t".
		    $name."\t".
			"homology\t".
			    ($p_start) . "\t" . 
				($p_end) . "\t" .
				    $align->score. "\t" . 
					$strand. "\t" . 
					    ".". "\t" . 
						$homology_region_id."\n";

	    
	    # put the align string into result
	    $count++;
	    
	} # end of ALIGN
    }   # end of FRAG
}     # end of TARGET


if( ! $alignments_found ) {
  print STDERR "Sorry, no alignments found in the specified region.\n";
  exit(0);
} 
else {
  print STDERR ($count - 1)." alignments found\n";
}

############################################################

sub usage {
  print STDERR <<EOF

Usage: dumpgff_compara_alignments [options]
  where options should be 
  
  -compara_host
  -compara_port
  -compara_name
  -compara_conf   : filename to specify compara database connection
  
  -input_id       : chr_name.start-end 
  -alignment_type : default WGA 
  -focus_species 

  to specify the region you want to dump with alignments to other
  available species.
      
  -focus_species can be one of these (between quotes)
    "Homo sapiens"
    "Mus musculus"
    "Rattus norvegicus"
    "Fugu rubripes"
    "Anopheles gambiae"
    "Drosophila melanogaster"
    "Caenorhabditis elegans"
    "Caenorhabditis briggsae"


EOF
}
