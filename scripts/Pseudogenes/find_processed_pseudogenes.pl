#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $dbhost    = 'ecs2d';
my $dbuser    = 'ensro';
my $dbname    = 'homo_sapiens_vega_12_31';
my $dbpass    = undef;
my $dnadbhost;
my $dnadbname;



# logic names of the protein evidence:
my $targetted_protein  = "rodent-protein";
my $similarity_protein = "other-protein";

my $genetype;

$dbuser = "ensro";
&GetOptions(
	    'dbname:s'    => \$dbname,
	    'dbhost:s'    => \$dbhost,
	    'dnadbname:s'    => \$dnadbname,
	    'dnadbhost:s'    => \$dnadbhost,
	    'genetype:s'  => \$genetype,
);

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost ){
    print STDERR "Usage: $0 -dbname -dbhost -dnadbhost -dnadbname [ -genetype]\n";
    exit(0);
}

my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dbpass,
					       );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					    );


print STDERR "connected to $dbname : $dbhost\n";

if ($genetype){
  print STDERR "looking for pseudogenes in $genetype genes\n";
}

my %list_per_hid;

my $gene_adaptor  = $db->get_GeneAdaptor;
my $slice_adaptor = $db->get_SliceAdaptor;

my $count = 0;


############################################################
# first classify the transcripts according to their evidence
############################################################

my %chr_name;


GENE:
foreach my $gene_id ( @{$db->get_GeneAdaptor->list_geneIds} ){
  
  my $gene;
  eval{
    $gene = $gene_adaptor->fetch_by_dbID($gene_id);
  };
  unless ($gene){
    print STDERR "could not retrieve gene\n";
    next GENE;
  }
  $count++;
  #last GENE if $count>100;
  
  if ( $genetype ){
    next unless ( $genetype eq $gene->type );
  }

 TRAN:
  foreach my $trans ( @{$gene->get_all_Transcripts} ) {
    
    my @evidence = &get_evidence( $trans );
    if ( @evidence ){
      my $slice    = $slice_adaptor->fetch_by_transcript_id($trans->dbID);
      my $fakegene = Bio::EnsEMBL::Gene->new();
      $fakegene->add_Transcript( $trans );
      my $tmp_gene = $fakegene->transform( $slice );
      my @transcripts = @{$tmp_gene->get_all_Transcripts};
      
      ############################################################
      # note that a transcript can appear in more than one list
      ############################################################
      $chr_name{$transcripts[0]} = $tmp_gene->chr_name;
      foreach my $evidence ( @evidence ){
	push ( @{$list_per_hid{$evidence}}, $transcripts[0] );
      }
    }
  }
  
}

#print STDERR "$count genes retrieved\n";

my %pseudos_per_id;
my %parent_transcript;

############################################################
# go through each list and get those transcripts that
# are not spliced and have a transcript based on the same evidence
# which are spliced
############################################################

EVIDENCE:
foreach my $protein_id ( keys %list_per_hid ){
  my @transcripts = 
    sort { scalar(@{$b->get_all_Exons}) <=> scalar( @{$a->get_all_Exons} ) } @{$list_per_hid{$protein_id}};
  
  # Is the first one spliced?
  my $first = shift @transcripts;
  if ( &is_spliced( $first ) ){
    
    $parent_transcript{ $protein_id } = $first;
    #take those that are not spliced
    foreach my $t ( @transcripts ){
      unless ( &is_spliced( $t ) ){
	push ( @{$pseudos_per_id{ $protein_id }}, $t );
      }
    }
  }
}


my %pseudogenes_per_id;

############################################################
# compare each potetial pseudogene translation with
# the translation of the parent transcript (the one which is spliced and based on the same
# evidence. Accept if they have similar translation.
############################################################

PSEUDOS:
foreach my $protein_id ( keys %parent_transcript ){
  my $parent = $parent_transcript{$protein_id};
  foreach my $pseudo ( @{$pseudos_per_id{ $protein_id } } ){
    #if ( &compare_translations( $parent, $pseudo ) ){
    push ( @{$pseudogenes_per_id{$protein_id}}, $pseudo );
    my $id = $pseudo->stable_id || $pseudo->dbID;
    print $protein_id."\t".$pseudo->stable_id."\t".$pseudo->dbID."\t".$chr_name{$pseudo}."\t".$pseudo->start."\t".$pseudo->end."\n";
    #}
  }
}




############################################################

sub compare_translations{
  my ( $t1, $t2 ) = @_;
  my $translation1;
  my $translation2;
  
  eval {
    $translation1 = $t1->translate;
  };
  if ($@) {
    print STDERR "Couldn't translate transcript ".$t1->dbID."\n";
    return 0;
  }
  eval {
    $translation2 = $t2->translate;
  };
  if ($@) {
    print STDERR "Couldn't translate transcript ".$t2->dbID."\n";
    return 0;
  }
  if ( $translation1 && $translation2 ){
    my $seq1  = $translation1->seq;
    my $seq2 = $translation2->seq;
    
    if($seq1 eq $seq1) {
      return 1;
    }
    elsif($seq1 =~ /$seq2/){
      return 1;
    }
    elsif($seq1 =~ /$seq2/){
      return 1;
    }
  }
  
  return 0;
}


############################################################

sub get_evidence{
  my $t = shift;
  my %evidence;
  foreach my $exon ( @{$t->get_all_Exons} ){
    foreach my $evi ( @{$exon->get_all_supporting_features} ){
      $evidence{  $evi->hseqname } = 1;
      #if ( $evi->analysis->logic_name eq $targetted_protein || 
      #	 $evi->analysis->logic_name eq $targetted_protein ){
      #	$evidence{  $evi->hseqname } = 1;
      #    }
    }
  }
  my @ids = keys %evidence;
  return @ids;
}

############################################################
# we consider the transcript not-spliced if it has only one
# exon or whether it has only small introns <9, that is, frameshifts
# which are also indicators of a pseudogene.
############################################################

sub is_spliced{
  my $t = shift;
  my @exons = @{$t->get_all_Exons};
  if ( scalar (@exons ) == 1 ){
    return 0;
  }
  elsif( scalar (@exons) > 1 ){
    
    # check that there are not funky frame shifts
    @exons = sort{ $a->start <=> $b->start } @exons;
    for(my $i=0; $i<$#exons; $i++){
      my $intron = $exons[$i+1]->start - $exons[$i]->end - 1;
      if ( $intron > 9 ){
	return 1;
      }
    }
    return 0;
  }
  else{
    return 0;
  }
}

############################################################
