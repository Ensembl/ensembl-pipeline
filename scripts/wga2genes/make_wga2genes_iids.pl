#!/usr/local/ensembl/bin/perl

# Generation of input ids for WGA2Genes. Segments the
# genome into regions that contain contain genes. 
# Nominally, each iid wil be a slice that contains a single 
# gene, but my contain more than one gene if genes overlap


use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $gene_dbname,
    $gene_dbhost,
    $gene_dbport,
    $gene_dbuser,
    $gene_dbpass,
    $source_type,
    $target_logic,
    $write,
);

$dbuser = 'ensro';
$dbport = 3306;
$gene_dbuser = 'ensro';
$gene_dbport = 3306;

$source_type = 'protein_coding';

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'genedbname=s' => \$gene_dbname,
            'genedbhost=s' => \$gene_dbhost,
            'genedbport=s' => \$gene_dbport,
            'genedbuser=s' => \$gene_dbuser,
            'genedbpass=s' => \$gene_dbpass,
            'targetlogic=s' => \$target_logic,
            'write' => \$write,
            );

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
);

my $gene_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
	'-dbname' => $gene_dbname,
	'-host' => $gene_dbhost,
	'-port' => $gene_dbport,
	'-user' => $gene_dbuser,
        '-pass' => $gene_dbpass,
        );

my $ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($target_logic);
if (not defined $ana) {
  die "Could not find analysis with logic name '$target_logic' in pipe db\n";
}


my @gene_clusters;

my @slices = @{$gene_db->get_SliceAdaptor->fetch_all('toplevel')};
foreach my $sl (@slices) {
  next if $sl->seq_region_name =~ /^MT$/;

  my @genes = @{$sl->get_all_Genes};
  @genes = grep { $_->biotype eq $source_type } @genes;
  @genes = sort { $a->start <=> $b->start } @genes;
  
  my $new_slice = 1;
  foreach my $g (@genes) {
    my $g_start = $g->start + $sl->start - 1;
    my $g_end   = $g->end   + $sl->start - 1;

    if ($new_slice or
        $gene_clusters[-1]->{end} < $g_start) {
      push @gene_clusters, {
        coordsys => $sl->coord_system->name,
        name  => $g->slice->seq_region_name,
        start => $g_start,
        end   => $g_end,
        genes => [$g],
      };
      $new_slice = 0;
    } else {
      push @{$gene_clusters[-1]->{genes}}, $g;
      if ($gene_clusters[-1]->{end} < $g_end) {
        $gene_clusters[-1]->{end} = $g_end;
      }
    } 
  }
}

foreach my $c (@gene_clusters) {
  my $iid = sprintf("%s::%s:%d:%d:", 
                    $c->{coordsys},
                    $c->{name},
                    $c->{start},
                    $c->{end});

 if ($write) {
   my $s_inf_cont = $db->get_StateInfoContainer;
   
   eval {
     $db->get_StateInfoContainer->
         store_input_id_analysis( $iid, $ana, '' );
   };
   if ($@) {
     print STDERR "Input id $iid already present\n";
   } else {
     print STDERR "Stored input id $iid\n";
   }       
 } else {
    printf("INSERT into input_id_analysis(input_id, input_id_type, analysis_id) " . 
           " values('%s', '%s', %d);\n", 
           $iid, 
           $ana->input_id_type,
           $ana->dbID);
    
  }
}


