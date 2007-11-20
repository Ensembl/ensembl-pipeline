#!/usr/local/ensembl/bin/perl

# "insert into seq_region_attrib (seq_region_id,attrib_type_id,value) values (seq_region_id_of_mitochondrion,11,5)"

# based on ensembl branch 25
use strict;
use FlyBaseGff;
use FlyBaseConf;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

$| = 1;

my $load_genes;
my $load_simple_features;

&GetOptions(
  'load_genes' => \$load_genes,  # flag
  'load_simple_features' => \$load_simple_features,  # flag
);

if (!$load_genes && !$load_simple_features) {
  throw("Must specify whether you are loading genes or simple_features");
}

print "reading config-file FlyBaseConf.pm\n";

# the data base that we will write to
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host => $FB_DBHOST,
                                            -user => $FB_DBUSER,
                                            -dbname => $FB_DBNAME,
                                            -pass  => $FB_DBPASS,
                                            -port  => $FB_DBPORT,
                                           );


# # #
# check if all files are there and readable
# # #
for my $conf(@{$FB_CHR_INFO}){
     unless ( -r  $conf->{'gff_file'}){
	print "Can't find file ".$conf->{'gff_file'}."\n";
	exit(0);
     }else{
	print "File ".$conf->{'gff_file'}." is OK\n";
     }
}

# # #
# loop through each chromosome's gff file as listed as in Configuration-file
# # #
foreach my $chr(@{$FB_CHR_INFO}) {
  print STDERR ">> handling ".$chr->{'chr_name'}." with file ".$chr->{'gff_file'}."\n" if($FB_DEBUG);
  my $slice= $db->get_SliceAdaptor()->fetch_by_region('chromosome',$chr->{chr_name});

  my $gff = FlyBaseGff->new(
                            -DB => $db,
                            -GFF=> $chr->{'gff_file'},
                            -SLICE => $slice
                           );

  # # #
  # process genes of gff
  # # #

  # note that the mt:ori gene is not stored because it has no transcripts
  if ($load_genes) {
  $gff->store_as_gene_object("gene","mRNA","protein_coding");
  $gff->store_as_gene_object("gene","ncRNA","ncRNA");
  $gff->store_as_gene_object("gene","snRNA","snRNA");
  $gff->store_as_gene_object("gene","tRNA","tRNA");
  $gff->store_as_gene_object("gene","rRNA","rRNA");
  $gff->store_as_gene_object("gene","pseudogene","pseudogene");
  $gff->store_as_gene_object("gene","snoRNA","snoRNA");
  $gff->store_as_gene_object("gene","miRNA","miRNA");

  } 
  # # #
  # store all simplefeatures as referenced in the FlyBaseConf.pm
  # # #
  if ($load_simple_features) {
    my $sa = $db->get_SimpleFeatureAdaptor();
  
    foreach my $feat (@{$SIMPLE_FEATURES}) {
      my $feature_type = $feat->{type};
      print "searching for simple features of type $feature_type ... ";
      my $sf_num = $gff->store_as_simple_feature($sa, $feature_type);
      
      if ($sf_num > 0 ) { 
        print "number found + stored in db: $sf_num \n";
      }else { 
        print " NOTHING found\n" ; 
      }
    }
  }
}



# # # 
# dump all loaded proteins
# # #

 # this script dumps out the translation of a transcript, given the transcript's stable_id or dbID
  my $cmd = "perl $FB_CVS_DIR/ensembl-pipeline/scripts/protein_pipeline/" .
            "dump_translations.pl -dbname $FB_DBNAME -dbhost $FB_DBHOST -dbport $FB_DBPORT -dbuser ensro " .
	    "-stable_id -db_id -file $FB_DUMPED_TRANSLATIONS_FILE"; 


  print "CMD is:\n$cmd\n";
  system ("$cmd") ;

print "flybase_to_ensembl is DONE ! \n";
