#!/usr/local/ensembl/bin/perl  -w

# based on ensembl branch 25
use strict;
use FlyBaseGff;
use FlyBaseConf;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor


$| = 1;



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host => $FB_DBHOST,
                                            -user => $FB_DBUSER,
                                            -dbname => $FB_DBNAME,
                                            -pass  => $FB_DBPASS,
                                            -port  => $FB_DBPORT,
                                           );


my $sa = $db->get_SimpleFeatureAdaptor();



# loop through each chr-gff in Configuration-file

foreach my $chr(@{$FB_CHR_INFO}) {
  print "handling ".$chr->{'chr_name'}." with file ".$chr->{'gff_file'}."\n" if($FB_DEBUG);

  my $slice= $db->get_SliceAdaptor()->fetch_by_region('chromosome',$chr->{chr_name});


  my $gff = FlyBaseGff->new(
                            -DB => $db,
                            -GFF=> $chr->{'gff_file'},
                            -SLICE => $slice
                           );



  #
  # process genes of gff
  #
  # ###############################################################################

  $gff->store_genes();







  #
  # store all simplefeatures as referenced in the FlyBaseConf.pm
  #
  # ###############################################################################

  foreach my $feat (@{$SIMPLE_FEATURES}) {
    my $feature_type = $feat->{type};
    my $feature_label = $feat->{label};
    my $logic_name = $feat->{logic_name};
    print "searching for simple features of type $feature_type ... ";

    my $sf_num = $gff->store_as_simple_feature($sa, $feature_type,  $logic_name, $feature_label);
    print "$sf_num features found and stored in db.\n";
  }


}


































