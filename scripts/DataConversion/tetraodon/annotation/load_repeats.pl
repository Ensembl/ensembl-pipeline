#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::RepeatConsensus;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $gff_feature, 
    $gff_source
);


&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'gff_source=s' => \$gff_source,
            'gff_feature=s' => \$gff_feature,
            );

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
	'-pass' => $dbpass
);



# make the analysis objects
my $ana_adaptor = $db->get_AnalysisAdaptor;

my $ana_obj;

my (%analyses, %repeat_consensi, %slices, @all_repeats);

while (<>) {
  /^\#/ and next; 

  my @line = split /\t/;
  if ((not defined($gff_source) or $line[1] eq $gff_source) and
      (not defined($gff_feature) or $line[2] eq $gff_feature)) {

    my $ana_obj;

    if (not exists($analyses{$line[1]})) {
      my $logic_name = $line[1];

      # need to make the analysis object; try fetching from the database;
      if (not defined($ana_obj = $ana_adaptor->fetch_by_logic_name($logic_name))) {
        $ana_obj  = Bio::EnsEMBL::Analysis->new(
                                                -logic_name      => $logic_name,
                                                -gff_source      => $logic_name,
                                                -gff_feature     => 'Repeat');
        # the repeat feature adaptor will not do this for us
        $db->get_AnalysisAdaptor->store($ana_obj);
      }
      $analyses{$logic_name} = $ana_obj;
    } else {
      $ana_obj = $analyses{$line[1]};
    }

    my ($repeat_name)  = $line[8] =~ /repeat_name\s+\"([^\"]+)\"/;
    my ($repeat_class) = $line[8] =~ /repeat_class\s+\"([^\"]+)\"/;

    my ($cons_obj, $slice);

    if (not ($cons_obj = $repeat_consensi{"$repeat_class.$repeat_name"})) {
      $cons_obj = Bio::EnsEMBL::RepeatConsensus->new;
      $cons_obj->name($repeat_name);
      $cons_obj->repeat_class($repeat_class);
      $repeat_consensi{"$repeat_class.$repeat_name"} = $cons_obj;
    }

    if (not ($slice = $slices{$line[0]})) {
      $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', $line[0]);
      $slices{$line[0]} = $slice;
    }
    
    my $rf = Bio::EnsEMBL::RepeatFeature->new;
    $rf->score($line[5]);
    $rf->start($line[3]);
    $rf->end($line[4]);
    $rf->strand($line[6] eq "-" ? -1 : 1);
    $rf->hstart(0);
    $rf->hend(0);
    $rf->repeat_consensus($cons_obj);
    $rf->slice($slice);
    $rf->analysis($ana_obj);

    push @all_repeats, $rf;


  }
}



print STDERR "Created ", scalar(@all_repeats), " features; writing...\n";
$db->get_RepeatFeatureAdaptor->store(@all_repeats);

