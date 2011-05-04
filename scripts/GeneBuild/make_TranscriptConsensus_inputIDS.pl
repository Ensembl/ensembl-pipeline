#!/usr/local/ensembl/bin/perl


=head1 NAME

ensembl-pipeline/scripts/GeneBuild/make_input_ids_for TranscriptConsensus

=head1 SYNOPSIS

this scripts will create input ids for the TranscriptConsensus analysis

=head1 DESCRIPTION

This scripts requires a TranscriptConsensus analysis and its corresponding Submit_ analysis
It will fetch the required input genes for each slice,
then create input ids of a specified size while making sure it doesn't cut any of the evidence used

=head1 OPTIONS

    -seq_region_names    list of seq_regions on which to run the script
                         if none specified, it runs on all the seq_regions

    -sr_names_file       the name of the file which contains a list of seq_region names on which
                         the script should be run.

    -dbname              name of the reference database containing the input_ids

    -dbhost              host of the reference database

    -dbpass              pass of the reference database

    -dbuser              user of the reference database

    -dbport              port of the reference database

    -outfile             output file, will contain the mysql insert command lines for the newly created inputs

    -coord_system        coord_system of the input seq_regions

    -path                coord_system version of the input seq_regions

    -slice_size          maximum size for the input_ids

    -analysis            logic_name of the TranscriptConsensus analysis
                         this will be used to fetch the corresponding configuration in the TranscriptConsensus config file

    -submit_analysis     logic_name of the associated Submit_analysis
                         the script will check the there is a rule that has logic_name has a goal and submit_analysis as a condition
                         the input_ids will then be added to the submit_analysis


=head1 EXAMPLES

perl make_TranscriptConsensus_inputIDS.pl -dbhost myhost -dbuser myuser -dbpass mypass -dbport 3306 -dbname mydb -submit_analysis SubmitConsensus -slice_size 1000000 -coord_system toplevel -outfile consensus_input_ids.sql -path mypath -analysis Sim_consensus

this commandline will write a file consensus_input_ids.sql containing the input_ids to be added to SubmitConsensus,
which is the dummy analysis for Sim_consensus

=cut


use strict;
use warnings;


use Getopt::Long;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;

use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus;
use Bio::EnsEMBL::Analysis::Config::Databases;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;

use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use ScriptUtils;




$| = 1 ;


my ( $outfile, $coordsystem, $dbhost, $dbname, $dbpass ) ;
my $dbuser = 'ensro' ;
my $dbport = 3306 ;
my $slice_size ;
my $analysis ;
my $submit_analysis ;
my @seq_region_names ;
my $sr_names_file;
my $path ;


&GetOptions(
            'seq_region_names:s'=> \@seq_region_names,
            'sr_names_file:s'   => \$sr_names_file,
            'dbname:s'          => \$dbname,
            'dbhost:s'          => \$dbhost,
            'dbpass:s'          => \$dbpass,
            'dbuser:s'          => \$dbuser,
            'dbport:s'          => \$dbport,
            'outfile:s'         => \$outfile,
            'coord_system:s'    => \$coordsystem,
            'path:s'            => \$path,
            'slice_size=i'      => \$slice_size,
            'analysis=s'        => \$analysis,
            'submit_analysis=s' => \$submit_analysis,
            );

if (scalar(@seq_region_names) && $sr_names_file) {
  throw("You cannot provide seq_region names via the command line and via a file at the same time.\n")
}


if (scalar(@seq_region_names)) {
  @seq_region_names = split(/,/, join(',',@seq_region_names) ) ;
} elsif ($sr_names_file) {
  open (IN, $sr_names_file) || die "Cannot open your seq_region_names file: $sr_names_file";
  my @names = <IN>;
  foreach my $name(@names) {
    chomp $name;
    push (@seq_region_names, $name);
  }
}

$slice_size = 100_000 unless $slice_size ;

my @all_biotypes ;  # array of all biotypes to cluster 


unless ($coordsystem) {
  print STDERR "you haven't supplied a coord_system_name using -coord_system --> I'll use 'chromosome' now\n" ;
  $coordsystem = 'chromosome' ;
}

unless ($outfile) {
  $outfile = "input_ids_" . $seq_region_names[0] . ".sql" ;
  print STDERR "you haven't supplied a name for output-fasta-file with -outfile--> I'll use '$outfile' now\n" ;
}

unless( $dbhost && $dbname &&  $outfile ) {
  print STDERR ("Usage: -dbname $dbname -dbhost $dbhost -dbuser $dbuser " . "-coord_system $coordsystem -outfile $outfile\n" ) ;
  exit(0) ;
}


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                           -host   => $dbhost,
                                           -user   => $dbuser,
                                           -port   => $dbport,
                                           -dbname => $dbname,
                                           -pass   => $dbpass,
                                          );

my $sa = $db->get_SliceAdaptor;


my $pa = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(
                                           -host   => $dbhost,
                                           -user   => $dbuser,
                                           -port   => $dbport,
                                           -dbname => $dbname,
                                           -pass   => $dbpass,
                                          );

my $rule = $pa->get_RuleAdaptor() ;
my $ana = $pa->get_AnalysisAdaptor->fetch_by_logic_name($analysis) ;
my $submit = $pa->get_AnalysisAdaptor->fetch_by_logic_name($submit_analysis) ;
my $input_id_type = $ana->input_id_type ;
my $analysis_id = $submit->dbID ;

throw ( "can't find analysis with logic_name $analysis in db $dbname \@ $dbhost\n")
  unless $ana ;
throw ( "can't find analysis with logic_name $submit_analysis in db $dbname \@ $dbhost\n")
  unless $submit ;
my $rule_id = $rule->fetch_by_goal($ana) ;
throw ( "Your submit analysis $submit_analysis doesn't match your logic_name $analysis\n")
  unless (${ $rule_id->list_conditions }[0] eq $submit_analysis) ;

my $chrhash = get_chrlengths_v20($db, $path, $coordsystem) ;

filter_to_chr_list(\@seq_region_names, $chrhash, $db->dbc->dbname) ;

my @input_ids ;

my %databases = %{$DATABASES} ;
my %database_hash ;

my $var_hash = $TRANSCRIPT_CONSENSUS_CONFIG_BY_LOGIC ;

my %check;
foreach my $k (keys %{$var_hash} ) {
  my $uc_key = uc($k);
  if (exists $check{$uc_key}) {
    throw("You have two entries in your config with the same name (ignoring case)\n");
  }
  $check{$uc_key} = $k;
}

foreach my $k (keys %check) {
  my $old_k = $check{$k};
  my $entry = $var_hash->{$old_k};
  delete $var_hash->{$old_k};
  $var_hash->{$k} = $entry;
}

my $uc_logic = uc($analysis);
my $entry ;
if (exists $var_hash->{$uc_logic}) {
  $entry = $var_hash->{$uc_logic};
}

my $input_dbs = $entry->{INPUT_GENES} ;
my $ab_initio = $entry->{AB_INITIO_LOGICNAMES} ;

# Check all databases declared in the Databases.pm config file exist
for my $category (keys %databases ) {

  if( !$databases{$category}{-host} || !$databases{$category}{-dbname} ) {
    print "Skipping category " . $category . " no dbinfo " . "defined\n" ;
    next ;
  }

  print "\n$category: $databases{$category}{-host} $databases{$category}{-dbname} :\n--------------------------------------\n" ;
  my %constructor_args = %{$databases{$category}} ;
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor ( %constructor_args, ) ;
  $database_hash{$category} = $dba ;
}

# Create input slices from the existing ones
SLICE:
foreach my $chr (@{sort_chr_names($chrhash)}) {
  my $chrstart = 1 ;
  my $chrend   = $chrhash->{$chr} ;

  print STDERR "Chr $chr from $chrstart to " . $chrend . "\n" ;

  my $slice ;
  my $slicename ;
  $slicename = "$coordsystem:$path:$chr:$chrstart:$chrend:1" ;
  $slice = $sa->fetch_by_name($slicename) ;

  my @all_genes ;

  my $coord_sys_name  =  $slice->coord_system_name ;
  my $coord_sys_version = $slice->coord_system->version ;
  my $seq_region_name = $slice->seq_region_name ;

  # dont construct anything if we have a slice shorter than the slice_size 
  if ($slice->length < $slice_size ) {
    push @input_ids , $slice->name ;
    next ;
  }

  for my $category ( keys %{$input_dbs} ) {
    my $dba = $database_hash{$category};
    my $tmp_slice = $dba->get_SliceAdaptor->fetch_by_name($slice->name) ;

    my $genes ;
    foreach my $biotype ( @{$input_dbs->{$category}} ) {
      push @all_biotypes, $biotype ;
      print "Fetching all genes from " . $tmp_slice->name . " by " . "biotype ".$biotype . " from " . $tmp_slice->adaptor->dbc->dbname . "\n";
      my @genes = @{ $tmp_slice->get_all_Genes_by_type($biotype,undef,1)} ;
      print "Have " . scalar(@genes) . " genes [$biotype]\n" ;
      push @all_genes, @genes ;
    }
    for ( @all_genes ) {
      print "gene_db_id : " . $_->dbID . "\t" . $_->biotype . "\n"  ;
    }

    foreach my $logic_name_becomes_biotype ( @{ $ab_initio } ) {
      push @all_biotypes, $logic_name_becomes_biotype ;
      my $pt = $tmp_slice->get_all_PredictionTranscripts( $logic_name_becomes_biotype, 1) ;
      print "Have " . scalar(@$pt) . " genes [$logic_name_becomes_biotype]\n" ;
      my $ab_initio_genes = convert_prediction_transcripts_to_genes($pt, $logic_name_becomes_biotype ) ;
      push @all_genes, @$ab_initio_genes ;
    }

  }

  if (@all_genes == 0) {
    warn("Slice ".$slice->name." has no genes ");
    next SLICE;
  }

  print "have " .scalar(@all_genes) . " genes for slice ".$slice->name." \n" ;
  print "having slice_size $slice_size\n" ;

  create_input_ids(\@all_genes, \@all_biotypes, $coord_sys_name, $coord_sys_version, $seq_region_name);
}


open(O,">$outfile") || die "cant write to file $outfile\n" ;
for my $id ( @input_ids ) {
  print O "insert into input_id_analysis (input_id, input_id_type, analysis_id) values ( \"$id\", \"$input_id_type\", $analysis_id ); \n" ;
}
close(O) ;
print "sql-statements written to $outfile - finished\n" ;


sub create_input_ids {
  my ($genes, $biotypes, $coord_sys_name, $coord_sys_version, $seq_region_name) = @_ ;
  my @all_genes = @$genes ;
  my @all_biotypes = @$biotypes ;

  @all_genes = sort { $a->seq_region_start <=> $b->seq_region_start }  @all_genes ;
  for ( @all_genes ) {
    print "gene_db_id : " . $_->dbID . "\t" . $_->biotype . "\t" .$_->seq_region_start . "\n"  ;
  }

  my %types_hash;
  $types_hash{all} = \@all_biotypes;
  my $clustered;
  my $unclustered;
  ($clustered,$unclustered)  =  cluster_Genes(\@all_genes, \%types_hash);
  my @cluster = @$clustered;
  push @cluster ,@$unclustered;

  @cluster = sort {$a->start <=> $b->start } @cluster ;
  print "Printing cluster info\n";
  for my $cluster (@cluster) {
    print $cluster->start . "---" . $cluster->end . "\n" ;
  }

  # Make blocks where clusters are

  my @cluster_blocks ;
  my $curblock = undef ;
  foreach my $cluster (@cluster) {
    if (defined($curblock) && $curblock->end >= $cluster->start) {
      if ($cluster->end > $curblock->end) {
        $curblock->end($cluster->end) ; 
      }
    } else {
      $curblock = Bio::EnsEMBL::SeqFeature->new(-START => $cluster->start, -END => $cluster->end) ;
      push @cluster_blocks,$curblock ;
    }
  }
  @cluster_blocks = sort {$a->start <=> $b->start} @cluster_blocks ;

  my $start ;
  my $nblock = scalar(@cluster_blocks) ;
  my $actual_max_size = 0 ;

  for (my $i=0 ; $i<@cluster_blocks ; $i++) {
    my $c = $cluster_blocks[$i] ;
    if (!defined($start)) {
      $start = $c->start ;
    }

    my $diff = $c->end - $start ;
    if ($diff > $slice_size || $i==$nblock-1) {
      my $end = $c->end ;
      my $id = $coord_sys_name . ":" . $coord_sys_version .":".$seq_region_name .  ":" .  $start . ":" . $end . ":1" ;
      print "$id\n" ;
      push @input_ids , $id ;
      if ($actual_max_size < $diff) {
        $actual_max_size = $diff ;
      }
      $start = undef;
    }
  }
  print "Maximum size of generated ids = $actual_max_size\n";
}


sub  convert_prediction_transcripts_to_genes {
  my ($pt, $logic_name_becomes_biotype ) = @_ ;
  my @new_genes ;
  for my $pt (@$pt) {
    # conversion 
    my $gene_from_pt = Bio::EnsEMBL::Gene->new(
                       -start    => $pt->start,
                       -end      => $pt->end,
                       -strand   => $pt->strand,
                       -slice    =>$pt->slice,
                       -biotype  => $logic_name_becomes_biotype,
                       -analysis =>$pt->analysis,
                       ) ;

    my $new_tr = Bio::EnsEMBL::Transcript->new(
                    -BIOTYPE  => $logic_name_becomes_biotype ,
                    -ANALYSIS => $pt->analysis ,
                 ) ;

    my @pt_exons  = @{ $pt->get_all_Exons} ;
    for (my $i = 0 ; $i < scalar(@pt_exons) ; $i++) {
      my $pte = $pt_exons[$i] ;
      bless $pte, "Bio::EnsEMBL::Exon" ;
      $pte->phase(0) ;
      $pte->end_phase(0) ;
    }
    for (@pt_exons) {
      $new_tr->add_Exon($_) ;
    }
    $gene_from_pt->add_Transcript($new_tr) ;

    push @new_genes, $gene_from_pt ;
  }
  return \@new_genes ;
}

