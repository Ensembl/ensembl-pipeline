#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;
use Bio::EnsEMBL::Pipeline::Flag;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::SeqIO;
use ncRNA_update_config ; 
my $dbhost;
my $dbuser;
my $pass;
my $dbname;
my $dbport;
my $chunk = 5;
my $help;
my %RFAM;

my @dbID;
my $count =0;
my $num=0;
my $start;
my $dln = 'DummyFlag';

my @multiexon_files;
my @species_list;
my $verbose;
my $norfam;


$| = 1; 

&GetOptions(
	    'pass=s'    => \$pass,
	    'verbose!'  => \$verbose, 
	    'species=s' => \@species_list,
	    'norfam!'      => \$norfam,
           );

if(!$pass || $help){
  die ("perl predict_ncRNA.pl 
-pass *(password) 
-species (species list) 
-verbose 
-norfam  (Only run miRNA annotation)
writes paths and rulemanager command lines to shell script species.csh(* required)\n");
  $help = 1;
}

my @speciess;
if (scalar(@species_list)) {
  @species_list = split(/,/,join(',',@species_list));
  foreach my $species (@species_list){
    if ($CONFIG->{$species}){
      push @speciess, $species;
    } else {
    print "Skipping species $species\n";
    }
  }
} else {
  @speciess = (keys %$CONFIG);
}

# might want to check some config first
SPECIES:  foreach my $species (@speciess){
   print "checking species $species\n";
  # check that all the dbs exist and have finished their run!
  my $sdb = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
    (
     -host   => $CONFIG->{$species}->{"WRITEHOST"},
     -user   => 'ensro',
     -port   => $CONFIG->{$species}->{"WRITEPORT"},
     -dbname => $CONFIG->{$species}->{"WRITENAME"},
    );
  my $dna_db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
    (
     '-host'   => $CONFIG->{$species}->{"DBHOST"},
     '-user'   => 'ensro',
     '-dbname' => $CONFIG->{$species}->{"DBNAME"},
     '-port'   => $CONFIG->{$species}->{"DBPORT"},
    );
  #add dna_db
   $sdb->dnadb($dna_db);
   my $sic = $sdb->get_StateInfoContainer;
   # assumes dbs were set up using the update_ncRNA script so the analysis ids aer correct.
   my $RFAM = $sic->list_input_ids_by_analysis('2');
   my $jobs = $sic->list_input_ids_by_analysis('7');
   my $miRBase = $sic->list_input_ids_by_analysis('3');
   my $lock = $sdb->get_meta_value_by_key('pipeline.lock');
   die "Pipeline is locked, rulemanager may still be running aborting...\n" if $lock;
   die "RFAM jobs not finished yet ".scalar(@$RFAM)." jobs finished out of ".scalar(@$jobs)." jobs submitted\n" unless (scalar(@$RFAM) == scalar(@$jobs) or $norfam);
   die "miRBase jobs not finished yet ".scalar(@$miRBase)." jobs finished out of ".scalar(@$jobs)." jobs submitted\n" unless (scalar(@$miRBase) == scalar(@$jobs));
   print " ok\n";
  }
 print "Checks complete - loading the next stage of input_ids...\n";

  SPECIES :foreach my $species (@speciess){
  my @input_ids;
  my @flags;
  print "Running species $species\n";
  print "opening connections to databases\n";
 
  my $sdb = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
    (
     -host   => $CONFIG->{$species}->{"WRITEHOST"},
     -user   => $WRITEUSER,
     -pass   => $pass,
     -port   => $CONFIG->{$species}->{"WRITEPORT"},
     -dbname => $CONFIG->{$species}->{"WRITENAME"},
    );
;
  my $fa = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($sdb);
  # dont rerun the filtering if it has already been done
  my @test_flags = @{$fa->fetch_all};
  unless (scalar(@test_flags > 0)){
    my $aa = $sdb->get_AnalysisAdaptor;
    my $sa = $sdb->get_SliceAdaptor;
    my $analysis = $aa->fetch_by_logic_name('ncRNA');
    my $mianalysis = $aa->fetch_by_logic_name('BlastmiRNA');
    my $RFanalysis = $aa->fetch_by_logic_name('RfamBlast');
    my $migoalanalysis = $aa->fetch_by_logic_name('miRNA');
    my $limit = 0;
    my %rfam_blasts;
    my %mirbase_blasts;
    my %rfam_threshold;
    my $dafa = $sdb->get_DnaAlignFeatureAdaptor;
    my @dafs;
    unless ($norfam){
      print "reading RFAM family data...\n";
      # fetch them by familly!!!!!!!!!
      my %thr;
      my @domcount =  @{sql("select count(*) , LEFT(hit_name,7) from dna_align_feature where  analysis_id = ". $RFanalysis->dbID .
			    " group by LEFT(hit_name,7) ;",$sdb)};
      foreach my $dom (@domcount){
	$thr{$dom->[1]} = $dom->[0];
      }
      
      my $total = scalar(keys (%thr));
      my $complete;
      my $last;
      print "\n0_________10________20________30________40________50________60________70________80________90_______100\n";
      foreach my $domain (sort keys %thr){
	$count ++;
	my @dafs;
	if ($thr{$domain} > 2000){
	  my @top_scores = @{sql("SELECT evalue FROM dna_align_feature WHERE analysis_id = ". $RFanalysis->dbID .
				 " AND LEFT(hit_name,7) = '".$domain."' ORDER BY  evalue ASC limit 2000;",$sdb)};
	  my $cutoff = pop (@top_scores)->[0];
	  @dafs = @{$dafa->generic_fetch("left(hit_name,7) = \"$domain\" AND evalue <= $cutoff")};
	} else {
	  @dafs = @{$dafa->generic_fetch("left(hit_name,7) = \"$domain\" ")};
	}
	$complete = int($count/$total*100);
	if ($complete > $last){
	  my $num = $complete -  $last;
	  foreach (my $i = 0; $i < $num; $i++){
	    print "=";
	  }
	}
	$last = $complete;
	@dafs = sort {$a->p_value <=> $b->p_value} @dafs if (scalar(@dafs) > 2000 );
      DAF:  foreach my $daf(@dafs){
	  next if ($daf->p_value > 0.01);
	  #    print $daf->p_value." ";
	  last DAF if ($rfam_blasts{$domain} &&  scalar(@{$rfam_blasts{$domain}}) >= 2000 );
	  push @{$rfam_blasts{$domain}},$daf;
	}
    }
      print "\nGenerating Flags\n";
      
      foreach my $domain(keys %rfam_blasts){
	my @hits = @{$rfam_blasts{$domain}};
	foreach my $hit (@hits){
	  my $flag = Bio::EnsEMBL::Pipeline::Flag->new
	    (
	     '-type'         => 'dna_align_feature',
	     '-ensembl_id'   => $hit->dbID,
	     '-goalAnalysis' => $analysis,
	    );
	  push @flags,$flag;
	}
      }
    }
    
    print "Filtering miRNAs";
    print "\n0_________10________20________30________40________50________60________70________80________90_______100\n";
    
    my %mi_types;
    my @domcount =  @{sql("select count(*) , hit_name from dna_align_feature where  analysis_id = ". $mianalysis->dbID .
			  " group by hit_name;",$sdb)};
    foreach my $dom (@domcount){
      $mi_types{$dom->[1]} = $dom->[0];
    }
   $count = 0;
    my $total = scalar(keys %mi_types);
    my $last = 0;
    foreach my $type (keys %mi_types){
      $count ++;
      if ($mi_types{$type} > 50){
	my @top_scores = @{sql("SELECT evalue FROM dna_align_feature WHERE analysis_id = ". $mianalysis->dbID .
			       " AND hit_name = '".$type."' ORDER BY  evalue ASC limit 50;",$sdb)};
	my $cutoff = pop (@top_scores)->[0];
	@dafs = @{$dafa->generic_fetch("hit_name = \"$type\" AND evalue <= $cutoff")};
      } else {
	@dafs = @{$dafa->generic_fetch("hit_name = \"$type\" ")};
      }
      my $complete = int($count/$total*100);
      if ($complete > $last){
	my $num = $complete -  $last;
	foreach (my $i = 0; $i < $num; $i++){
	  print "=";
	}
      }
      $last = $complete;
      @dafs = sort {$a->p_value <=> $b->p_value} @dafs if (scalar(@dafs) > 50 );
    DAF:  foreach my $daf(@dafs){
	#    print $daf->p_value." ";
	last DAF if ($mirbase_blasts{$type} &&  scalar(@{$mirbase_blasts{$type}}) >= 50 );
	push @{$mirbase_blasts{$type}},$daf;
      }
    }
    print "\nGenerating Flags\n";
    
    foreach my $domain(keys %mirbase_blasts){
      my @hits = @{$mirbase_blasts{$domain}};
      foreach my $hit (@hits){
	my $flag = Bio::EnsEMBL::Pipeline::Flag->new
	  (
	   '-type'         => 'dna_align_feature',
	   '-ensembl_id'   => $hit->dbID,
	   '-goalAnalysis' => $migoalanalysis,
	  );
	push @flags,$flag;
      }
    }

    print "Storing and randomizing flags so the input ids dont run in families\n";
    while (scalar(@flags >= 1)){
      my $random =  rand($#flags);
      # make sure you have an integer
      $random = int $random;
      $fa->store($flags[$random]);
      splice(@flags,$random,1);
    }
    $count =0;
    print "Making input ids for ncRNA\n";

    die("analysis object not found ncRNA\n") unless ($analysis);
    my @ids = @{$fa->fetch_by_analysis($analysis)};
    @ids = sort {$a->dbID <=> $b->dbID} @ids;
    foreach my $id (@ids){
      push @input_ids,$id->dbID;
    }

    my $inputIDFactory = new Bio::EnsEMBL::Pipeline::Utils::InputIDFactory
      (
       -db => $sdb,
       -top_level => 'top_level',
       -slice     => 1,
       -logic_name => 'DummyFlag',
      );
    $inputIDFactory->input_ids(\@input_ids);
    $inputIDFactory->store_input_ids;
  }
  # add the single input id to start the miRNA job
  # check its not already there first...
  my @check =  @{sql("select count(*) from input_id_analysis where input_id = \"BlastmiRNA\"",$sdb)};
  unless ($check[0]->[0]){
    sql_write(" insert into input_id_analysis(input_id,input_id_type,analysis_id,created) ".
	      "values(\"BlastmiRNA\",\"GENOME\",6,now());",$sdb);
  }



  print "Finished species $species\n";
  # do I want to start the rulemanager again at this point?
  my $perlpath = $ENV{"PERL5LIB"};

  open (SPEC,">$species.csh") or die ("Cannot open file $species.csh");
  print SPEC "#!/bin/csh\n\n";
  $ENV{"PERL5LIB"} = "$DATADIR/$species:$CVSDIR/ensembl-analysis/modules:$CVSDIR/ensembl-analysis/scripts:".
     "$CVSDIR/ensembl-pipeline/scripts:$CVSDIR/ensembl-pipeline/modules:".
      "$CVSDIR/ensembl/scripts:$CVSDIR/ensembl/modules:".
        "$BIOPERL_LIVE_PATH:$BIOPERL_RUN_PATH";
  print SPEC "setenv PERL5LIB ".$ENV{"PERL5LIB"}."\n";
  
  system ("perl $CVSDIR/ensembl-pipeline/scripts/setup_batchqueue_outputdir.pl"); 
  # if all be well, run the rulemanager
  my $cmd_rulemanager = "bsub -o $species.out -q normal perl $CVSDIR/ensembl-pipeline/scripts/rulemanager.pl ".
    "-dbname  $CONFIG->{$species}->{\"WRITENAME\"} ".
      "-dbport $CONFIG->{$species}->{\"WRITEPORT\"} ".
	"-dbhost $CONFIG->{$species}->{\"WRITEHOST\"} ".
	  "-dbuser ensadmin -dbpass $pass -once\n";

  print SPEC "$cmd_rulemanager\n";
  print "Monitor:\n";
  my $cmd = "perl $CVSDIR/ensembl-pipeline/scripts/monitor ".
    "-dbname  $CONFIG->{$species}->{\"WRITENAME\"} ".
      "-dbport $CONFIG->{$species}->{\"WRITEPORT\"} ".
	"-dbhost $CONFIG->{$species}->{\"WRITEHOST\"} ".
	  "-dbuser ensadmin -dbpass $pass -current";
  
  print "$cmd\n";
  # set it back to previous
  $ENV{"PERL5LIB"}= $perlpath;
}

exit;

sub sql {
  my ($query,$db) = @_;
  my @array;
#  print "QUERY\n$query\n";
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  while ( my @row = $sth->fetchrow_array ) {
    push @array , \@row;
  }
  return \@array;
}

sub sql_write {
  my ($query,$db) = @_;
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  return ;
}
__END__
