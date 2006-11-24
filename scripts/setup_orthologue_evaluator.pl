#!/usr/local/ensembl/bin/perl -w

=pod

=head1 NAME

  setup_orthologue_evaluator.pl 



=head1 SYNOPSIS

perl setup_orthologue_evaluator.pl -dbname <dbname> -dbuser -write 


=head1 DESCRIPTION

 
 This script can be usesd to setup the 3 analysis in an easy way. 
 It adds the initial Analysis you specify in OrthologueEvaluator.pm configuration file and 
 adds the post-analysis and rules /conditions as well. 

 The script will re-write your Exonerate2Genes.pm file if you use the -write option. 

 This script depends on PPI and List::MoreUtils 

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'

=cut

 

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry; 

use FindBin;  
use lib "$FindBin::Bin" ;   
use Data::Dumper;  
use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluatorExonerate;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;
use Bio::EnsEMBL::Utils::Exception qw(info verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::ConfigUtils; 
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Analysis; 
use Bio::EnsEMBL::Pipeline::Rule; 
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory; 
verbose('WARNING') ; 
$! = 1;

my %opt = (
    dbname => '',
    dbhost => '',
    dbuser => 'ensadmin',
    dbpass => 'ensembl',
    dbport => '3306',
    exonerate_file => 'exonerate-1.0.0',
);

&GetOptions(
    \%opt,
    'dbname=s',
    'dbuser=s',
    'dbhost=s',
    'dbport=s',
    'dbpass=s',
    'verbose+',  
    'exonerate_file=s', # specify the name of the exonerate program_file for analysis table 
                        # in /usr/local/bin/ensembl/bin if don't want to use the default 
                        # (exonerate-1.0.0)

    'write',   # if you want to write the Exonerate2Genes file  
); 

#
# The names of the initial analysis are hard-coded in sub get_pre_analysis.  
#  
# The are :  
#
#  MissingOrtholouges
#  FindPartialGenes
#  FindSplitGenes 
#

unless ($opt{dbname} && $opt{dbhost} ) { 
  throw("You need to run the script with db-parameters to connect to your reference database\n".
        "Please specify -dbport as well if it's not 3306\n".
        "perl setup_orthologue_analysis -dbname <REFDB> -dbhost <SERVER> -dbuser <USER> -dbpass <PASS> -dbport <3306>\n"  ) ;   
}


Bio::EnsEMBL::Registry->load_all($$MAIN_CONFIG{LOCATION_OF_COMPARA_REGISTRY_FILE});



# check_analysis_to_configure ( $MAIN_CONFIG, $EXONERATE_CONFIG_BY_LOGIC ) ; 

my $e2g_conf = Bio::EnsEMBL::Analysis::Tools::ConfigUtils->new(
                 -config => "Bio::EnsEMBL::Analysis::Config::Exonerate2Genes");  

my $oa_conf = Bio::EnsEMBL::Analysis::Tools::ConfigUtils->new(
                 -config => "Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluator"); 

my $oa_ex = Bio::EnsEMBL::Analysis::Tools::ConfigUtils->new(
                 -config => "Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluatorExonerate");   

my $dbs = Bio::EnsEMBL::Analysis::Tools::ConfigUtils->new(
                 -config => "Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases"); 

#
# extract the configuration for the default exonerate run out of 
# Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluatorExonerate 
# by it's name in the config   
#

my $basic_xrate_param = $oa_ex->get_config_by_name("EXONERATE_PROTEIN_CONF") ; 


# check which analysis need configuration and configure them  
my %analysis_to_configure = %{$oa_conf->get_config_by_name("MAIN_CONFIG" )};   

my %main_analysis_setup ;  

# configure LOCATE_MISSING_ORTHOLOGUES   
my @initial_analysis_to_run ; 

if ( $analysis_to_configure{"RUN_LOCATE_MISSING_ORTHOLOGUES"}) {  
   push @{$main_analysis_setup{LOCATE_MISSING_ORTHOLOGUES}}, 
    @{ setup_config("LOCATE_MISSING_ORTHOLOGUES",$oa_conf,$basic_xrate_param,$dbs,$e2g_conf)};  
    push @initial_analysis_to_run, 'pre_MissingOrtholouges';
}     

 
# configure RUN_FIND_PARTIAL_GENES  

if ( $analysis_to_configure{"RUN_FIND_PARTIAL_GENES"}) {  
  push @{$main_analysis_setup{FIND_PARTIAL_GENES}}, 
   @{setup_config("FIND_PARTIAL_GENES",$oa_conf,$basic_xrate_param,$dbs,$e2g_conf)}; 
   push @initial_analysis_to_run, 'pre_FindPartialGenes';
}      

if ( $analysis_to_configure{"RUN_FIND_SPLIT_GENES"}) {  
  push @{$main_analysis_setup{FIND_SPLIT_GENES}}, 
   @{setup_config("FIND_SPLIT_GENES",$oa_conf,$basic_xrate_param,$dbs,$e2g_conf)}; 
   push @initial_analysis_to_run, 'pre_FindSplitGenes';
}      



if ($opt{write}){  
   print "writing Exonerate2Genes-configuration and backing up your old one\n"; 
   $e2g_conf->write_config ; 
} else { 
   warning("You haven't used the -write option so i won't write a configuration file for Exonerate2Genes\n");
}

   # ---> post_analysis will be wrapped Exonerate2Genes 
 
my $pa = pipeline_adaptor(\%opt) ;  
 
my $out_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $$DATABASES{"ORTHOLOGUE_DB"} } );    

my @test_runnable_statements ; 

for my $analysis_type ( keys %main_analysis_setup ) {  

    print "Prepare setup for pre_analysis type: $analysis_type\n" ; 
     
      # setup INITIAL analysis   
      my ($pre_ana, $submit ) = get_pre_analysis($analysis_type) ;  
    
      # store in ref-db and out-db
      check_and_store_analysis ($pa, $pre_ana) ; 
      check_and_store_analysis ($pa, $submit) ;   
      check_and_store_analysis ($out_db, $pre_ana) ;  
      add_rule($pa, $pre_ana, $submit) ; 

      # make input_ids for initial analysis 
      # logic_name:slicename for LOCATE_MISSING_ORTHOLOGUES and FIND_PARTIAL_GENES 
      # logic_name for FIND_SPLIT_GENES  
      # "logic_name" refers to the analysis which will run in the second stage 
      #
      
      for my $logic_name ( @{$main_analysis_setup{$analysis_type}} ) {    
        my $dba ;
        my $input_ids ;   

        # get correct core db-adaptor   
        if (  $analysis_type eq "LOCATE_MISSING_ORTHOLOGUES" ) {   

          my $species_alias = ${$$LOCATE_MISSING_ORTHOLOGUES{"ANALYSIS_SETS"}{$logic_name}}[0]; 
          $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species_alias,'core') ; 
          $input_ids = generate_input_ids($dba,$logic_name) ;    

        } elsif (  $analysis_type eq "FIND_PARTIAL_GENES" ) {  
          $dba = $pa ;   
          $input_ids = generate_input_ids($dba,$logic_name) ;      

        } elsif (  $analysis_type eq "FIND_SPLIT_GENES" ) {   
          $dba = $pa ;    
          $input_ids = [ keys %{$$FIND_SPLIT_GENES{"ANALYSIS_SETS"}} ];  
        } 
        push @test_runnable_statements, "perl ensembl-analysis/scripts/test-RunnableDB" .
         " -dbname $opt{dbname} -dbhost $opt{dbhost}\\\n -dbuser $opt{dbuser} " . 
         "-dbpass $opt{dbpass} -dbport $opt{dbport} -analysis $logic_name -input_id $$input_ids[0]\n" ; 
        upload_input_ids ( $input_ids, $pa, $submit) ; 
    }
  

     # setup POST analysis


    for my $logic_name ( @{$main_analysis_setup{$analysis_type}} ) {    

      my ($post_analysis,$submit) = get_analysis_set($analysis_type, $logic_name,$opt{exonerate_file});      
      check_and_store_analysis ( $pa, $post_analysis ) ;  
      check_and_store_analysis ( $pa, $submit ) ;  
      check_and_store_analysis ($out_db, $post_analysis) ; 
      add_rule($pa,$post_analysis,$submit) ;  

    }  
}  

print "\n*** setup done ***\n" ; 
print "You can now run some testRunnables or the rulemanager\n"; 
print "="x80 . "\n\n" ; 
print join ("\n" , @test_runnable_statements) ; 

print "\n\n.... or the rulemanger :\n" . "="x80 . "\n\n" ; 
print "\nperl rulemanager.pl -reread_input_ids -dbhost $opt{dbhost} -dbname $opt{dbname}\\\n -dbport $opt{dbport} -dbpass $opt{dbpass} -dbuser $opt{dbuser} -analysis ".
join(",", @initial_analysis_to_run) . "\n\n" ;  







#
# some subroutines follow .....
#


sub upload_input_ids { 
   my ($input_ids, $dba, $submit_analysis)= @_;  

    my $if = Bio::EnsEMBL::Pipeline::Utils::InputIDFactory->new(
                 -db => $dba , 
                 -logic_name => $submit_analysis->logic_name ,
                 -input_id_type => $submit_analysis->input_id_type,
                 -file=> 1 , 
               ); 

      my $saved_ana = $dba->get_AnalysisAdaptor->fetch_by_logic_name($submit_analysis->logic_name);

      unless ($saved_ana) { 
        throw( "Cant find analysis $submit_analysis in db ") ;
      }

     # check first if input_id is not already stored ...

       my $ia_db = $dba->get_StateInfoContainer->list_input_id_by_Analysis($saved_ana) ;
       my @input_ids_not_stored ;

       my %tmp ;
       @tmp{@$ia_db} = 1;
       for my $i ( @$input_ids ) {
         push @input_ids_not_stored, $i unless (exists $tmp{$i}) ;
       }
      
      print scalar(@$input_ids) - scalar(@input_ids_not_stored) .
      " input ids already stored in db ". $dba->dbname . "\n" ; 

      $if->input_ids(\@input_ids_not_stored) ;
      $if->store_input_ids;
      print scalar(@input_ids_not_stored) . " input-ids uploaded to ".$dba->dbname . "\n" ; 
     
} 


sub generate_input_ids {
   my ($db, $logic_name ) = @_ ; 
   
   my $toplevel_slices = $db->get_SliceAdaptor->fetch_all('toplevel') ;    
   my @input_ids = map { "$logic_name:".$_->name } @$toplevel_slices ;  
   return \@input_ids ; 
}

# calculate input_ids 
#  MissingOrthologues.pm :  logic_name:chromosome:NCBI36:1:3000000:5000000:1 
#   
#   Tried to reply to all but it didn't let me. Thats' why you are getting this e-mail yourself\! 


# this is for later, if we want to have differnt modules for differnt 
# post_analysis ( now it's all exonerate2genes but the plan is to 
# wrap this analysis ...  
#

sub check_and_store_analysis { 
  my ($pa, $analysis) = @_ ;  

  unless ( $pa->get_AnalysisAdaptor->fetch_by_logic_name($analysis->logic_name) ) {   

#      warning("Can't find analysis " . $analysis->logic_name ." in database " . 
#               "\nCreating my very own analysis with hard-coded values ".
#               "out of setup script and rules as well\n") ;  

    print "\nAnalysis ".$analysis->logic_name . " does not exist -storing\n" ; 
    $pa->get_AnalysisAdaptor->store($analysis) ; 
  }
}  


sub get_pre_analysis { 
   my ( $analysis_type ) = @_ ; 

   my ( $sname, $input_id_type ) ; 

   if ( $analysis_type eq "FIND_PARTIAL_GENES") {         

       $input_id_type = "fpg_slice"; 
       $sname = "pre_FindPartialGenes";  

   } elsif ( $analysis_type eq "LOCATE_MISSING_ORTHOLOGUES") {    

       $input_id_type = "mo_slice";  
       $sname = "pre_MissingOrthologues" ; 

   }elsif ( $analysis_type eq "FIND_SPLIT_GENES" ) {    

       $input_id_type = "fsg_slice";   
       $sname = "pre_FindSplitGenes" ;  

   }else{  
     throw("unknown type - can't setup initial analysis") ; 
   }
 
   my $pre_ana = new Bio::EnsEMBL::Pipeline::Analysis ( 
              -logic_name => $sname , 
              -module      => $sname , 
              -input_id_type => $input_id_type , 
            );      

   my $submit_pre = new Bio::EnsEMBL::Pipeline::Analysis ( 
                  -logic_name => "Submit_$sname" , 
                  -module      => "Dummy" , 
                  -input_id_type => $input_id_type, 
                 );      

   return ( $pre_ana, $submit_pre ) ; 
}

sub get_analysis_set { 
   my ( $analysis_type, $logic_name,$exonerate_file ) = @_ ;  

     my ( $input_id_type ) ;  

     if ( $analysis_type eq "FIND_SPLIT_GENES" ) {   
        $input_id_type = "file_$logic_name" ; 
     } elsif ( $analysis_type eq "FIND_PARTIAL_GENES") {       
        $input_id_type = "file_$logic_name" ; 
     } elsif ( $analysis_type eq "LOCATE_MISSING_ORTHOLOGUES") {  
        $input_id_type = "file_$logic_name" ; 
     } else { 
       throw ("Can't find the analysis type specified\n") ; 
    }    


    my $post_analysis = new Bio::EnsEMBL::Pipeline::Analysis ( 
              -logic_name => $logic_name , 
              -program    => 'Exonerate' , 
              -program_file => $exonerate_file , 
              -module      => 'Exonerate2Genes' ,
              -input_id_type => "file_".$logic_name,
            );      

    my $post_submit_analysis = new Bio::EnsEMBL::Pipeline::Analysis
               ( 
                  -logic_name => "Submit_". $logic_name , 
                  -module      => 'Dummy', 
                  -input_id_type => lc("file_".$logic_name) , 
                )  ; 

    return ( $post_analysis ,$post_submit_analysis) ; 
}

    

# SUBS # 


sub pipeline_adaptor{  
  my ($opt )= @_ ; 

  my $dba = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new
            (
               '-dbname' => $$opt{dbname},
               '-host' => $$opt{dbhost},
               '-user' => $$opt{dbuser},
               '-port' => $$opt{dbport},
               '-pass' => $$opt{dbpass},
             ); 
  return $dba ; 
}

 


sub setup_config { 
  my ($analysis_name,$oa_conf,$basic_xrate_param,$dbs,$e2g_conf ) = @_ ;  

  my %conf=%{$oa_conf->get_config_by_name($analysis_name)};
  my @logic_names_to_use = keys %{$conf{"ANALYSIS_SETS"}} ;

  for my $logic_name ( @logic_names_to_use ) {
    print "\n\nSetting up configuration in Exonerate2Genes.pm for analysis : $logic_name\n" ;

    my %merged_config = %{$basic_xrate_param} ;  

    #
    # get OUTDB parameters out of Bio::EnsEMBL::A*::Config::GeneBuild::Databases
    # and add it to the default OrthologueEvaluatorExonerate Config 
    #

    if ( $dbs->exists_in_config("ORTHOLOGUE_DB")){   
      $dbs->get_config_by_name("ORTHOLOGUE_DB"); 
      my %outdb_parameters = %{$dbs->result} ;  
      $merged_config{OUTDB}=$dbs->result ; 
    }else{ 
      throw("Can't find entry for ORTHOLOGUE_DB in ".$dbs->file."\n") ; 
    }  

    # get SEQUENCE_DUMP_BASE_DIR  out of OrthologueEvaluator and merge it 

    $oa_conf->get_config_by_name("SEQUENCE_DUMP_BASE_DIR");
    $merged_config{QUERYSEQS} = $oa_conf->result;

    # change the e2g config 
    $e2g_conf->append_href(\%merged_config, "EXONERATE_CONFIG_BY_LOGIC" , $logic_name ) ;
   } 
   return \@logic_names_to_use ; 
}


sub get_logic_names {
   my ($h)=@_;
   return keys %{$$h{ANALYSIS_SETS}};
} 


sub analysis_already_in_config_file {
   my ($lg_to_check,$config) = @_ ;  
   if ( exists $$config{$lg_to_check}) { 
        return 1; 
   }
   return 0 ; 
}

sub add_rule { 
  my ( $pa, $ana , $condition )= @_;  
  
  my $ra = $pa->get_RuleAdaptor; 
  unless ( $ra->fetch_by_goal($ana)) {  

    print "\nStoring rule :".$ana->logic_name." [condition: ".$condition->logic_name ." ]\n"; 

    my $rule = Bio::EnsEMBL::Pipeline::Rule->new(); 

    $rule->goalAnalysis($pa->get_AnalysisAdaptor->fetch_by_logic_name($ana->logic_name));
    $rule->add_condition($condition->logic_name) ;  

    $pa->get_RuleAdaptor->store($rule) ; 
    print "*stored*\n";
  }else { 
   print "Not storing rule because it's already stored\n"; 
  }
}




