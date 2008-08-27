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

   The setup script requires some additonal libs :  

       Data::Dumper FindBin  PPI  Clone  List::MoreUtils 

   If you like to install these libs be aware that MoreUtils v. 1.18 is not compatible with 
   Utils v. 1.14. 

   Sanger-internals can add these modules (compiled for linux) to their PERL5LIB by 
   sourcing this file :

                     /nfs/acari/jhv/lib/set_orth_path.sh  
  
   Run setup script to set up the analysis and the input_ids in the reference database and 
   to store rules and conditions. This script also re-writes the Exonerte2Genes config 
   in Bio/EnsEMBL/Analysis/Config and backs up your old config - you need to run the setup
   script with the -write option to rewrite the Exonerate2Genes.pm config ( this is for the 
   Exonerate runs in the second stage of the analysis). Your existing Exonerate2Genes-file 
   file will backed up to Exonerate2Genes.pm.bak.0 and there will be a timestamp added 
   in the top when the file has been backed up - Sorry to say but all the comments in your 
   original file can't be transferred to the new file. 


   You find the whole documentation for the OrthologueEvaluation in 

        ensembl-doc/pipeline_docs/orthologue_evaluator_setup.txt 



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
use Bio::EnsEMBL::Analysis::Config::Databases;
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
# The names of the initial analysis are hard-coded in subroutine : get_pre_analysis.   
# I know that hard_coding stuff is not really ideal ... feel free to change this if you
# fancy other logic_names ...
#  
# The are :  
#
#  FindMissingOrtholouges
#  FindPartialGenes
#  FindSplitGenes 
#  FindFalseIntrons

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
                 -config => "Bio::EnsEMBL::Analysis::Config::Databases"); 

#
# extract the configuration for the default exonerate run out of 
# Bio::EnsEMBL::Analysis::Config::GeneBuild::OrthologueEvaluatorExonerate 
# by it's name in the config   
#

my $basic_xrate_param = $oa_ex->get_config_by_name("EXONERATE_PROTEIN_CONF") ; 


# check which analysis need configuration and configure them  
my %analysis_to_configure = %{$oa_conf->get_config_by_name("MAIN_CONFIG" )};   

my %main_analysis_setup ;  

# configure FIND_MISSING_ORTHOLOGUES 
my @initial_analysis_to_run ; 

if ( $analysis_to_configure{"RUN_FIND_MISSING_ORTHOLOGUES"}) {  
 
   my @analysis_sets_logic_names =
     @{ setup_config("FIND_MISSING_ORTHOLOGUES",$oa_conf,$basic_xrate_param,$dbs,$e2g_conf,1)};  

   push @{$main_analysis_setup{FIND_MISSING_ORTHOLOGUES}}, @analysis_sets_logic_names;
   push @initial_analysis_to_run, 'pre_FindMissingOrthologues';
}     

 
# configure FIND_PARTIAL_GENES  

if ( $analysis_to_configure{"RUN_FIND_PARTIAL_GENES"}) {   

  push @{$main_analysis_setup{FIND_PARTIAL_GENES}}, 
    @{setup_config("FIND_PARTIAL_GENES",$oa_conf,$basic_xrate_param,$dbs,$e2g_conf,1)}; 

  push @initial_analysis_to_run, 'pre_FindPartialGenes';
}      

# configure FIND_SPLIT_GENES 
if ( $analysis_to_configure{"RUN_FIND_SPLIT_GENES"}) {   

  push @{$main_analysis_setup{FIND_SPLIT_GENES}}, 
   @{setup_config("FIND_SPLIT_GENES",$oa_conf,$basic_xrate_param,$dbs,$e2g_conf,1)};  

  push @initial_analysis_to_run, 'pre_FindSplitGenes';
}       

# configure FIND_FALSE_INTRONS
if ( $analysis_to_configure{"RUN_FIND_FALSE_INTRONS"}) {  
  push @{$main_analysis_setup{FIND_FALSE_INTRONS}}, 
   @{setup_config("FIND_FALSE_INTRONS",$oa_conf,$basic_xrate_param,$dbs,$e2g_conf,0)};  

   push @initial_analysis_to_run, 'FindFalseIntrons';
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

      # make and upload input_ids for initial analysis  
      
       # - logic_name:slicename for FIND_MISSING_ORTHOLOGUES 
       # - logic_name:slicename for FIND_PARTIAL_GENES 
       # - slicename for FIND_FALSE_INTRONS
       # - logic_name for FIND_SPLIT_GENES   
      
      # "logic_name" refers to the analysis which will run in the second stage 
       
      
      for my $logic_name ( @{$main_analysis_setup{$analysis_type}} ) {    
        my $dba ;
        my $input_ids ;   

        # get correct core db-adaptor    
        
        if (  $analysis_type eq "FIND_MISSING_ORTHOLOGUES" ) {   

          my $species_alias = ${$$FIND_MISSING_ORTHOLOGUES{"ANALYSIS_SETS"}{$logic_name}}[0]; 
	  $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species_alias,'core') ;  

	  $input_ids = generate_input_ids($dba,$logic_name,1) ;    

        } elsif (  $analysis_type eq "FIND_PARTIAL_GENES" ) {  
          $dba = $pa ;   
          $input_ids = generate_input_ids($dba,$logic_name,1) ;      

        } elsif (  $analysis_type eq "FIND_SPLIT_GENES" ) {      

          
          # FindSplitGenes queries the compara database with raw sql, 
          # it does not require a proper input-id ( slice format or whatever ) 
          # any random string can be used !!!! 
          
          $dba = $pa ;    
          $input_ids = [ keys %{$$FIND_SPLIT_GENES{"ANALYSIS_SETS"}} ];   

        } elsif (  $analysis_type eq "FIND_FALSE_INTRONS" ) {     
          
          # FindFalseIntrons gets all genes out of the QUERY_SPECIES core db on 
          # the toplevels and invstigates the genes of configured biotypes i
          # input_ids are slice-formatted toplevel regions 
          # chromosome:BROADD2:X:12345:1  
          
          $dba = $pa ;    
          $input_ids = generate_input_ids($dba,$logic_name,0) ;      
        }  

        push @test_runnable_statements, "perl ensembl-analysis/scripts/test_RunnableDB" .
         " -dbname $opt{dbname} -dbhost $opt{dbhost}\\\n -dbuser $opt{dbuser} " . 
         " -dbpass $opt{dbpass} -dbport $opt{dbport} -analysis $logic_name ".
         " -input_id $$input_ids[0]\n" ;   

        upload_input_ids ( $input_ids, $pa, $submit) ; 
     }
  

     # setup POST analysis 

    for my $logic_name ( @{$main_analysis_setup{$analysis_type}} ) {    

      my ($post_analysis,$submit)= get_post_analysis_set($analysis_type, $logic_name,$opt{exonerate_file});      
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
      print  "\n" . $submit_analysis->logic_name .  " : uploading input_ids !\n" ; 
     # check first if input_id is not already stored ...

       my $ia_db = $dba->get_StateInfoContainer->list_input_id_by_Analysis($saved_ana) ;
       my @input_ids_not_stored ;

       my %tmp ;
       @tmp{@$ia_db} = 1;
       for my $i ( @$input_ids ) {
         push @input_ids_not_stored, $i unless (exists $tmp{$i}) ;
       }
     
       
      print $submit_analysis->logic_name . " : " ; 
      print  scalar(@$input_ids) - scalar(@input_ids_not_stored) .
        " input ids already stored in db ". $dba->dbname . "\n" ; 

      $if->input_ids(\@input_ids_not_stored) ;
      $if->store_input_ids;
      print $submit_analysis->logic_name . " : " . scalar(@input_ids_not_stored) .
       " input-ids uploaded to ".$dba->dbname . "\n" ; 
     
} 


sub generate_input_ids {
   my ($db, $logic_name,$use_logic_name_in_input_id ) = @_ ; 
   
   my $toplevel_slices = $db->get_SliceAdaptor->fetch_all('toplevel') ;    
   my @clean_slices;
   foreach my $t (@$toplevel_slices){
     if ($t->seq_region_name =~/MT/){
       warning( "Genes on slice ".$t->seq_region_name." will not be used because slice is MT region\n");
     }else{
       push @clean_slices, $t;
     }  
   } 
   my @input_ids; 


   
   if ( $use_logic_name_in_input_id ) {  

      # if we run FIND_PARTIAL_GENES or FIND_MISSING_ORTHOLOGUES there's the option to 
      # configure multiple analysis-sets so the input_id contains the name of the analysis-set  
     
      #     hum_mus_set:chromosome:NCBI36:1:3000000:5000000:1
     
      @input_ids = map { "$logic_name:".$_->name } @clean_slices ;  

   } else {   

     # FIND_FALSE_INTRONS does only have normal slices without the logic-name 
     # as input-id
    
     # chromosome:BROADD2:X:1:234545:1  
     @input_ids = @clean_slices ; 
   } 

   return \@input_ids ; 
}

# calculate input_ids 
#  FindMissingOrthologues.pm :  logic_name:chromosome:NCBI36:1:3000000:5000000:1 
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

    print "\nAnalysis ".$analysis->logic_name . " does not exist ...."; 
    $pa->get_AnalysisAdaptor->store($analysis) ;  
    print " **stored**\n"; 
  }
}  


sub get_pre_analysis { 
   my ( $analysis_type ) = @_ ; 

   my ( $sname, $input_id_type, $module ) ; 

   if ( $analysis_type eq "FIND_PARTIAL_GENES") {         

       $input_id_type = "fpg_slice"; 
       $sname = "pre_FindPartialGenes";  
       $module = "FindPartialGenes";  

   } elsif ( $analysis_type eq "FIND_MISSING_ORTHOLOGUES") {    

       $input_id_type = "mo_slice";  
       $sname = "pre_FindMissingOrthologues" ; 
       $module = "FindMissingOrthologues" ; 

   }elsif ( $analysis_type eq "FIND_SPLIT_GENES" ) {    
       
       $input_id_type = "fsg_slice";   
       $sname = "pre_FindSplitGenes" ; 
       $module = "FindSplitGenes" ;  

   }elsif ( $analysis_type eq "FIND_FALSE_INTRONS" ) {    

       # runs on whole toplevel-regions / chromosomes 
       # chromosome:BROADD1:1:34567:1 
       # pre-analysis  : FindFalseIntrons.pm 
       # post-analysis : RecoverFalseIntrons.pm 
       # ( post analysis is set up by RunnableDB we might want to change this )  
       
       $input_id_type = "ffi_slice";   
       $sname = "FindFalseIntrons" ; 
       $module = "FindFalseIntrons" ;  

   }else{  
     throw("unknown type - can't setup initial analysis") ; 
   }
 
   my $pre_ana = new Bio::EnsEMBL::Pipeline::Analysis ( 
              -logic_name => $sname , 
              -module      => $module, 
              -input_id_type => $input_id_type , 
            );      

   my $submit_pre = new Bio::EnsEMBL::Pipeline::Analysis ( 
                  -logic_name => "Submit_$sname" , 
                  -module      => "Dummy" , 
                  -input_id_type => $input_id_type, 
                 );      

   return ( $pre_ana, $submit_pre ) ; 
}




# configure the post analysis :  

# -  for FIND_MISSING_ORTHOLOGUES / FIND_SPLIT_GENES / FIND_PARTIAL_GENES 
#    an exoneate run with module Exonerate2Genes. the logic_name is the key of the 
#    ANALYSIS_SETS - hash in the OrthologueEvaluator.pm config file. 

# - FIND_FALSE_INTRONS has an RecoverFalseIntrons run as post analysis 


sub get_post_analysis_set { 
   my ( $analysis_type, $logic_name,$exonerate_file ) = @_ ;  

     my ( $input_id_type,$post_analysis,$post_submit_analysis ) ;  


     if ( $analysis_type eq "FIND_SPLIT_GENES" ) {    

        $input_id_type = "file_$logic_name" ;  

     } elsif ( $analysis_type eq "FIND_PARTIAL_GENES") {        

        $input_id_type = "file_$logic_name" ;  

     } elsif ( $analysis_type eq "FIND_MISSING_ORTHOLOGUES") {   

        $input_id_type = "file_$logic_name" ;  

     } elsif ( $analysis_type eq "FIND_FALSE_INTRONS") {   

        # FindFalseIntrons is a bit different to setup 
        $input_id_type = "prot_slice_form" ;  
        $logic_name    = "RecoverFalseIntrons";
 
        $post_analysis = new Bio::EnsEMBL::Pipeline::Analysis 
           ( 
              -logic_name => $logic_name ,
              -module      => "RecoverFalseIntrons",
              -input_id_type => $input_id_type ,
            );      

        $post_submit_analysis = new Bio::EnsEMBL::Pipeline::Analysis
               ( 
                  -logic_name => "Submit_". $logic_name , 
                  -module      => 'Dummy', 
                  -input_id_type => lc($input_id_type), 
                )  ; 

     } else { 
       throw ("Can't find the analysis type specified\n") ; 
     }    

    unless ( $post_analysis ) { 
       # if post_analysis is not already set up for FIND_FALSE_INTRONS above 
       $post_analysis = new Bio::EnsEMBL::Pipeline::Analysis ( 
              -logic_name => $logic_name , 
              -program    => 'Exonerate' , 
              -program_file => $exonerate_file , 
              -module      => 'Exonerate2Genes' ,
              -input_id_type => "file_".$logic_name,
            );      

       $post_submit_analysis = new Bio::EnsEMBL::Pipeline::Analysis
               ( 
                  -logic_name => "Submit_". $logic_name , 
                  -module      => 'Dummy', 
                  -input_id_type => lc("file_".$logic_name) , 
                )  ; 
    }
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

 

# checks each config hash for keys in ANALYSIS_SETS ( sub-analyses to configure ) 
# - sets up config in Exonerate2Genes.pm  

sub setup_config { 
  my ($analysis_name,$oa_conf,$basic_xrate_param,$dbs,$e2g_conf,$write_config ) = @_ ;  

  my %conf=%{$oa_conf->get_config_by_name($analysis_name)}; 

  my @logic_names_to_use = keys %{$conf{"ANALYSIS_SETS"}} ;

  if ( $write_config ) { 
    for my $logic_name ( @logic_names_to_use ) { 
       print "\n\nSetting up configuration in Exonerate2Genes.pm for analysis : $logic_name\n" ; 
       
       my %merged_config = %{$basic_xrate_param} ;  
  
       #
       # get OUTDB parameters out of Bio::EnsEMBL::Analysis::Config::Databases
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
   } else {  
     print "no need to write Exonerate2Genes-config\n";
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

    # check if the rule already exists in the database  
    
    print "\nStoring rule : ".$ana->logic_name." [condition: ".$condition->logic_name ." ]\n"; 

    my $rule = Bio::EnsEMBL::Pipeline::Rule->new(); 

    $rule->goalAnalysis($pa->get_AnalysisAdaptor->fetch_by_logic_name($ana->logic_name)); 
    $rule->add_condition($condition->logic_name) ;   
    $pa->get_RuleAdaptor->store($rule) ;  
  }else { 
   print "Not storing rule because it's already stored\n"; 
  }
}




