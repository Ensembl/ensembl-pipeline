#!/usr/local/ensembl/bin/perl -w


use strict;

use Bio::SeqIO;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster; 
use Bio::EnsEMBL::Analysis::Runnable::TranscriptCoalescer; 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptCoalescer; 
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor; 
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

$| = 1;

# get a contig with a piece-of/entire  chromosome

my ($outfile, $coord_system, $dbhost, $dbname, $dbpass  ) ; 
#my @biotypes ; 
my $dbuser = 'ensro';
my $dbport = 3306;
my ($name ,$slice_size, $analysis_id, $input_id_type,$logic_name )  ;  

&GetOptions(
            'seq_region_name:s' => \$name,  # use 'all' if you want all 
            'dbname:s'       => \$dbname,
            'dbhost:s'       => \$dbhost,
            'dbpass:s' => \$dbpass,
            'dbuser:s' => \$dbuser,
            'dbport:s' => \$dbport,
            'outfile:s'      => \$outfile,
            'coord_system:s' => \$coord_system,
            #'biotypes=s' => \@biotypes ,
            'slice_size=i' => \$slice_size ,  
            'analysis_id=i' => \$analysis_id , 
            'input_id_type=s' => \$input_id_type , 
            'logic_name=s' => \$logic_name, 
	    );

$slice_size = 100_000 unless $slice_size ; 
#@biotypes = split (/,/,join(',',@biotypes)) ; 

my @all_biotypes ;  # array of all biotypes to cluster 

#unless ($input_id_type) { 
# print STDERR " YOU have to supply -input_id_type 1MSLICE or whatever type your analysis has\n" ; 
# exit(0) ; 
#}
#unless ($analysis_id) { 
# print STDERR " YOU have to supply -analysis_id 60  whatever the analysis_id of the Submit-analysis uses\n" ; 
# exit(0) ; 
#} 

unless ($coord_system) {
  print STDERR "you haven't supplied a coord_system_name using   -coord_system --> I'll use 'chromosome' now\n" ; 
  $coord_system = 'chromosome' ; 
}
unless ($outfile) { 
  print STDERR "you haven't supplied a name for output-fasta-file with -outfile--> I'll use 'outfile.txt' now\n" ; 
  $outfile= "outfile_" . $name . ".txt" ;
}
unless ($name) { 
  print "Use option -seq_region_name to provide seq_region_name or use \'all\'\n" ;  
  $name = 'all' ; 
}
unless( $dbhost && $dbname &&  $outfile ){
  print STDERR ("Usage: -dbname $dbname -dbhost $dbhost -dbuser $dbuser ".
                "-coord_system $coord_system  -outfile $outfile\n");
  exit(0);
}





my $db= new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                           -host  => $dbhost,
                                           -user  => $dbuser,
                                           -port => $dbport,
                                           -dbname=> $dbname,
                                           -pass => $dbpass,
                                          );


my $dbp= new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(
                                           -host  => $dbhost,
                                           -user  => $dbuser,
                                           -port => $dbport,
                                           -dbname=> $dbname,
                                           -pass => $dbpass,
                                          );





# new code 
my $analysis;
if (!$analysis_id) {
  $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name) ; 
  $analysis_id = $analysis->dbID() ;   
} else {
  $analysis = $db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
}

if (!$input_id_type) {
  $input_id_type = $dbp->get_AnalysisAdaptor->fetch_analysis_input_id_type($analysis) ; 
}
my @slices; 
my @input_ids ; 


if ($name =~ m/all/i) {
   @slices =  @{ $db->get_SliceAdaptor->fetch_all($coord_system) } ; 
   print STDERR "Having " . scalar ( @slices  ) . " regions to process\n" ;  
   print STDERR "sequences of all slices of type $coord_system will be used ( $outfile )\n" ; 
} else {
  my $slice = $db->get_SliceAdaptor->fetch_by_region($coord_system, $name );
  push @slices, $slice ; 
}


my %database_hash;
my %coalescer_hash = %{$TRANSCRIPT_COALESCER_DB_CONFIG};
my %databases = %{$DATABASES};
for my $category (keys %databases ) {
  
  if(!$databases{$category}{-host} || !$databases{$category}{-dbname}){
    print "Skipping category ".$category." no dbinfo ".
      "defined\n";
    next;
  }
  print "\n$category: $databases{$category}{-host} $databases{$category}{-dbname} :\n--------------------------------------\n" ;
  my %constructor_args = %{$databases{$category}};
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     %constructor_args,
    );
  $database_hash{$category} = $dba;
}

foreach my $slice ( @slices) {
  my @all_genes ; 

  my $coord_sys_name  =  $slice->coord_system_name ; 
  my $coord_sys_version = $slice->coord_system->version ;
  my $seq_region_name = $slice->seq_region_name ; 
 
  # dont construct anything if we have a slice shorter than the slice_size 
  if ($slice->length < $slice_size ) { 
    push @input_ids , $slice->name ; 
    next ; 
  } 

  #
  # getting genes out of different databases defiend by config file 
  #
 
  
   for my $category (keys %coalescer_hash ) {
      my $dba = $database_hash{$category};
     # use slice out of different db   
     my $tmp_slice = $dba->get_SliceAdaptor->fetch_by_name($slice->name) ;
     

     # getting simgw and est genes 
     my $genes ; 
     for my $biotype ( @{$coalescer_hash{$category}{BIOTYPES} }) {
       push @all_biotypes, $biotype ; 
       print "Fetching all genes from ".$tmp_slice->name." by ".
         "biotype ".$biotype."from ".$tmp_slice->adaptor->dbc->dbname."\n";
       my @genes =@{ $tmp_slice->get_all_Genes_by_type($biotype)} ; 
       print "Have " . scalar(@genes) . " genes [$biotype]\n" ; 
       push @all_genes , @genes ; 
     }
     
     # PREDICTION  TRANSCRIPTS 
     for my $logic_name_becomes_biotype ( @{$coalescer_hash{$category}{AB_INITIO_LOGICNAMES} }) { 
       
       push @all_biotypes, $logic_name_becomes_biotype ; 
       
       my $pt = $tmp_slice->get_all_PredictionTranscripts( $logic_name_becomes_biotype ) ;
       print "Have " . scalar(@$pt) . " genes [$logic_name_becomes_biotype]\n" ; 
       my $ab_initio_genes = convert_prediction_transcripts_to_genes(
                                                                     $pt,$logic_name_becomes_biotype ) ;
       push @all_genes, @$ab_initio_genes ; 
     }
   }   
  throw("Slice ".$slice->name." has no genes ") if(@all_genes == 0);
    print "have " .scalar(@all_genes) . " genes for slice ".$slice->name." \n" ; 
  print "having slice_size $slice_size\n" ; 
  
  my $id = create_input_id(\@all_genes, \@all_biotypes, $coord_sys_name, $coord_sys_version, $seq_region_name);
  my @values = split /\:/, $id;
  if(!$values[0] || !$values[2]){
    throw("created id ".$id." isn't right");
  }
  push @input_ids , $id if($id);   
}




open(O,">$outfile") || die "cant write to file $outfile\n" ;
for my $id ( @input_ids ) { 
  print O "insert into input_id_analysis (input_id, input_id_type, analysis_id) values ( \"$id\", \"$input_id_type\", $analysis_id ); \n" ; 
}
close(O) ; 
print "sql-statements written to $outfile - finished\n" ; 





sub create_input_id{
  my ($genes, $biotypes, $coord_sys_name, $coord_sys_version, $seq_region_name) = @_;
  my @all_genes = @$genes;
  my @all_biotypes = @$biotypes;
   @all_genes = sort { $a->seq_region_start <=> $b->seq_region_start }  @all_genes ; 

  my  @cluster  =  @{cluster_Genes(\@all_genes, \@all_biotypes)}  ;  
  @cluster = sort {$a->start <=> $b->start } @cluster ;  
  print "Printing cluster info\n";
   for my $cluster (@cluster) { 
     print $cluster->start . "---" . $cluster->end . "\n" ; 
   }


  my $start = 1 ; 

  for (my $i=0 ; $i<@cluster ; $i++) { 
    my $c = $cluster[$i] ;  
    my $diff = $c->end - $start ; 
    
    if ($diff > $slice_size ) { 
      my $end = $cluster[$i-1]->end ; 
      
      my $id = $coord_sys_name . ":" . $coord_sys_version .":".$seq_region_name .  ":" .  $start . ":" . $end . ":1" ; 
      #print "$id\n" ;
      push @input_ids , $id ;  
      $start = $cluster[$i]->start ; 
      next ;  
    }
    # recover last one 
  }
  my $end = $cluster[-1]->end if(@cluster > 0);  
  my $id = $coord_sys_name . ":" . $coord_sys_version . ":$seq_region_name:"   . $start . ":" . $end . ":1" if(@cluster > 0); 
  return $id;
}



sub cluster_Genes {
  my ($genes, $all_biotypes ) = @_ ; 
  my %tmp ; 
  $tmp{'genes'} = $all_biotypes ; 
  my $types_hash = \%tmp ; 
  #
  # steves old cluster-routine clusters genes of two types : 'ncbi' and 'hinxton' 
  # ( see get_twoay_cluster.pl) 
  # he uses  two gene-sets : genes and compare_genes (each set may contain differnt biotypes)
  #   
  # he uses the sets to see if a cluster only consists of genes out of one set (ncbi) or hinxton 
  # and retrieves all sets of a cluster with "get_sets_included" 
  #
  # $types_hash{'est'} = [est_100, est_200] ... 
  #
  # we do something 'nearly similar : we are clustering genes of diffrent sets (simgw, est, abinitio) 
  # and have methods to access these sets 
  # --> GeneCluster has methods get_Genes_of_Type / get_Genes_by_Type / get_Genes_by_Set
  # all genes on slice are handed over and a %types_hash which holds the setname and the  
  
  return ([],[]) if (!scalar(@$genes));

  # sorting of ALL genes
  my @sorted_genes = 
          sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end }  @$genes;
  
  print "Clustering ".scalar( @sorted_genes )." genes on slice\n" ;

  my @clusters;
  GENE: foreach my $gene (@sorted_genes) {
  
    my @matching_clusters;
    
    ## 
    ## if there are Clusters (initialisation below) than check  
    ## if gene lies in the boundaries of the cluster and has at least 
    ## one exon which overlaps with an exon of a gene which already 
    ## belongs to the cluster
    ##
            
    CLUSTER: foreach my $cluster (@clusters) {
    
    # 
    # if gene lies in the boundaries of the cluster......
    #
    
      if ($gene->end  >= $cluster->start && $gene->start <= $cluster->end) {
    
        # search for a gene in the cluster which overlaps the new gene, 
      
        foreach my $cluster_gene ($cluster->get_Genes){
      
        # check if clustered gene overlaps 
      
          if ($gene->end  >= $cluster_gene->start && $gene->start <= $cluster_gene->end) {
      
            #                             CASE 1: 
            #
            #         START----------------$cluster_gene---------------END 
            # START-------------$gene-----------------------END
            #
            #                             CASE 2 : 
            #                              
            #         START----------------$cluster_gene---------------END 
            #               START-------------$gene-----------END
            #
            #                             CASE 3 : 
            #
            #         START----------------$cluster_gene----------END 
            #                                              START------$gene-------END
            #
            # add gene target-gene to cluster if it has at least
            # one gene wich overlaps with an exon of the clustered gene 
            # and add to cluster  
            #
                  
            if (_compare_Genes( $gene, $cluster_gene)) {
              push (@matching_clusters, $cluster);
              next CLUSTER;
            }
          }
        }
      }
    } # CLUSTER 
        
    ##
    ## Initialization of we have no matching cluster (above) 
    ###############################################################
    
    #
    # if above was found NO matching cluster
    # than make a new one 
    # 
    
    if (scalar(@matching_clusters) == 0) {
      my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster->new();

      foreach my $set_name (keys %$types_hash) {
        $newcluster->gene_Types($set_name,$types_hash->{$set_name});
      }
      
      #$newcluster->gene_Types() ; 

      $newcluster->put_Genes($gene);
      push(@clusters,$newcluster);
    
      #
      # if above was found ONE matching cluster
      #
    } elsif (scalar(@matching_clusters) == 1) {
      $matching_clusters[0]->put_Genes($gene);
  
    } else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster->new();
  
      foreach my $set_name (keys %$types_hash) {
        $merged_cluster->gene_Types($set_name,$types_hash->{$set_name});
      }
      
      my %match_cluster_hash;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_Genes($clust->get_Genes);
        $match_cluster_hash{$clust} = $clust;
      }
      $merged_cluster->put_Genes($gene);
      push @new_clusters,$merged_cluster;
      
      # Add back non matching clusters
      foreach my $clust (@clusters) {
        if (!exists($match_cluster_hash{$clust})) {
          push @new_clusters,$clust;
        }
      }
      @clusters =  @new_clusters;
    }
  }
  
  # Seperate genes which are UNclustered (only one gene in cluster ) and
  # from clusteres which hold more than one gene 

  my (@new_clusters, @unclustered);
  foreach my $cl (@clusters){
    if ( $cl->get_Gene_Count == 1 ){
      push @unclustered, $cl;
    } else{
      push( @new_clusters, $cl );
    }
  }
  print "All Genes clustered\nGot " . scalar(@new_clusters) . " new Clusters\n"  ; 
  push @new_clusters, @unclustered ; 
  #return (\@new_clusters, \@unclustered);
  return \@new_clusters ; 
}



sub _compare_Genes {         
  my ($gene1,$gene2,$translate) = @_;
  # quit if genes do not have genomic overlap 
  #
  # start-------gene1------end   start--------gene2----------end
  #  
  
  if ($gene1->end < $gene2->start || $gene1->start > $gene2->end) {
    print "Gene 1  " . $gene1->start . " " . $gene1->end . " \n";
    print "Gene 2  " . $gene1->start . " " . $gene1->end . " \n";
    print "Failed extents check - returning 0\n";
    return 0;
  }
  
  
  # $overlaps = ( $exon1->end >= $exon2->start && $exon1->start <= $exon2-> end );  
  
  if ($translate) {
    # exon-overlap only on coding exons !
    my $exons1 = get_coding_exons_for_gene($gene1);
    my $exons2 = get_coding_exons_for_gene($gene2);
    foreach my $exon1 (@$exons1) {
      foreach my $exon2 (@$exons2) {
        if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
          #print "Passed CDS overlap check - returning 1\n";
          return 1;
        }
      }
    }
  } else {
    #
    # overlap check based on all (noncoding + coding) Exons 
    #
    foreach my $exon1 (@{$gene1->get_all_Exons}){
      foreach my $exon2 (@{$gene2->get_all_Exons}){
        if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
          #print "Passed exon overlap check (noncod. + cod. exons checked)  - returning 1\n";
          return 1;
        }
      }
    }
  } 
   #print "Failed overlap check (translate = $translate) - returning 0\n";
  return 0;
}      


#


#
sub  convert_prediction_transcripts_to_genes {
  my ($pt,$logic_name_becomes_biotype ) = @_ ;
  my @new_genes ;
  for my $pt (@$pt) {
    # conversion 
    my $gene_from_pt = Bio::EnsEMBL::Gene->new(
                       -start => $pt->start ,
                       -end => $pt->end ,
                       -strand => $pt->strand ,
                       -slice =>$pt->slice ,
                       -biotype => $logic_name_becomes_biotype,
                       -analysis=>$pt->analysis,
                       ) ;

    my $new_tr = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended->new(
                    -BIOTYPE => $logic_name_becomes_biotype ,
                    -ANALYSIS => $pt->analysis ,
                 ) ;

    my @pt_exons  = @{$pt->get_all_Exons} ;

    for (my $i=0 ; $i<scalar(@pt_exons) ; $i++) {

      # converting Bio::EnsEMBL::PredictionExon into ExonExtened (ISA Bio::EnsEMBL::Exon)  
      my $pte =$pt_exons[$i] ;
      bless $pte,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended" ;
      $pte->end_phase(0);
      $pte->phase(0);
      $pte->next_exon($pt_exons[$i+1]) ;
      $pte->prev_exon($pt_exons[$i-1]) ;
      $pte->transcript($new_tr) ;
      $pte->analysis($pt->analysis) ;
    } ;

    # my $new_tr = Bio::EnsEMBL::Transcript->new( -exons => $pt_exons  ) ;  
    # bless $new_tr, "Bio::EnsEMBL::Analysis::Runnable::Condense_EST::TranscriptExtended" ; 
    #
    # Extending the Bio::EnsEMBL::Transcript object by ev_set, ev_rank methods 
    #
    for (@pt_exons) {
      $new_tr->add_Exon($_);
    }

    $gene_from_pt->add_Transcript($new_tr) ;

    push @new_genes , $gene_from_pt ;
  }
  return \@new_genes ;
}





