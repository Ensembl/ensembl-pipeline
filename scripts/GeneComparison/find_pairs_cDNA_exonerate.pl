#!/usr/local/ensembl/bin/perl -w

=head1 NAME

  find_pairs_cDNA_exonerate.pl

=head1 SYNOPSIS
 
  /path/to/exonerate/output/> find_pairs_cDNA_exonerate.pl -file ../best_pairs
  
=head1 DESCRIPTION
  
  script to process the exonerate output from the comparison of two sets of cDNA sequences. 
  It assumes the output is in one or more gzip'ed files.
  It should be run in the directory where those files exist.

  It takes all matches and does an 'stable marriage' algorithm to find the best pairs. 
  It prints the results in a file given with the option -file
  
=head1 OPTIONS

  -file ../filename
  where teh output goes. It is convenient to send it to another directory not to confuse the script.

=cut

use strict;
use Bio::EnsEMBL::Pipeline::ESTConf;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor;

#$| = 1;


# contig ids are of the form >AB015355.1.1.43999
# look for instance /data/blastdb/Ensembl/NCBI_28_dusted_masked_golden_contigs.fa

# in mouse they are of the form internal_id ( a number )  stable_id ( e.g. 1.77500001-78000000 )

### parameters ###
my $refdb_host  = $EST_REFDBHOST;
my $refdb_name  = $EST_REFDBNAME;
my $refdb_user  = $EST_REFDBUSER;
my $refdb_path  = $EST_GOLDEN_PATH;

my $count = 0;

my $skip_count = 0;
my $feat_count = 0;
my $bad_count  = 0;

my $file;

&GetOptions( 
	    'file:s'  => \$file,
	   );

unless ($file){
  print STDERR "Usage: $0 -file pairs > & log_file\n";
  print STDERR "Run it on the directory where you have the gzip'ed exonerate output\n";
  exit(0);
}

my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $refdb_host,
					       -user   => $refdb_user,
					       -dbname => $refdb_name,
					      );

$refdb->static_golden_path_type( $refdb_path );

# make sure you've copied the analysisprocess table of the previous EST release
# into the current one, then set the analysis_id to be the same as in that table for exonerate:

open LS, 'ls |' or die "Can't open pipe from ls : $!";

my $current_target_id;
my $current_cdna_id;

# each target and cdna id has a list of candidates:
my %target_list;
my %cdna_list;

# keep a score matrix
my %score_matrix;

CHUNK:
while(<LS>){
  chomp;
  my $file = $_;
  my $newfile;
  next CHUNK if ( $file =~/(\.)+\// );
  
  if ( $file =~ /(\S+).gz/){
    $newfile = $1;
    print STDERR "UNCompressing $file\t";
    system("gunzip $file");
  }
  else{
    $newfile = $file;
  }
  
  print STDERR "PROCessing $newfile\n";
  open FILE, "<$newfile";
  
 LINE:
  while(<FILE>){
    
    # output looks like 
    #ENST00000219596 2836    3164    1145.00 86      -1      AF113127        1001    1323    1
    #ENST00000219596 3222    3447    748.00  82      -1      AF113127        1025    1247    1
    #ENST00000219596 3139    3449    1062.00 82      1       AF113129        2498    2809    1
    #ENST00000219596 2871    3157    981.00  84      1       AF113129        2542    2824    1
    #ENST00000301749 2100    2399    988.00  82      -1      AF113127        1026    1321    1
    #ENST00000301749 2542    2800    845.00  81      -1      AF113127        1027    1283    1
    #ENST00000301749 2093    2422    1119.00 82      1       AF113129        2501    2829    1
    #ENST00000301749 2536    2800    974.00  85      1       AF113129        2542    2805    1
    #ENST00000219478 3068    3364    1070.00 86      1       AF113127        1027    1321    1
    #ENST00000219478 4086    4327    881.00  85      -1      AF113127        1025    1266    1
    #ENST00000219478 4066    4329    1021.00 87      1       AF113129        2544    2809    1
    #ENST00000219478 3068    3329    957.00  87      -1      AF113129        2548    2805    1
    
    # the entries are:
    my ( $target_id, 
	 $target_start, 
	 $target_end, 
	 $score, 
	 $perc_id, 
	 $target_strand, 
	 $cdna_id, 
	 $cdna_start, 
	 $cdna_end, 
	 $cdna_strand) = split();
    
    # store the score:
    unless ( $score_matrix{$target_id}{$cdna_id} ){
      $score_matrix{$target_id}{$cdna_id} = $score;
    }
    if (  $score_matrix{$target_id}{$cdna_id}  && $score_matrix{$target_id}{$cdna_id} < $score ){
      $score_matrix{$target_id}{$cdna_id} = $score;
    }

    # Note: Exonerate puts the beginning of the feature at 0, not at 1!!
    if ( $cdna_start == 0 ){
      $cdna_start = 1;
    }
    if ( $target_start == 0 ){
      $target_start = 1;
    }
    
    # first time in, create current ids:
    unless ( $current_target_id){
      $current_target_id = $target_id;
    }
    unless ( $current_cdna_id ){
      $current_cdna_id = $cdna_id;
    }
    
    # if we change target or cdna, update:
    unless ( $target_id eq $current_target_id ){
      $current_target_id = $target_id;
    }
    unless( $cdna_id eq $current_cdna_id){
      $current_cdna_id = $cdna_id;
    }
    push ( @{ $target_list{$target_id} }, $cdna_id );
    push ( @{ $cdna_list{$cdna_id} }, $target_id );
    
  }  # end of LINE

  close(FILE);
  system("gzip $newfile");

}    # close loop over LS


close(LS);

#### find best reciprocal matches:
print STDERR "finding best reciprocal matches\n";

# condition for solution:
# there is no two elements target_id, cdna_id such that they're not paired-up to each other but
# they have better score with each other ( $score_matrix is higher ) than with their corresponding partners

# don't forget to make sure that the lists are sorted by score (use schwartz transform for this)

my @unmarried_targets = keys( %target_list  );

# keep track of which one is married and the partners:
my %married_target;
my %married_cdna;
my %partner;

# go over every unpaired target
MARRIAGE:
while ( @unmarried_targets ){
  
  my $this_target = shift @unmarried_targets;
  #$married_target{ $this_target } = undef;
  
  # sort the potential partners by score in descending order
  @{ $target_list{ $this_target } } =  
    map  { $_->[1] }
  sort { $b->[0] <=> $a->[0] }
  map  { [ $score_matrix{ $this_target }{$_}, $_] } @{ $target_list{ $this_target } };
						       
  # go over the partners until you get married or run out of partners  
 PARTNER:
  while( @{ $target_list{ $this_target }} && !defined($married_target{$this_target}) ){
    
    #print STDERR "checking partner list for $this_target\n";
    my $potential_cdna_partner = shift( @{ $target_list{ $this_target } } );
    
    #print STDERR "looking at $potential_cdna_partner\n";
    # check whether first is already married
    if ( $married_cdna{ $potential_cdna_partner } && $married_cdna{ $potential_cdna_partner } == 1 ){
      
      # is it married to another target?
      if ( !( $partner{ $potential_cdna_partner } eq $this_target ) ){
      
	# is it a 'worse' marriage?
	if ( $score_matrix{ $partner{ $potential_cdna_partner }}{ $potential_cdna_partner } 
	     < $score_matrix{ $this_target }{ $potential_cdna_partner } ){
	  
	  # put the divorced target back into the pool only it it has moer potential partners
	  if ( @{ $target_list{ $partner{ $potential_cdna_partner } } } ){
	    push ( @unmarried_targets, $partner{ $potential_cdna_partner } );
	  }
	  # divorce the target
	  #print STDERR "divorcing ".$partner{ $potential_cdna_partner }."\n";
	  delete $married_target{ $partner{ $potential_cdna_partner } };
	  delete $partner{ $partner{ $potential_cdna_partner } };
	  delete $partner{ $potential_cdna_partner };

	  # let be happier marriage
	  $married_target{ $this_target } = 1;
	  $married_cdna{ $potential_cdna_partner } = 1;
	  $partner{ $potential_cdna_partner } = $this_target;
	  $partner{ $this_target } = $potential_cdna_partner;
	  #print STDERR "new marriage: $this_target --  $potential_cdna_partner\n";
	  next MARRIAGE;
	}
	else{
	  # look at the next potential cdna partner
	  next PARTNER;
	}
      }
      # hmm, this cdna is married, to whom?
      elsif ( $partner{ $potential_cdna_partner } eq $this_target ) {
	# we have already a couple
	$partner{ $this_target } = $potential_cdna_partner;
	next MARRIAGE;
      }
      elsif ( !defined( $partner{ $potential_cdna_partner } ) ){
	# we have a cheater!
	$married_cdna{ $potential_cdna_partner } = 0;
	next PARTNER;
      }
    }
    else{
      
      # this cdna is still single, let be marriage:
      $married_target{ $this_target } = 1;
      $married_cdna{ $potential_cdna_partner } = 1;
      $partner{ $potential_cdna_partner } = $this_target;
      $partner{ $this_target } = $potential_cdna_partner;
      #print STDERR "setting partner{ $this_target } = $potential_cdna_partner\n";
      next MARRIAGE;
    }
    
  } # end of PARTNER
  
}   # end of MARRIAGE
    

open( OUT,">$file" ) or die("cannot open file $file");

foreach my $target ( keys %married_target ){
  print OUT $target."\t".$partner{ $target }."\t".$score_matrix{ $target }{ $partner{$target} }."\n";
}

# since we take all the target_list keys, we will find at least one partner per target_id.
# this is not necessarily good, since matches may be very bad. 
# Bad matches can be filtered afterwards?

close(OUT);

