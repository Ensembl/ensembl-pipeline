#!/usr/local/ensembl/bin/perl  -w

############################################################
# script to run the ClusterMerge Algorithm (Eduardo Eyras)
#
# written by Eduardo Eyras (eae@sanger.ac.uk)
#
# You may distribute this module under the same terms as perl itself
#
############################################################

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge;
use Bio::EnsEMBL::Utils::GFF2Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Exon;

my $gff_format = 1;
my $input;
my $output;

############################################################
# deafult options

my $comp_level        = 5;
my $splice_mismatch   = 10;
my $intron_mismatch   = 0;
my $min_order         = 1;
my $restrict_ext_site = 0;
my $sets;
my $help;

############################################################

&GetOptions( 
	     'gff_format:s'        => \$gff_format,
	     'input:s'             => \$input,
	     'output:s'            => \$output,
	     'comp_level:n'        => \$comp_level,
	     'splice_mismatch'     => \$splice_mismatch,
	     'intron_mismatch'     => \$intron_mismatch,
	     'restrict_ext_site:s' => \$restrict_ext_site, 
	     'sets'                => \$sets,
	     'help'                => \$help,
	     );

if ( $help || !$input || !$output){
  &usage;
  exit(0);
}

  
############################################################

open ( IN, "<$input") || die("could not open input file $input");

my $gff2exon = Bio::EnsEMBL::Utils::GFF2Exon->new( -gff_version => $gff_format );

my @transcripts;
my $trans_tag;
my %trans;

INPUT:
while(<IN>){
  
  ############################################################
  # parse gff
  #print STDERR "gff: $_";
  #chomp;
  
  chomp;
  my %trans_tag;
  my $exon = $gff2exon->exon_from_gff( $_ , \%trans_tag );
  #print STDERR "trans_tag[ exon ] = $trans_tag{$exon}\n";
  if ( $exon && $exon->start && $exon->end ){
    push( @{$trans{$trans_tag{$exon}}}, $exon );
  }
}

foreach my $tag ( keys %trans ){
  my $transcript = Bio::EnsEMBL::Transcript->new();
  $transcript->dbID( $tag );
  foreach my $exon ( @{ $trans{$tag} } ){
    $transcript->add_Exon($exon);
  }
  push( @transcripts, $transcript );
}

#print STDERR "transcripts created\n";
#foreach my $t ( @transcripts ){
#    print STDERR "transcript: $t\n";
#    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($t);
#}


############################################################

print STDERR "running ClusterMerge with parameters:\n";
print STDERR "comp_level        = $comp_level\n";
print STDERR "splice_mismatch   = $splice_mismatch\n";
print STDERR "intron_mismatch   = $intron_mismatch\n"; 
print STDERR "min_order         = $min_order\n";
print STDERR "restrict_ext_site = $restrict_ext_site\n";

############################################################

my $cluster_merge = 
    Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge->new(
							-transcripts                   => \@transcripts,
							-comparison_level              => $comp_level,
							-splice_mismatch               => $splice_mismatch,
							-intron_mismatch               => $intron_mismatch,
							-minimum_order                 => $min_order,
							_restrict_external_splice_site => $restrict_ext_site,
						       );

  


$cluster_merge->run;

open ( OUT, ">$output") || die("could not open input file $output");

############################################################
# can retrieve the non-redundant sets:
if ( $sets ){
  
  # list of listrefs, each one cointaining a list of transcript objects
  my @sets = $cluster_merge->sub_clusters;
  
  my $count = 0;
  foreach my $set ( @sets ){
    $count++;
    print OUT "set $count:\n";
    foreach my $transcript ( @$set ){
      foreach my $exon ( @{$transcript->get_all_Exons} ){
	print OUT $gff2exon->gff_string($exon,$transcript);
      }
    }
  } 
}

############################################################
# or the merged transcripts:
else{
  my @merged_transcripts = $cluster_merge->output;
  foreach my $transcript ( @merged_transcripts ){
    foreach my $exon ( @{$transcript->get_all_Exons} ){
      print OUT $gff2exon->gff_string($exon,$transcript)."\n";
    }
  }
}

close (OUT);



sub usage {
    print STDERR <<EOF
    script to run the ClusterMerge algorithm
    Usage: cluster_merge.pl options
    Where options are:

    -input           : name of the input file in gff format

    -output          : name of the output file

    -gff_format      : 1 and 2 supported. (default = 1)

    -comp_level      : comparison level (default = 5 )
                       1 --> strict: exact exon matching. 
                       Does not use any other parameteres passed in. Example:

                       #####-----#####-----#####
                       #####-----#####-----#####

                       2 --> allow edge exon mismatches. 
                       Uses the parameter 'exon_match' if is defined. Example:
    
                         ###-----#####-----#######
                       #####-----#####-----#####-----#####

                       3 ---> allow internal mismatches. 
                       Uses the parameters 'exon_match' and 'splice_mismatch' if defined. Example:

                       #####-----######----#######
                       #####-----######----####------#####

                       4 ---> allow intron mismatches. 
                       This one can use all the parameters if they have been defined.

                       ################----#######
                       #####-----######----####------#####
  
                       5 ---> loose match. It allows intron mismatches if so desired. There is no limitation on
                       the number of mismatches at the splice-sites. Example

                       #################----#######          OR          #######------####----#######  
                       #####-----######----####------#####                #####-----######----####------#####


     -splice_mismatch: maximum number of non-opverlapping nucleotides allowed in splice sites 
                       ( not used at comp_level = 5 )

     -intron_mismatch: maximum number of non-opverlapping nucleotides allowed in introns (default = 0 )

     _rest_ext_site  : restrict extermal splice site; if we want to restrict how much can exceed
                       an external exon an internal splice site:
    
                                       |--d--|
                       ######-------##########
                      #######-------####-------------#######
    
		       If TRUE 'd' must be <= splice_mismatch, else 'd' can be anything

     -exon_match     : TRUE if we want both transcripts to match 1-to-1 all their exons

     -min_order      : minimum number of transcripts required to be in a cluster to create a merged transcript 
                       ( default = 1 )

     -sets           : outputs the non redundant lists instead of the merged transcript 
                       (switched off by default)

     -help           : outputs this help
EOF
}
