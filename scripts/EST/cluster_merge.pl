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
my $internal_splice_overlap = 0;
my $sets;
my $help;

############################################################

&GetOptions( 
	    #'gff_format:s'        => \$gff_format,
	     'input:s'                   => \$input,
	     'output:s'                  => \$output,
	     'comp_level:n'              => \$comp_level,
	     'splice_mismatch:n'         => \$splice_mismatch,
	     'intron_mismatch:n'         => \$intron_mismatch,
	     'internal_splice_overlap:n' => \$internal_splice_overlap,
	     'sets'                      => \$sets,
	     'help'                      => \$help,
	     );

if ( $help || !$input || !$output){
  &usage;
  exit(0);
}

  
############################################################

open ( IN, "<$input") || die("could not open input file $input");

#my $gff2exon = Bio::EnsEMBL::Utils::GFF2Exon->new( -gff_version => $gff_format );

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
  my $exon = &exon_from_gff( $_ , \%trans_tag );
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
print STDERR "comp_level              = $comp_level\n";
print STDERR "splice_mismatch         = $splice_mismatch\n";
print STDERR "intron_mismatch         = $intron_mismatch\n"; 
print STDERR "min_order               = $min_order\n";
print STDERR "internal_splice_overlap = $internal_splice_overlap\n";

############################################################

my $cluster_merge = 
    Bio::EnsEMBL::Pipeline::Runnable::ClusterMerge->new(
							-transcripts                   => \@transcripts,
							-comparison_level              => $comp_level,
							-splice_mismatch               => $splice_mismatch,
							-intron_mismatch               => $intron_mismatch,
							-minimum_order                 => $min_order,
							-internal_splice_overlap       => $internal_splice_overlap,
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
	print OUT &gff_string($exon,$transcript);
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
      print OUT &gff_string($exon,$transcript)."\n";
    }
  }
}

close (OUT);

############################################################

sub exon_from_gff {
  my ($gffstring, $trans_hash) = @_;
    
  #print STDERR "gff_string: $gffstring\n";

  #39896	human_cdna	exon	1030942	1031591	100	+	0	6604
  #39897	human_cdna	exon	1034227	1034285	100	+	0	6604
  #39898	human_cdna	exon	1035280	1035339	100	+	0	6604
  #39899	human_cdna	exon	1035741	1035820	100	+	0	6604
  #39900	human_cdna	exon	1036477	1036629	100	+	0	6604
  #243299	human_cdna	exon	1037161	1037229	100	+	0	40150
  #243301	human_cdna	exon	1043180	1043315	100	+	0	40150
  #243303	human_cdna	exon	1048194	1048340	100	+	0	40150
  #243305	human_cdna	exon	1053724	1053855	100	+	0	40150
  #243307	human_cdna	exon	1054346	1054451	100	+	0	40150
  #243309	human_cdna	exon	1055667	1055833	100	+	0	40150
  #243310	human_cdna	exon	1057048	1057221	100	+	0	40150
  #243311	human_cdna	exon	1058529	1058630	100	+	0	40150
  #243313	human_cdna	exon	1061487	1062080	100	+	0	40150
  my ($seqname, 
      $source, 
      $primary, 
      $start, 
      $end, 
      $score, 
      $strand, 
      $frame, 
      @group) 
    = split(/\s+/, $gffstring);
  
  $frame = 0 unless( $frame =~ /^\d+$/);
    
  my $exon = Bio::EnsEMBL::Exon->new();
  $exon->seqname($seqname);
  #$exon->source_tag($source);
  $exon->start($start);
  $exon->end($end);
  my $phase = ( 3 - $frame )%3;
  $exon->phase($phase);
  $exon->end_phase( ( $exon->phase + $exon->length)%3 );
  if ( $score ){
    $exon->score( $score );
  }
  if ( $strand eq '-' ) { $exon->strand(-1); }
  if ( $strand eq '+' ) { $exon->strand(1); }
  if ( $strand eq '.' ) { $exon->strand(0); }

  ############################################################
  # warning: it parses only the first element of the group
  
  $trans_hash->{ $exon } = $group[0];
  return $exon;
}

############################################################

sub gff_string{
  my ($exon, $transcript) = @_;
  my ($str,$source,$primary_tag,$score,$frame,$name,$strand);
  
  if( $exon->can('score') ) {
    $score = $exon->score();
  }
  $score = '.' unless defined $score;
  
  if( $exon->can('frame') ) {
    $frame = $exon->frame();
  }
  $frame = '.' unless defined $frame;
  
  $strand = $exon->strand();
  if(! $strand) {
    $strand = ".";
  } elsif( $strand == 1 ) {
    $strand = '+';
  } elsif ( $exon->strand == -1 ) {
    $strand = '-';
  }
  
  $name        = $exon->seqname();
  $source      = "merged";
  $primary_tag = "exon";
  
  $str = join("\t",
	      $name,
	      $source,
	      $primary_tag,
	      $exon->start(),
	      $exon->end(),
	      $score,
	      $strand,
	      $frame);
  
  my $tag_str = $transcript->type;
  $str .= "\t$tag_str";
  return $str;
}


############################################################

sub usage {
    print STDERR <<EOF
    script to run the ClusterMerge algorithm
    Usage: cluster_merge.pl options
    Where options are:

    -input           : name of the input file in gff format

    -output          : name of the output file

    -comp_level      : comparison level (default = 3 )
                       1 --> strict: exact exon matching. 
                       Does not use any other parameteres passed in. Example:

                       #####-----#####-----#####
                       #####-----#####-----#####

                       2 --> allow edge exon mismatches. 
                       Uses the parameter 'exon_match' and 'internal_splice_overlap' if defined. 
                       Example:
    
                         ###-----#####-----#######
                       #####-----#####-----#####-----#####

                       3 ---> allow internal mismatches. 
                       Uses the parameters 'exon_match', 'splice_mismatch' and 'internal_splice_overlap' 
                       if defined. Example:

                       #####---########----######
                       #####-----######----####------#####

                       4 ---> allow intron mismatches. 
                       This one can use all the parameters if they have been defined.

                       ################----#######
                       #####-----######----####------#####
  
                       5 ---> loose match. It allows intron mismatches if so desired. There is no limitation on
                       the number of mismatches at the splice-sites. Examples:

                       #################----#######           and      #######------####----#######  
                       #####-----######----####------#####               #####-----######----####------#####
                    
                       would be merged as redundant.

     -splice_mismatch: maximum number of non-opverlapping nucleotides allowed in splice sites 
                       ( not used at comp_level = 5 )

     -intron_mismatch: maximum number of non-opverlapping nucleotides allowed in introns (default = 0 )

     -internal_splice_overlap: (default = 0 ) 
                       number of base pairs (N) we allow an external exon overlap
                       an intron in another transcript:
                                       |--N--|
                       ######-------##########
                      #######-------####-------------#######
    
     -exon_match     : TRUE if we want both transcripts to match 1-to-1 all their exons

     -min_order      : minimum number of transcripts required to be in a cluster to create a merged transcript 
                       ( default = 1 )

     -sets           : outputs the non redundant lists instead of the merged transcript 
                       (switched off by default)

     -help           : outputs this help
EOF
}

