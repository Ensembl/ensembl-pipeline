#!/usr/local/ensembl/bin/perl -w

=head1 NAME

make_input_ids_for_similarity_build 

=head1 SYNOPSIS

make_input_ids_for_similarity_build.pl

=head1 DESCRIPTION


=head1 OPTIONS

=head1 EXAMPLES

=cut

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(parse_config);
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BlastMiniGenewise qw(GENEWISE_CONFIG_BY_LOGIC);
use strict;
use Getopt::Long;
use SimilarityInputIdConfig;

my $logic_name;
my $write;
my $help;
my $max_slice_size;
my $verbose;
my $coord_system = 'toplevel';
my $slice_name_file;
my $protein_count = 20;
my @restricted_prots;
my $pipeline_dbname = 'REFERENCE_DB';
my $output_logic_name;

&GetOptions(
            'output_logic_name=s' => \$output_logic_name,
            'bmg_logic_name=s'      => \$logic_name,
            'write!'            => \$write,
            'help'              => \$help,
            'max_slice=s'       => \$max_slice_size,
            'verbose'           => \$verbose,
            'coord_system=s'    => \$coord_system,
            'slice_name_file:s' => \$slice_name_file,
            'number_of_proteins_per_job:i' => \$protein_count , 
            'pipeline_dbname:s' => \$pipeline_dbname,
           );

exec( 'perldoc', $0 ) if $help;


my $config_object = SimilarityInputIdConfig->new();

parse_config($config_object, $GENEWISE_CONFIG_BY_LOGIC, $logic_name);

if(!$protein_count){
  $protein_count = 20;
  print "Using default protein number per job 20\n";
}

my $pipeline_db = create_db($pipeline_dbname, 1);
my $aa = $pipeline_db->get_AnalysisAdaptor;

my $analysis = $aa->fetch_by_logic_name($logic_name);
if ( not $analysis ) {
  throw "Could not find analysis for $logic_name\n";
}

my ( $ana_id, $ana_type ) = ( $analysis->dbID, $analysis->input_id_type );
if ( not $ana_type ) {
  throw "Could not find dbID/input_id_type for $logic_name\n";
}

my $output_analysis = $aa->fetch_by_logic_name($output_logic_name);
throw("Failed to fetch an analysis for ".$output_logic_name." from ".$pipeline_db->dbc->dbname) if(!$output_analysis);
throw("Can't store input ids on ".$output_analysis." for ".$analysis->logic_name.
      " because input id types are mismatched") 
  if($output_analysis->input_id_type ne $analysis->input_id_type);

if(!$max_slice_size){
  $max_slice_size = 1000000;
  print "Using default max slice size of 1000000\n";
}

print "Generating initial input ids\n" if($verbose);


my @slice_names;
if ($slice_name_file) {
  open( SLICES, "<$slice_name_file" ) or throw "could not open file with slice $slice_name_file";
  while (<SLICES>) {
    /^(\S+):(\S*):(\S+):(\d*):(\d*):/ and do {
      my $inputIDFactory = Bio::EnsEMBL::Pipeline::Utils::InputIDFactory->
          new(
              -db                   => $pipeline_db,
              -slice                => 1,
              -slice_size           => $max_slice_size,
              -coord_system         => $1,
              -coord_system_version => $2,
              -seq_region_name      => $3,
              -seq_region_start     => $4, 
              -seq_region_end       => $5,
              -logic_name    => $logic_name
              );
      
      push @slice_names, @{$inputIDFactory->generate_input_ids};
    }
  }
} else {
  my $inputIDFactory = Bio::EnsEMBL::Pipeline::Utils::InputIDFactory->
      new(
          -db           => $pipeline_db,
          -slice        => 1,
          -slice_size   => $max_slice_size,
          -coord_system => $coord_system,
          -logic_name   => $logic_name
          );
  @slice_names = @{ $inputIDFactory->generate_input_ids }
}



my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => 'protein');
my $kill_list = $kill_list_object->get_kill_list;

my @input_ids;

foreach my $slice_name(@slice_names){
  my @iids = get_input_ids_from_slice($slice_name, $config_object, $protein_count, $kill_list);
  print "Have ".@iids." from ".$slice_name."\n";
  foreach my $input_id(@iids){
    if(!$input_id){
      throw("Seem to have an undefined input_id");
      print $input_id."\n";
    }
  }
  push(@input_ids, @iids);
  print "Have ".@iids." from ".$slice_name."\n" if($verbose);
}

foreach my $input_id(@input_ids){
  my $sic = $pipeline_db->get_StateInfoContainer;
  if($write){
    eval{
      print "Storing ".$input_id." with ".$output_analysis->logic_name."\n" if($verbose);
      $sic->store_input_id_analysis($input_id, $output_analysis, '');
    };
    if($@){
      print "Problem with ".$input_id." store $@\n";
    }
  }
}

sub create_db{
  my ($string, $pipeline) = @_;
  
  my $constructor_args = $DATABASES->{$string}; 
  
  foreach my $arg ( qw ( -user -port -host -dbname) ) {  
    unless ( $$constructor_args{$arg}){ 
      throw ("Database-connection-details not properly configured : Arguemnt : $arg missing in Databases.pm\n") ; 
    }
  }
  
  my $db;
  if($pipeline){
    $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
                                                        %$constructor_args,
                                                       );
  }else{
    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                                              %$constructor_args,
                                             );
  }
  return $db;
}


sub get_input_ids_from_slice{
  my ($name, $config, $protein_count, $kill_list) = @_;

  my $paf_db = create_db($config->PAF_SOURCE_DB);
  my $gw_db = create_db($config->GENE_SOURCE_DB);

  my $paf_slice = $paf_db->get_SliceAdaptor->fetch_by_name($name);
  my $gw_slice = $gw_db->get_SliceAdaptor->fetch_by_name($name);
  my @iids;
  my @mask_exons;
  foreach my $type(@{$config->BIOTYPES_TO_MASK}){
    print "Getting genes of type ".$type." from ".$gw_slice->name."\n" if($verbose);
    foreach my $gene(@{$gw_slice->get_all_Genes_by_type($type)}){
      foreach my $mask_exon(@{$gene->get_all_Exons}){
        push(@mask_exons, $mask_exon);
      }
    }
  }
  my @mask_regions;
  foreach my $mask_exon ( sort { $a->start <=> $b->start } @mask_exons ) {
    if ( @mask_regions and $mask_regions[-1]->{'end'} > $mask_exon->start ) {
      if ( $mask_exon->end > $mask_regions[-1]->{'end'} ) {
        $mask_regions[-1]->{'end'} = $mask_exon->end;
      }
    } else {
      push @mask_regions, { start => $mask_exon->start, end => $mask_exon->end }
    }
  }
  my $num_seeds = 0;
   foreach my $logicname(@{$config->PAF_LOGICNAMES}) {
    my %features;
    my @features = @{$paf_slice->get_all_ProteinAlignFeatures($logicname)};
      FEATURE:foreach my $f(@features){
        next FEATURE if($config->PAF_MIN_SCORE_THRESHOLD && $f->score < $config->PAF_MIN_SCORE_THRESHOLD);
        next FEATURE if($config->PAF_UPPER_SCORE_THRESHOLD && $f->score > $config->PAF_UPPER_SCORE_THRESHOLD);
        push(@{$features{$f->hseqname}}, $f);
    }
    
    my @ids_to_ignore;
  SEQID: foreach my $sid ( keys %features ) {
      my $ex_idx = 0;
      my $count  = 0;
      
      #print STDERR "Looking at $sid\n";
    FEAT: foreach my $f ( sort { $a->start <=> $b->start } @{ $features{$sid} } ) {
        
        #printf STDERR "Feature: %d %d\n", $f->start, $f->end;
        for ( ; $ex_idx < @mask_regions ; ) {
          my $mask_exon = $mask_regions[$ex_idx];
          
          #printf STDERR " Mask exon %d %d\n", $mask_exon->{'start'}, $mask_exon->{'end'};
          if ( $mask_exon->{'start'} > $f->end ) {
            
            # no exons will overlap this feature
            next FEAT;
          } elsif ( $mask_exon->{'end'} >= $f->start ) {
            
            # overlap
            push @ids_to_ignore, $f->hseqname;
            
            #printf STDERR "Ignoring %s\n", $f->hseqname;
            next SEQID;
          } else {
            $ex_idx++;
          }
        }
      }
    }
    print "Ignoring ".@ids_to_ignore." features\n";
    foreach my $dud_id ( @ids_to_ignore, keys %{$kill_list} ) {
      if ( exists $features{$dud_id} ) {
        delete $features{$dud_id};
      }
    }
    
    $num_seeds += scalar( keys %features );
  }
  # rule of thumb; split data so that each job constitutes one piece of
  # genomic DNA against ~20 proteins.
  #
  return () if($num_seeds == 0);

  my $num_chunks = int( $num_seeds / $protein_count ) 
    + 1;
  for ( my $x = 1 ; $x <= $num_chunks ; $x++ ) {
    
    #generate input id : $chr_name.1-$chr_length:$num_chunks:$x
    my $new_iid = $paf_slice->name . ":$num_chunks:$x";
    push @iids, $new_iid;
  }
  return @iids;
}


package SimilarityInputIdConfig;

sub PAF_LOGICNAMES{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_LOGICNAMES} = $arg;
  }
  return $self->{PAF_LOGICNAMES}
}

sub PAF_MIN_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_MIN_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_MIN_SCORE_THRESHOLD}
}

sub PAF_UPPER_SCORE_THRESHOLD{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_UPPER_SCORE_THRESHOLD} = $arg;
  }
  return $self->{PAF_UPPER_SCORE_THRESHOLD}
}

sub PAF_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PAF_SOURCE_DB} = $arg;
  }
  return $self->{PAF_SOURCE_DB}
}

sub GENE_SOURCE_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_SOURCE_DB} = $arg;
  }
  return $self->{GENE_SOURCE_DB}
}

sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_DB} = $arg;
  }
  return $self->{OUTPUT_DB}
}

sub OUTPUT_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{OUTPUT_BIOTYPE} = $arg;
  }
  return $self->{OUTPUT_BIOTYPE}
}

sub GENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENEWISE_PARAMETERS} = $arg;
  }
  return $self->{GENEWISE_PARAMETERS}
}

sub MINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MINIGENEWISE_PARAMETERS}
}

sub MULTIMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{MULTIMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{MULTIMINIGENEWISE_PARAMETERS}
}

sub BLASTMINIGENEWISE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BLASTMINIGENEWISE_PARAMETERS} = $arg;
  }
  return $self->{BLASTMINIGENEWISE_PARAMETERS}
}


sub EXONERATE_PARAMETERS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXONERATE_PARAMETERS} = $arg;
  }
  return $self->{EXONERATE_PARAMETERS}
}


sub FILTER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_PARAMETERS} = $arg;
  }
  return $self->{FILTER_PARAMETERS}
}

sub FILTER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FILTER_OBJECT} = $arg;
  }
  return $self->{FILTER_OBJECT}
}

sub BIOTYPES_TO_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{BIOTYPES_TO_MASK} = $arg;
  }
  return $self->{BIOTYPES_TO_MASK}
}

sub EXON_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{EXON_BASED_MASKING} = $arg;
  }
  return $self->{EXON_BASED_MASKING}
}

sub GENE_BASED_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{GENE_BASED_MASKING} = $arg;
  }
  return $self->{GENE_BASED_MASKING}
}


sub POST_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{POST_GENEWISE_MASK} = $arg;
  }
  return $self->{POST_GENEWISE_MASK}
}

sub PRE_GENEWISE_MASK{
  my ($self, $arg) = @_;
  if($arg){
    $self->{PRE_GENEWISE_MASK} = $arg;
  }
  return $self->{PRE_GENEWISE_MASK}
}

sub REPEATMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REPEATMASKING} = $arg;
  }
  return $self->{REPEATMASKING}
}

sub SEQFETCHER_OBJECT{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_OBJECT} = $arg;
  }
  return $self->{SEQFETCHER_OBJECT}
}

sub SEQFETCHER_PARAMS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SEQFETCHER_PARAMS} = $arg;
  }
  return $self->{SEQFETCHER_PARAMS}
}

sub USE_KILL_LIST{
  my ($self, $arg) = @_;
  if($arg){
    $self->{USE_KILL_LIST} = $arg;
  }
  return $self->{USE_KILL_LIST}
}

sub LIMIT_TO_FEATURE_RANGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{LIMIT_TO_FEATURE_RANGE} = $arg;
  }
  return $self->{LIMIT_TO_FEATURE_RANGE}
}


sub FEATURE_RANGE_PADDING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{FEATURE_RANGE_PADDING} = $arg;
  }
  return $self->{FEATURE_RANGE_PADDING}
}

sub WRITE_REJECTED{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{WRITE_REJECTED} = $arg;
  }
  return $self->{WRITE_REJECTED};
}

sub REJECTED_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{REJECTED_BIOTYPE} = $arg;
  }
  return $self->{REJECTED_BIOTYPE};
}

sub SOFTMASKING{
  my ($self, $arg) = @_;
  if($arg){
    $self->{SOFTMASKING} = $arg;
  }
  return $self->{SOFTMASKING}
}

1;
