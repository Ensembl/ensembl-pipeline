#module for running best in genome analysis for pmatchs docs in ensembl-docs

package Bio::EnsEMBL::Pipeline::RunnableDB::BestPmatch;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General;
use Bio::EnsEMBL::Pipeline::Runnable::Pmatch;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::Second_PMF;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::PmatchFeature;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
   

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);



sub fetch_input{
  my ($self) = @_;

  my @pmatch_features = $self->db->get_PmatchFeatureAdaptor->fetch_by_logic_name($GB_INITIAL_PMATCH_LOGICNAME);

  $self->pmatch_features(@pmatch_features);
}

sub run{
  my ($self) = @_;

  my @hits = $self->pmatch_features;
  my $hits = \@hits;
  my $pmf2 = Bio::EnsEMBL::Pipeline::Tools::Pmatch::Second_PMF->new( 
							     '-phits' => $hits,
							    );
  
  $pmf2->run;

  my @output = $pmf2->output;
  my @unique = $self->uniquify(@output);
  $self->output(@unique);
  @unique = undef;
  
}



sub write_output{
  my ($self) = @_;
  my $pmfa = $self->db->get_PmatchFeatureAdaptor();
  $pmfa->write_PmatchFeatures($self->output);
}



##################
#accessor methods#
##################

sub pmatch_features{
  my ($self, @args) = @_;
  if(!$self->{'_pmatch_features'}){
    $self->{'_pmatch_features'} = [];
  }
  if(@args){
    push(@{$self->{'_pmatch_features'}}, @args);
  }
  return @{$self->{'_pmatch_features'}} if($self->{'_pmatch_features'});
}


sub output{
  my ($self, @output) = @_;

  if(!$self->{'_output'}){
    $self->{'_output'} = [];
  }
  if(@output){
    push(@{$self->{'_output'}}, @output);
  }
  
  return @{$self->{'_output'}};
}
###############
#other methods#
###############




sub uniquify{
  my ($self, @features) = @_;
  my $file = "/tmp/$$.pm.out";
  my @output;
 
 # print STDERR "UNIQUFY\n";
  open (OUT, ">".$file);
  foreach my $hit(@features) {
    print OUT $hit->chr_name  . ":" .$hit->start . "," .$hit->end   . ":" . $hit->protein_id . ":" .$hit->coverage . "\n";
  #  print STDERR $hit->chr_name  . ":" .$hit->start . "," .$hit->end   . ":" . $hit->protein_id . ":" .$hit->coverage . "\n";
  }
  #print STDERR "UNIQUFY\n";
  close(OUT);
  my $sorted_file = "/tmp/$$.pm.out.sorted";
  my $command = "sort -u ".$file." > ".$sorted_file;
  system($command);
  
  open(PMATCHES, "<$sorted_file") or die "couldn't open ".$sorted_file." $!";
  while(<PMATCHES>){
    #print;
    chomp;
    # eg chr22:10602496,10603128:Q9UGV6:99.4
    # cb25.fpc0002:136611,136883:CE00763:47.7
    if(!/(\S+):(\d+),(\d+):(\S+):(\S+)/){
      die "Cannot parse [$_]\nClean out protein and pmatch_feature tables before rerunning script!\n";
    }
    my $chr = $1;
    my $start = $2;
    my $end = $3;
    my $protein = $4;
    my $coverage = $5;
    #print STDERR "chr ".$chr." start ".$start." end ".$end." protein ".$protein." coverage ".$coverage."\n";
    if($start > $end){
      $start = $3;
      $end = $2;
    }

   
    my $cdna_id = undef;
    my $pmf = Bio::EnsEMBL::Pipeline::Tools::Pmatch::PmatchFeature->new(-protein_id  => $protein,
								 -start       => $start,
								 -end         => $end,
								 -chr_name    => $chr,
								 -cdna_id     => $cdna_id,
								 -coverage    => $coverage,
								 -analysis => $self->analysis,
								);
    push(@output, $pmf);
  }
  unlink $file;
  unlink $sorted_file;
  return @output;
}






1;
