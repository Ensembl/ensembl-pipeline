#module for running pmatch analysis

package Bio::EnsEMBL::Pipeline::RunnableDB::Pmatch;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General;
use Bio::EnsEMBL::Pipeline::Runnable::Pmatch;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::PmatchFeature;
use Bio::SeqIO;   

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($prot_file, $max_intron_length ) =  $self->_rearrange([qw(
											     PROTEIN_FILE
											     MAX_INTRON_LENGTH
											    )], @args);

  $prot_file = $GB_PFASTA unless($prot_file);
  $self->throw("Can't run pmatch without a protein file : $!") unless($prot_file);
  $self->protein_file($prot_file);
  $max_intron_length = $GB_PMATCH_MAX_INTRON unless($max_intron_length);
  $self->throw("Can't run pmatch without a directory of max intron length : $!") unless($max_intron_length);
  $self->max_intron($max_intron_length);

  return $self;
}




sub fetch_input{
  my ($self) = @_;
  #print STDERR "Fetching input\n";
  $self->make_protlist;
  #print STDERR "Trying to parse ".$self->input_id." with ".$GB_INPUTID_REGEX."\n";
  $self->fetch_sequence($GB_PMATCH_MASKING);
  #print STDERR "Have sequence ".$self->query->name."\n";
  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Pmatch->new(
							       '-query' => $self->query,
							       '-protein_file' => $self->protein_file,
							       '-options' => $self->analysis->parameters,
							       '-max_intron_size' => $self->max_intron,
							       '-protein_lengths' => $self->prot_lengths,
							      );
  $self->runnable($runnable);
  #print STDERR "have created runnable\n";
}
 


sub run{
  my ($self) = @_;
  $self->runnable->run;
  my @hits = $self->runnable->output;
  #prin STDERR "Have ".@hits." output from pmatch\n";
  my @output;
  foreach my $hit(@hits){
    #print STDERR "have ".$hit."\n";
    my ($name, $start, $end) = $self->convert_coords_to_genomic($hit->query,
                                                                $hit->qstart,
                                                                $hit->qend);
    $hit->query($name);
    $hit->qstart($start);
    $hit->qend($end);
    push(@output, $hit);
  }
  #print STDERR "Have ".@output." output from conversion\n";
  $self->convert_output(\@output);
}



sub write_output{
  my ($self) = @_;
  my $pmfa = $self->db->get_PmatchFeatureAdaptor();
  #print STDERR "Writing to ".$self->pipeline_db->dbname." on ".$self->pipeline_db->host."\n";
  $pmfa->write_PmatchFeatures($self->output);
}



##################
#accessor methods#
##################

sub prot_lengths{
  my ($self, $arg) = @_;

  if(!$self->{'_prot_len'}){
    $self->{'_prot_len'} = {};
  }
  if($arg){
    if(ref($arg) ne 'HASH'){
      $self->throw("the protein lengths must be in a hash reference\n");
    }
    $self->{'_prot_len'} = $arg;
  }
  #print "have ".$self->{'_prot_len'}." for prot lengths\n";
  return $self->{'_prot_len'};
}

sub protein_file{
  my ($self, $file) =@_;

  if(!$self->{'_prot_file'}){
    $self->{'_prot_file'} = undef;
  }
  if($file){
    $self->{'_prot_file'} = $file;
  }

  return $self->{'_prot_file'};
}


sub max_intron{
  my ($self, $file) =@_;

  if(!$self->{'_max_intron'}){
    $self->{'_max_intron'} = undef;
  }
  if($file){
    $self->{'_max_intron'} = $file;
  }

  return $self->{'_max_intron'};
}


sub runnable{
  my ($self, $runnable) = @_;

  if(!$self->{'_runnables'}){
    $self->{'_runnables'} = undef;
  }
  if($runnable){
    $self->{'_runnables'} = $runnable;
  }

  return $self->{'_runnables'};
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




sub make_protlist{
  my ($self) = @_;
  ##print STDERR "making prot_length list from $protfile\n";
  my %plengths;
#  open (INFILE, "<".$self->protein_file) or die "Can't open ".$self->protein_file."\n";
#  $/ = '>'; # looks for fields separated by > instead of "\n"
# SEQFETCH:
#  while(<INFILE>) {
    
#    my @lines = split /\n/;
#    next SEQFETCH unless scalar(@lines) > 1;
#    my $description = shift (@lines);
#    $description =~ /(\S+)/; # find the first block of non-whitespace

#    next SEQFETCH unless defined $description and $description ne '';
#    my $bs = new Bio::Seq;
#    $bs->display_id($1);
    
#    my $seq;
#    foreach my $line(@lines) {
#      chomp $line;
#      $seq .= $line;
#    }
    
#    $bs->seq($seq);
#    #print STDERR $bs->display_id ." = ".$bs->length."\n";
#    $plengths{$bs->display_id} = $bs->length;
#  }

#  # note double quotes ...    
#  $/ = "\n";
#  close INFILE;

  my $seqio = new Bio::SeqIO(-format => 'Fasta', 
			     -file => $self->protein_file);
  while(my $seq = $seqio->next_seq){
    #print STDERR "have id ".$seq->id." ".$seq->length."\n";
    $plengths{$seq->id} = $seq->length;
  }
  my $ref = \%plengths;
  #print STDERR "have ".$ref." hash \n";
  $self->prot_lengths($ref);
}


sub convert_coords_to_genomic{
  my ($self, $name, $start, $end) = @_;
  #print STDERR "Have ".$name." ".$start." ".$end."\n";
  my $chr_name = $self->query->seq_region_name;
  my $chr_start = $self->query->start;
  my $chr_end = $self->query->end; 
  if($start - 1 == 0){
    return ($chr_name, $start, $end);
  }
  my $genomic_start = $start+$chr_start-1;
  my $genomic_end = $end+$chr_start-1;
  # print STDERR "Have ".$chr_name." ".$genomic_start." ".$genomic_end."\n";
  return($chr_name, $genomic_start, $genomic_end);
}


sub convert_output{
  my ($self, $output) = @_;

  my @results = @$output;
  
  
  my @out = $self->uniquify(@results);
  foreach my $f(@out){
     $f->chr_name($self->query->seq_region_name);
  }
  $self->output(@out);
}


sub uniquify{
  my ($self, @features) = @_;
  my $file = "/tmp/$$.pm.out";
  my @output;
 
 # print STDERR "UNIQUFY\n";
  open (OUT, ">".$file);
  foreach my $hit(@features) {
    #print OUT $hit->query  . ":" .$hit->qstart . "," .$hit->qend   . ":" . $hit->target . ":" .$hit->coverage . "\n";
  #  print STDERR $hit->query  . ":" .$hit->qstart . "," .$hit->qend   . ":" . $hit->target . ":" .$hit->coverage . "\n";
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
