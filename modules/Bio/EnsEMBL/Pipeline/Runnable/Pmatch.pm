# this runnable will pmatch proteins against genomic both the commandline construction and the parsing assumes this is what
# is being run for more information on how we run the code see running_pmatchs.txt in ensembl-docs

package Bio::EnsEMBL::Pipeline::Runnable::Pmatch;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::Pmatch::First_PMF;
use Bio::EnsEMBL::Root;
use File::Copy;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  #print STDERR "PMATCH args : @args\n";
  my ($sequence, $proteins, $program, $options, $max_intron_size, $protein_lengths) = $self->_rearrange([qw( 
												     QUERY
												     PROTEIN_FILE
												     PROGRAM
												     OPTIONS
												     MAX_INTRON_SIZE
												     PROTEIN_LENGTHS
											 )], @args);

  #print "sequence ".$sequence."\n";
  $self->throw("Can't run with out a query sequence") unless($sequence);
  $self->query($sequence);
  $self->throw("Can't run with out a peptide file ") unless($proteins);
  $self->protein_file($proteins);
  unless($program){
    $program = 'pmatch';
  }
  $self->program($self->find_executable($program));
  $self->options($options);
  $self->max_intron_size($max_intron_size);
  if((!$protein_lengths) ||  ref($protein_lengths) ne "HASH"){
    $self->throw("Can't run without a hash of protein lengths\n");
  }
  $self->prot_lengths($protein_lengths);
  return $self;
}

##########
#analysis#
##########

sub run{
   my ($self) = @_;
   #print STDERR "Running pmatch on ".$self->sequence_file."\n";
   $self->workdir('/tmp') unless $self->workdir();
   $self->checkdir();
   $self->writefile;
   $self->run_analysis;
   $self->parse_results;
   $self->deletefiles;

   1;
}



sub run_analysis{
  my ($self) = @_;

  
  my $command = $self->program." -D ".$self->protein_file." ".$self->filename." > ".$self->results;
  print STDERR $command."\n";
  $self->throw("Error running pmatch on " . $self->filename) if system($command);
}

sub parse_results{
  my ($self) = @_;
  
  $self->sort_results;
  
}




sub sort_results{
  my ($self) = @_;
  my $sorted_file = $self->query->name.".$$.sorted_results";
  print STDERR "Sorted file ".$sorted_file."\n";
  my $command = "sort -k6,6 -k3,3n ".$self->results ." > ".$sorted_file;
  $self->throw("Error [sorting results] running pmatch on " . $self->filename) if system($command);
  rename($sorted_file, $self->results);
  my $current_pmf;
  my $prot_id;
  open (PM, "<".$self->results) or die "couldn't open ".$self->results."\n";
  
  #print STDERR "Prot lengths ".$self->prot_lengths."\n";
  my $line_count = 0;
  my $pmf_count = 0;
 PMATCH:  
  while(<PM>){
    $line_count++;
    #print STDERR;
    my @cols = split;
    # dump out line to file just in case this turns into a problem
    
    if(!$prot_id || $cols[5] ne $prot_id){

      $self->add_pm_filter($current_pmf);
      $pmf_count++;
      # start a new PMF
      $current_pmf = new Bio::EnsEMBL::Pipeline::Tools::Pmatch::First_PMF(
								   -plengths => $self->prot_lengths,
								   -maxintronlen => $self->max_intron_size,
				  );
      $prot_id = $cols[5];
      #print STDERR "working with ".$prot_id."\n";
      $current_pmf->make_coord_pair($_);
    }else{
      # add this hit into current PMF
      $current_pmf->make_coord_pair($_);
    }
  }
  #print STDERR "Have ".$line_count." lines from results\n";
  $pmf_count++;
  #print STDERR "Have ".$pmf_count." First_PMF's\n";
  $self->add_pm_filter($current_pmf);
  # make sure we at least try to proces the last one!
  
  my @filters = @{$self->get_pm_filters};
  #print STDERR "Have ".@filters." filters\n";
  foreach my $f(@filters){
    my @hits = $f->merge_hits;
    #print STDERR "adding ".$hits[0]." to output\n";
    $self->add_merged_hits(@hits);
  }
}







###########
#accessors#
###########

sub query {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_query'} = $seq ;
        $self->filename($seq->id.".$$.seq");
        $self->results($self->filename.".pmatch");
	$self->file($self->results);
    }
    return $self->{'_query'};
}



sub add_pm_filter{
  my ($self, $arg) = @_;
  
  if(!$self->{'_pm_filters'}){
    $self->{'_pm_filters'} = [];
  }

  if($arg){
    push(@{$self->{'_pm_filters'}}, $arg);
  }
}

sub get_pm_filters{
  my ($self) = @_;

  return $self->{'_pm_filters'};
}

sub add_merged_hits{
  my ($self, @arg) = @_;
  
  if(!$self->{'_output'}){
    $self->{'_output'} = [];
  }
  if(@arg){
    if(!$arg[0]->isa("Bio::EnsEMBL::Pipeline::Tools::Pmatch::MergedHit")){
      $self->throw("the first element of arg isn't a MergedHit but a ".$arg[0]." this can't be used for output\n");
    }
    push(@{$self->{'_output'}}, @arg);
  }
}

sub output{
  my ($self) = @_;

  my @output = @{$self->{'_output'}} if($self->{'_output'});
  return @output;
}

sub sequence_file{
  my ($self, $arg) = @_;

  if(!$self->{'_sequence'}){
    $self->{'_sequence'} = undef;
  }
  if($arg){
    if(! -e $arg){
      $self->throw($arg." should exist : $!");
    }else{
      $self->{'_sequence'} = $arg;
    }
  }

  return $self->{'_sequence'};
}

sub protein_file{
  my ($self, $arg) = @_;

  if(!$self->{'_protein'}){
    $self->{'_protein'} = undef;
  }
  if($arg){
    if(! -e $arg){
      $self->throw($arg." should exist : $!");
    }else{
      $self->{'_protein'} = $arg;
    }
  }

  return $self->{'_protein'};
}

sub program{
  my($self, $arg) = @_;

  if(!$self->{'_program'}){
    $self->{'_program'} = undef;
  }

  if($arg){
    $self->throw("Pmatch not found at $arg: $!\n") unless (-e $arg);
    $self->{'_program'} = $arg;
  }

  return $self->{'_program'};
}

sub options{
  my($self, $arg) = @_;

  if(!$self->{'_options'}){
    $self->{'_options'} = undef;
  }

  if($arg){
    $self->{'_options'} = $arg;
  }

  return $self->{'_options'};
}

sub max_intron_size{
  my($self, $arg) = @_;

  if(!$self->{'_intron'}){
    $self->{'_intron'} = undef;
  }

  if($arg){
    $self->{'_intron'} = $arg;
  }

  return $self->{'_intron'};
}


sub filename{
  my($self, $arg) = @_;

  if(!$self->{'_filename'}){
    $self->{'_filename'} = undef;
  }

  if($arg){
    $self->{'_filename'} = $arg;
  }

  return $self->{'_filename'};
}


sub prot_lengths{
  my ($self, $arg) = @_;

  if(!$self->{'_prot_len'}){
    $self->{'_prot_len'} = undef;
  }
  if($arg){
    if(ref($arg) ne 'HASH'){
      $self->throw("the protein lengths must be in a hash reference\n");
    }
    $self->{'_prot_len'} = $arg;
  }

  return $self->{'_prot_len'};
}



1;
