package Bio::EnsEMBL::Pipeline::Runnable::STS_GSS;


use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use FileHandle;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use Bio::Tools::BPlite;
use Bio::Tools::BPlite::HSP;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);



sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_query'}     = undef;     # location of Bio::Seq object
    $self->{'_program'}   = undef;     # location of Blast
    $self->{'_database'}  = undef;     # name of database
    $self->{'_threshold'} = undef;     # Threshold for hit filterting
    $self->{'_options'}   = undef;     # arguments for blast
    $self->{'_filter'}    = 1;         # Do we filter features?
    $self->{'_fplist'}    = [];        # an array of feature pairs (the output)

    $self->{'_workdir'}   = undef;     # location of temp directory
    $self->{'_filename'}  = undef;     # file to store Bio::Seq object
    $self->{'_results'}   = undef;     # file to store results of analysis
    $self->{'_seqfetcher'} = undef;
    $self->{'_prune'}     = 1;         # 
    $self->{'_coverage'}  = 10;
    $self->{'_merged_hsps'} =[];
    # Now parse the input options and store them in the object
    #print "@args\n";
    my( $query, $program, $database, $threshold, $threshold_type, $filter,$coverage,$prune,$options, $seqfetcher) = 
	    $self->_rearrange([qw(QUERY 
				  PROGRAM 
				  DATABASE 
				  THRESHOLD
				  THRESHOLD_TYPE
				  FILTER 
				  COVERAGE
				  PRUNE
				  OPTIONS
				  SEQFETCHER)], 
			      @args);

    
    if ($query) {
      $self->clone($query);
    } else {
      $self->throw("No query sequence input.");
    }
    #print STDERR "find executable output ".$self->find_executable($program)."\n";
    $self->program($self->find_executable($program));
    if ($database) {
      $self->database($database);
    } else {
      $self->throw("No database input");
    }
    
    if ($options) {
#this option varies the number of HSP displayed proportionally to the query contig length
          if ($::pipeConf{'B_factor'}){
                my $b_factor = $::pipeConf{'B_factor'};
                my $b_value = int ($query->length / 1000 * $b_factor); 
                if ($::pipeConf{'blast'} eq 'ncbi'){
                        $options .= " -b $b_value" unless ($b_value < 250);
                }
                else {
                        $options .= " B=$b_value" unless ($b_value < 250); 
                }
          } 
      $self->options($options);
    } else {
      $self->options(' -p1 ');  
    }
    #print "options = ".$self->options."\n";
    if (defined($threshold)) {
      $self->threshold($threshold);
    }

    if (defined($threshold_type)) {
      $self->threshold_type($threshold_type);
    }

    if (defined($filter)) {
        $self->filter($filter);
    }
    
    if (defined($prune)) {
      $self->prune($prune);
    }
    if (defined($coverage)) {
      $self->coverage($coverage);
    }

    $self->throw("No seqfetcher provided")           
    unless defined($seqfetcher);

    $self->seqfetcher($seqfetcher) if defined($seqfetcher);
    return $self; # success - we hope!
}

#################
#get/set methods#
#################

sub clone {
    my ($self, $seq) = @_;
    if ($seq) {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
        $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
      }

      $self->{'_query'} = $seq ;

      $self->filename($self->clone->id.".$$.seq");
      $self->results($self->filename.".blast.out");
      
    }
    return $self->{'_query'};
}


sub program {
  my ($self, $location) = @_;
  
  if ($location) {
    $self->throw("executable not found at $location: $!\n")     unless (-e $location && -x $location);
    $self->{'_program'} = $location ;
  }
  return $self->{'_program'};
}


sub database {
    my ($self, $db) = @_;

    if (defined($db)) {
      $self->{'_database'} = $db ;
    }
    return $self->{'_database'};
  }



sub options {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_options'} = $args ;
  }
  return $self->{'_options'};
}

sub filter {
    my ($self,$args) = @_;

    if (defined($args)) {
        if ($args != 0 && $args != 1) {
            $self->throw("Filter option must be 0 or 1");
        }
        $self->{'_filter'} = $args;
    }
    return $self->{'_filter'};
}

sub get_threshold_types {
  my ($self) = @_;

  return ("PID","PVALUE");
}

sub threshold_type {
  my ($self,$type) = @_;

  my @types = $self->get_threshold_types;
  
  if (defined($type)) {
    my $found = 0;
    foreach my $allowed_type ($self->get_threshold_types) {
      if ($type eq $allowed_type) {
        $found = 1;
      }
    }
    if ($found == 0) {

      $self->throw("Type [$type] is not an allowed type.  Allowed types are [@types]");
    } else {
      $self->{_threshold_type} = $type;
    }
  }
  return $self->{_threshold_type} || $types[0];
}

sub get_pars {
  my ($self) = @_;

  if (!defined($self->{_hits})) {
     $self->{_hits} = [];
  }
	
  return @{$self->{_hits}};

}

sub seqfetcher {
  my( $self, $value ) = @_;    
  if ($value) {
    #need to check if passed sequence is Bio::DB::RandomAccessI object
#    $value->isa("Bio::DB::RandomAccessI") || $self->throw("Input isn't a Bio::DB::RandomAccessI");
    $self->{'_seqfetcher'} = $value;
  }
  return $self->{'_seqfetcher'};
}

sub add_merged_hsp{
  my ($self, $hsp) = @_;
  if(!$hsp){
    $self->throw("no hsp passed : $!\n");
  }else{
    push(@{$self->{'_merged_hsps'}}, $hsp);
  }
}

sub each_merged_hsp{

  my ($self) = @_;

  return@{$self->{'_merged_hsps'}};

}

###########
# Analysis methods
##########



sub run {
  my ($self, $dir) = @_;

  my $seq = $self->clone || $self->throw("Query seq required for Blast\n");

  $self->workdir('/tmp') unless ($self->workdir($dir));
  $self->checkdir();
  #print STDERR "checked directories\n";
  #write sequence to file
  $self->writefile(); 
  $self->run_analysis();
  #print STDERR "ran analysis\n";
  #parse output and create features
  $self->parse_results;
  $self->deletefiles();

  #return @hsps;
}



sub run_analysis {
   
  my ($self) = @_;
  
  # This routine expands the database name into $db-1 etc for
  # split databases

  my @databases = $self->fetch_databases;
  
  $self->database(@databases);

  foreach my $database (@databases) {
    my $db = $database;
    print "database ".$database."\n";
    $db =~ s/.*\///;
    #allow system call to adapt to using ncbi blastall. defaults to WU blast.	
    my $command = $self->program ;
    $command .= ($::pipeConf{'blast'} eq 'ncbi') ? ' -d '.$database : ' '.$database;
    $command .= ($::pipeConf{'blast'} eq 'ncbi') ? ' -i ' .$self->filename :  ' '.$self->filename;
    $command .= ' '.$self->options. ' > '.$self->results . ".$db";
    print $command."\n";
    $self->throw("Failed during blast run $!\n") unless (system ($command) == 0) ;
  }
}



sub fetch_databases {
  my ($self) = @_;
  
  my @databases;
    
  my $fulldbname;

  # If we have passed a full path name don't append the $BLASTDB
  # environment variable.
  
  if ($self->database =~ /\//) {
    $fulldbname = $self->database;
  } else {
    $fulldbname = $ENV{BLASTDB} . "/" .$self->database;
  }
  
  # If the expanded database name exists put this in
  # the database array.
  #
  # If it doesn't exist then see if $database-1,$database-2 exist
  # and put them in the database array
  
  if (-e $fulldbname) {
    push(@databases,$self->database);
  } else {
    my $count = 1;
    
    while (-e $fulldbname . "-$count") {
      push(@databases,$fulldbname . "-$count");
      $count++;
    }
  }
  
  if (scalar(@databases) == 0) {
    $self->throw("No databases exist for " . $self->database);
  }

  return @databases;

}

sub get_parsers {
  my ($self)  = @_;

  my @parsers;

  foreach my $db ($self->database) {
    $db =~ s/.*\///;

    my $fh = new FileHandle;
    $fh->open("<" . $self->results . ".$db");
    
    my $parser = new Bio::Tools::BPlite ('-fh' => $fh);
    
    push(@parsers,$parser);
  } 

  return @parsers;
}


sub parse_results {
  my ($self) = @_;

  my %ids;

  my @parsers;

  
  @parsers = $self->get_parsers;
  
  my @hsps;
  foreach my $parser (@parsers) {
    while (my $sbjct = $parser->nextSbjct) {
      while (my $hsp = $sbjct->nextHSP) {
	push(@hsps,$hsp);
      }
    }
  }
  
  #return @hsps;
  print "there are ".scalar(@hsps)."\n";
  my %unique_hids;

  foreach my $hsp(@hsps){
    my @hseqname = split /\|/, $hsp->hseqname();
    $unique_hids{$hseqname[1]} ||= [];
    push(@{$unique_hids{$hseqname[1]}}, $hsp);
  }
  
  my @hids = keys(%unique_hids);

  foreach my $hid(@hids){
    #  print "hsp ".$hid." has ".scalar(@{$unique_hids{$hid}})." hits\n";
    my @hsp_array = @{$unique_hids{$hid}};
    my $hid_seq = $self->seqfetcher->get_Seq_by_acc($hid);
    my $hid_len = $hid_seq->length;
    foreach my $hsp(@hsp_array){
#      #print "before seqname = ".$hid." length = ".$hid_len." \nstart ".$hsp->start." end ".$hsp->end." strand ".$hsp->strand." hstart ".$hsp->hstart." hend ".$hsp->hend."\n";
      my $genomic_start = $hsp->start;
      my $genomic_end = $hsp->end;
      my $hstart = $hsp->hstart;
      my $hend = $hsp->hend;
     if($hsp->strand == -1){
	$hsp->start($genomic_start-($hid_len+100));
	if($hsp->start <= 0){
	  $hsp->start(1);
	}
	$hsp->end($genomic_end+($hstart+100));
	if($hsp->end >= $self->clone->length){
	  $hsp->end($self->clone->length);
	}
	$hsp->hstart(1);
	$hsp->hend($hid_len);
      }elsif($hsp->strand == 1){
	$hsp->start($genomic_start-($hstart+100));
	if($hsp->start <= 0){
	  $hsp->start(1);
	}
	$hsp->end($genomic_end+($hid_len+100));
	if($hsp->end >= $self->clone->length){
	  $hsp->end($self->clone->length);
	}
	$hsp->hstart(1);
	$hsp->hend($hid_len);
      }else{
	$self->throw($hsp->hseqname." has got ".$hsp->strand." as a stand : $!\n");
      }
  
      #print "after seqname = ".$hid."\nstart ".$hsp->start." end ".$hsp->end." hstart ".$hsp->hstart." hend ".$hsp->hend."\n";
    }
   ####needs to merge overlapping sequences#### 
    my @sorted_forward_hsps;
    my @sorted_reverse_hsps;
    foreach my $hsp(@hsp_array){
      if($hsp->strand == -1){
	push(@sorted_reverse_hsps, $hsp)
      }elsif($hsp->strand == 1){
	push(@sorted_forward_hsps, $hsp)
      }else{
	$self->throw($hsp->hid." hasn't got a stand : $!\n");
      }
    }
    @sorted_forward_hsps = sort{$a->start <=> $b->start || $a->end <=> $b->end} @sorted_forward_hsps;
    @sorted_reverse_hsps = sort{$a->start <=> $b->start || $a->end <=> $b->end} @sorted_reverse_hsps;
    $self->fuse_overlapping_hsps(\@sorted_forward_hsps);
    $self->fuse_overlapping_hsps(\@sorted_reverse_hsps);
  }

  my @merged = $self->each_merged_hsp;

  print "there are ".scalar(@merged)." hsps\n";
}

sub fuse_overlapping_hsps {
    
  my ($self, $hsp_ref) = @_;
  
  my @hsps = @$hsp_ref;
  
    for (my $i = 1; $i < @hsps;) {
        my $f_a = $hsps[$i - 1];
        my $f_b = $hsps[$i];
        
        # If coordinates overlap or abut
        if ($f_a->end + 1 >= $f_b->start) {
            # Transfer the end of b to a if b ends beyond a
            if ($f_a->end < $f_b->end) {
                $f_a->end($f_b->end);
            }
            
            # Remove b from the list
            splice(@hsps, $i, 1);
            
            # Critical error check!
            if ($f_a->start > $f_a->end) {
                die "Feature start (", $f_a->start,
                    ") greater than end (", $f_a->end, ")";
            }
        } else {
            # Haven't spliced, so move pointer
            $i++;
        }
    }
  
  foreach my $hsp (@hsps){

    $self->add_merged_hsp($hsp);

  }
  
  


}



1;
