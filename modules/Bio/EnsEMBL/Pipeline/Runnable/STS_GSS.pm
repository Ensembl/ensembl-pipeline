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
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);



sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_query'}     = undef;     # location of Bio::Seq object
    $self->{'_unmasked'} = undef;
    $self->{'_blast_program'}   = undef;     # location of Blast
    $self->{'_est_program'}   = undef;     # location of Blast
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
    $self->{'_merged_features'} =[];
    # Now parse the input options and store them in the object
    #print "@args\n";
    my( $query, $unmasked, $program, $database, $threshold, $threshold_type, $filter,$coverage,$prune,$options, $seqfetcher) = 
	    $self->_rearrange([qw(QUERY
				  UNMASKED
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
    if ($unmasked) {
      $self->unmasked($unmasked);
    } else {
      $self->throw("No unmasked query sequence input.");
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


sub unmasked {
    my ($self, $seq) = @_;
    if ($seq) {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
        $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
      }

      $self->{'_unmasked'} = $seq ;

     
    }
    return $self->{'_unmasked'};
}

sub program {
  my ($self, $location) = @_;
  
  if ($location) {
    $self->throw("executable not found at $location: $!\n")     unless (-e $location && -x $location);
    $self->{'_blast_program'} = $location ;
  }
  return $self->{'_blast_program'};
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
sub prune {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_prune} = $arg;
  }
  return $self->{_prune};
}

sub coverage {
  my($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_coverage} = $arg;
  }
  return $self->{_coverage};
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

sub add_merged_feature{
  my ($self, $feature) = @_;
  if(!$feature){
    $self->throw("no feature passed : $!\n");
  }else{
    push(@{$self->{'_merged_features'}}, $feature);
  }
}

sub each_merged_feature{

  my ($self) = @_;

  return@{$self->{'_merged_features'}};

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
  my @blast_output = $self->run_blasts();
  #print STDERR "ran analysis\n";
  #parse output and create features
  $self->parse_features(\@blast_output);
  my @features = $self->each_merged_feature;
  my @results;
  foreach my $feature(@features){
   # print "running on ".$feature->hseqname." with score ".$feature->score." and evalue ".$feature->p_value."\n";
    my @genes = $self->run_est2genome($feature);
    push(@results, @genes);
  }
  $self->deletefiles();
  $self->update_analysis(\@results);
  $self->output(\@results);
  
}



sub run_blasts {
   
  my ($self) = @_;

  my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(-query => $self->clone,
							   -program => $self->program,
							   -database => $self->database,
							   -threshold => $self->threshold,
							   -threshold_type => $self->threshold_type,
							   -filter => $self->filter,
							   -coverage => $self->coverage,
							   -prune => $self->prune,
							   -options => $self->options,
							   
							  );

  $blast->run;
  my @features = $blast->output;

  return @features;
							    
}




sub parse_features {
  my ($self, $features) = @_;

 
  
  my @features = @$features;
 
  
  #return @features;
  #print "there are ".scalar(@features)."\n";
  my %unique_hids;

  foreach my $feature(@features){
    my $hseqname = $feature->hseqname();
    $unique_hids{$hseqname} ||= [];
    push(@{$unique_hids{$hseqname}}, $feature);
  }
  
  my @hids = keys(%unique_hids);

  foreach my $hid(@hids){
    #  print "feature ".$hid." has ".scalar(@{$unique_hids{$hid}})." hits\n";
    my @feature_array = @{$unique_hids{$hid}};
    my $hid_seq = $self->seqfetcher->get_Seq_by_acc($hid);
    my $hid_len = $hid_seq->length;
    foreach my $feature(@feature_array){
      #print "before seqname = ".$hid." length = ".$hid_len." start ".$feature->start." end ".$feature->end." strand ".$feature->strand." hstart ".$feature->hstart." hend ".$feature->hend."\n";
      my $genomic_start = $feature->start;
      my $genomic_end = $feature->end;
      my $hstart = $feature->hstart;
      my $hend = $feature->hend;
     if($feature->strand == -1){
	$feature->start($genomic_start-($hid_len+200));
	if($feature->start <= 0){
	  $feature->start(1);
	}
	$feature->end($genomic_end+($hstart+200));
	if($feature->end >= $self->clone->length){
	  $feature->end($self->clone->length);
	}
	$feature->hstart(1);
	$feature->hend($hid_len);
      }elsif($feature->strand == 1){
	$feature->start($genomic_start-($hstart+200));
	if($feature->start <= 0){
	  $feature->start(1);
	}
	$feature->end($genomic_end+($hid_len+200));
	if($feature->end >= $self->clone->length){
	  $feature->end($self->clone->length);
	}
	$feature->hstart(1);
	$feature->hend($hid_len);
      }else{
	$self->throw($feature->hseqname." has got ".$feature->strand." as a stand : $!\n");
      }
  
      #print "after seqname = ".$hid."start ".$feature->start." end ".$feature->end." hstart ".$feature->hstart." hend ".$feature->hend." strand ".$feature->strand."\n";
    }
   ####needs to merge overlapping sequences#### 
    my @sorted_forward_features;
    my @sorted_reverse_features;
    foreach my $feature(@feature_array){
      if($feature->strand == -1){
	push(@sorted_reverse_features, $feature)
      }elsif($feature->strand == 1){
	push(@sorted_forward_features, $feature)
      }else{
	$self->throw($feature->hid." hasn't got a stand : $!\n");
      }
    }
    @sorted_forward_features = sort{$a->start <=> $b->start || $a->end <=> $b->end} @sorted_forward_features;
    @sorted_reverse_features = sort{$a->start <=> $b->start || $a->end <=> $b->end} @sorted_reverse_features;
    $self->fuse_overlapping_features(\@sorted_forward_features);
    $self->fuse_overlapping_features(\@sorted_reverse_features);
  }

  my @merged = $self->each_merged_feature;

  #print "there are ".scalar(@merged)." features\n";
}

sub fuse_overlapping_features {
    
  my ($self, $feature_ref) = @_;
  
  my @features = @$feature_ref;
  
    for (my $i = 1; $i < @features;) {
        my $f_a = $features[$i - 1];
        my $f_b = $features[$i];
        
        # If coordinates overlap or abut
        if ($f_a->end + 1 >= $f_b->start) {
            # Transfer the end of b to a if b ends beyond a
            if ($f_a->end < $f_b->end) {
                $f_a->end($f_b->end);
            }
            
            # Remove b from the list
            splice(@features, $i, 1);
            
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
  
  foreach my $feature (@features){

    $self->add_merged_feature($feature);

  }
  
  


}

sub run_est2genome{

  my ($self, $feature) = @_;
  
  my $est = $self->seqfetcher->get_Seq_by_acc($feature->hseqname);
  my $seq = $self->unmasked->subseq($feature->start, $feature->end);
  my $start = $feature->start;
  my $strand = $feature->strand;
  my $end = $feature->end;
  my $query_seq = Bio::Seq->new (-seq => $seq,
				 -id => $self->clone->id,
				 -accession => $self->clone->id);
  print "running est2genome on ".$feature->hseqname." strand ".$feature->strand."\n";
  my $est2genome = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new(-genomic => $query_seq,
								     -est=> $est,
								    );
  $est2genome->run;
  my @features = $est2genome->output;
  foreach my $f(@features){
    $f->strand($strand);
    print "ran est2genome on ".$feature->hseqname." strand ".$feature->strand."\n";
    
      $f->start($f->start+$start-1);
      $f->end($f->end+$end-1);
    
      
  }
  print "there are ".scalar(@features)." outputted\n";
								     
  return @features;
}

sub update_analysis{
  my($self, $feature) = @_;

  my $analysis =  Bio::EnsEMBL::Analysis->new (-db        => 'dbGSS',
					       -db_versions => 1,
					       -program        => 'wublastn',
					       -program_version   => 1,
					       -module         => 'STS_GSS',
					       -module_version => 1,
					       -gff_source     => 'STS_GSS',
					       -gff_feature    => 'similarity', 
					       -logic_name     => 'Full_dbGSS',
					       -parameters     => '-p1 B=100000 V=500 E=0.001 Z=500000000',
				      );


  my @features = @$feature;
		 
  foreach my $f(@features){
    $f->source_tag("STS_GSS");
    $f->feature1->analysis($analysis);
    $f->feature1->source_tag("STS_GSS");
    $f->feature2->analysis($analysis);
    $f->feature2->source_tag("STS_GSS");
    foreach my $e($f->feature1->sub_SeqFeature){
      $e->feature1->analysis($analysis);
      $e->feature1->source_tag("STS_GSS");
      $e->feature2->analysis($analysis);
      $e->feature2->source_tag("STS_GSS");
      foreach my $s($e->feature1->sub_SeqFeature){
	$s->feature1->analysis($analysis);
	$s->feature1->source_tag("STS_GSS");
	$s->feature2->analysis($analysis);
	$s->feature2->source_tag("STS_GSS");
      }
    }
  }
    
}

sub output{

  my ($self, $output) = @_;
  if($output){
    my @outputs = @$output;
    push(@{$self->{'_fplist'}}, @outputs); 
  }

  return @{$self->{'_fplist'}};
}

1;
