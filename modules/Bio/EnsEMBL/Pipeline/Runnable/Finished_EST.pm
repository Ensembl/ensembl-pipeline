
### Bio::EnsEMBL::Pipeline::Runnable::Finished_EST

package Bio::EnsEMBL::Pipeline::Runnable::Finished_EST;

use strict;
use Data::Dumper;

use vars qw(@ISA);
use Bio::EnsEMBL::Pipeline::RunnableI;

use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
use Bio::EnsEMBL::Pipeline::MiniSeq;
use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::DB::RandomAccessI;

use Bio::PrimarySeqI;
use Bio::SeqIO;


@ISA = ('Bio::EnsEMBL::Pipeline::RunnableI');

=head2 new

  Arg [1]   : hash of parameters 
  Function  : to create a new STS_GSS object
  Returntype: an STS_GSS object
  Exceptions: It will throw is it recieves no unmasked sequence, it will throw if it isn''t running a blast an no blast features are passed to it. If it is running a blast it will throw if it receives no blast db location or masked sequence. It will also throw if it gets no seqfetcher 
  Caller    : 
  Example   : 

=cut

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
    $self->{'_no_blast'} = 0;
    $self->{'_features'} = [];
    $self->{'_workdir'}   = undef;     # location of temp directory
    $self->{'_filename'}  = undef;     # file to store Bio::Seq object
    $self->{'_results'}   = undef;     # file to store results of analysis
    $self->{'_seqfetcher'} = undef;
    $self->{'_prune'}     = 1;         # 
    $self->{'_coverage'}  = 10;
    $self->{'_merged_features'} =[];
    $self->{'_percent_id'} = undef;
    $self->{'_padding'} = undef;
    $self->{'_percent_filter'} = undef;
    $self->{'_tandem'} = undef;
    # Now parse the input options and store them in the object
    #print "@args\n";
    my( $query, $unmasked, $program, $database, $threshold, $threshold_type, $filter,$coverage,$prune,$options, $seqfetcher, $features, $no_blast, $percent_id, $padding, $tandem) = 
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
				  SEQFETCHER
				  FEATURES
				  NO_BLAST
				  PERCENT_ID
				  PADDING
				  PERCENT_FILTER
				  TANDEM_CHECK)], 
			      @args);
    
    if($no_blast == 1){
      $self->{'_no_blast'} = $no_blast;
    } 
    
    if ($unmasked) {
      $self->unmasked($unmasked);
    } else {
      $self->throw("No unmasked query sequence input.");
    }
   
    if($program){
      $self->program($self->find_executable($program));
    }
    
      if (defined($features)) {
	if (ref($features) eq "ARRAY") {
	  push(@{$self->{'_features'}},@$features);
	} else {
	  $self->throw("[$features] is not an array ref.");
	}
      }else{
	$self->throw("should pass in features haven't \n");
      }
    
    
    if ($options) {
#this option varies the number of HSP displayed proportionally to the query contig length
      $self->options($options);
    } else {
      $self->options(' -p1 ');  
    }
    if ($percent_id) {

      $self->percent_id($percent_id);
    } else {
      $self->percent_id(95);  
    }
     if ($percent_id) {

      $self->padding($padding);
    } else {
      $self->padding(200);  
    }
    #print "options = ".$self->options."\n";
    if (defined($threshold)) {
      $self->threshold($threshold);
    } else {
          $self->threshold(95);  
    }

    if (defined($threshold_type)) {
      $self->threshold_type($threshold_type);
    }else {
      $self->threshold_type('PID');
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
    if (defined($tandem)) {
      $self->tandem_check($tandem);
    }
    $self->throw("No seqfetcher provided")           
    unless defined($seqfetcher);

    $self->seqfetcher($seqfetcher) if defined($seqfetcher);
    return $self; # success - we hope!
}

#################
#get/set methods#
#################


=head2 accessor methods

  Arg [1]   : variable to be set 
  Function  : if passed set the varible to the argument passed the return the argument
  Returntype: the varible to be set
  Exceptions: some throw exceptions if not passed the correct varible
  Caller    : $self
  Example   : $clone = $self->clone;

=cut


sub clone {
    my ($self, $seq) = @_;
    if ($seq) {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
        $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
      }

      $self->{'_query'} = $seq ;

      
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
      $self->filename($self->unmasked->id.".$$.seq");
      $self->results($self->filename.".blast.out");
     
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

sub no_blast {
    my ($self) = @_;

    return $self->{'_no_blast'};
}

sub options {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_options'} = $args ;
  }
  return $self->{'_options'};
}

sub percent_filter {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_percent_filter'} = $args ;
  }
  return $self->{'_percent_filter'};
}

sub percent_id {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_percent_id'} = $args ;
  }
  return $self->{'_percent_id'};
}

sub padding {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_padding'} = $args ;
  }
  return $self->{'_padding'};
}

sub prune {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_prune} = $arg;
  }
  return $self->{_prune};
}

sub tandem_check {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_tandem} = $arg;
  }
  return $self->{_tandem};
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
   $value->isa("Bio::DB::RandomAccessI") || $self->throw("Input isn't a Bio::DB::RandomAccessI");
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


sub features{
  my ($self) = @_;
  return @{$self->{'_features'}};
}

sub run {
    my( $self ) = @_;
    
    my $seq = $self->clone;
    my @raw_hits = $self->features;
    my %blast_ests;

    foreach my $raw_hit ( @raw_hits ) {
        
        my $seqname = $raw_hit->hseqname;
        
        push(@{$blast_ests{$seqname}}, $raw_hit);
        
 
    }
    
    while (my ($est, $features) = each %blast_ests) {
      
      # Is there a match below the P value threshold?
      my $good_match = 0;
      my $est_strand = 0;
      foreach my $res (@$features) {
        if ($res->p_value <= $self->threshold) {
            $est_strand = $res->strand;
            $good_match = 1;
            last;
        }
      }
      
      # Skip EST matches that don't have a match below the threshold
      next unless $good_match;
        
      #my $miniseq = $self->make_miniseq(@$features);
#       
#      print @$features[0]->hseqname,"\nMiniseq_sequence \n",$miniseq->get_cDNA_sequence->seq,"\n";
#      
#      print "EST Sequence\n",$self->seqfetcher->get_Seq_by_acc(@$features[0]->hseqname)->seq,"\n\n";
      
       #make MiniEst2Genome runnables            
      my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome('-genomic'  => $self->unmasked,
								     '-features' => $features,
								     '-seqfetcher' => $self->seqfetcher);

      # run runnable
      $e2g->run;
      
      # sort out output
      my @e2g_genes = $e2g->output;
      
      my @output;
      
        foreach my $e2g_gene (@e2g_genes){
            my @exons = $e2g_gene->sub_SeqFeature;
            foreach my $exon (@exons){
                
                
                my @supp_evidence = $exon->sub_SeqFeature;
                    #display(@sub_Feats);                    
                    foreach my $supp_evidence (@supp_evidence) {
                        
                        $supp_evidence->feature2->source_tag('TEST');
                        $supp_evidence->feature2->primary_tag('TEST');
                        
                        $supp_evidence->feature1->strand($est_strand);
                        $supp_evidence->feature2->strand($est_strand);
                        
                        warn $supp_evidence->feature2->source_tag;
                        #$supp_evidence->source_tag('TEST');
                        push (@{$self->{'_output'}},$supp_evidence);
                    }
            }

        }
      
      #print Dumper($self->output);
      
      
                    
    }
                   
}


=head2 make_miniseq

  Title   : make_miniseq
  Usage   : 
  Function: makes a mini genomic from the genomic sequence and features list
  Returns : 
  Args    : 

=cut




=head2 minimum_introm

  Title   : minimum_intron
  Usage   : 
  Function: Defines minimum intron size for miniseq
  Returns : 
  Args    : 

=cut

sub minimum_intron {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_minimum_intron'} = $arg;
    }

    return $self->{'_minimum_intron'} || 1000;
}

=head2 exon_padding

  Title   : exon_padding
  Usage   : 
  Function: Defines exon padding extent for miniseq
  Returns : 
  Args    : 

=cut
   
sub exon_padding {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_padding'} = $arg;
    }

#    return $self->{'_padding'} || 100;
    return $self->{'_padding'} || 100;


}

sub output {
    my ($self, @arg) = @_;
    
    if (@arg) {
      @{$self->{'_output'}} = @arg;
    }

    if (!defined($self->{'_output'})) {
       $self->{'_output'} = [];
    }
  
    return @{$self->{'_output'}};
  }


1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

