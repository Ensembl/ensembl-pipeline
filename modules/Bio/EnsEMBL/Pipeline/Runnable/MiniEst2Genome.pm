#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome->new('-genomic'    => $genseq,
								    '-features'   => $features,
								    '-seqfetcher' => $seqfetcher
								   )

    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;
use Bio::EnsEMBL::Pipeline::MiniSeq;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::DB::RandomAccessI;

use Bio::PrimarySeqI;
use Bio::SeqIO;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI );

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_fplist'} = []; #create key to an array of feature pairs
  
  my( $genomic, $features, $seqfetcher ) = $self->_rearrange([qw(GENOMIC
								 FEATURES
								 SEQFETCHER)], @args);
  
  $self->throw("No genomic sequence input")           
    unless defined($genomic);
  $self->throw("[$genomic] is not a Bio::PrimarySeqI") 
    unless $genomic->isa("Bio::PrimarySeqI");
  $self->genomic_sequence($genomic) if defined($genomic);

  $self->throw("No seqfetcher provided")           
    unless defined($seqfetcher);
  $self->throw("[$seqfetcher] is not a Bio::DB::RandomAccessI") 
    unless $seqfetcher->isa("Bio::DB::RandomAccessI");
  $self->seqfetcher($seqfetcher) if defined($seqfetcher);
  
  
  if (defined($features)) {
    if (ref($features) eq "ARRAY") {
      my @f = @$features;
      
      foreach my $f (@f) {
	$self->addFeature($f);
      }
    } else {
      $self->throw("[$features] is not an array ref.");
    }
  }
  
  return $self; # success - we hope!
}

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::DB::RandomAccessI object
    Args    :   Bio::DB::RandomAccessI object

=cut

sub seqfetcher {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::DB::RandomAccessI object
        $value->isa("Bio::DB::RandomAccessI") || $self->throw("Input isn't a Bio::DB::RandomAccessI");
        $self->{'_seqfetcher'} = $value;
    }
    return $self->{'_seqfetcher'};
}

=head2 addFeature 

    Title   :   addFeature
    Usage   :   $self->addFeature($f)
    Function:   Adds a feature to the object for realigning
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   Bio::EnsEMBL::FeaturePair

=cut

sub addFeature {
    my( $self, $value ) = @_;
    
    if(!defined($self->{'_features'})) {
	$self->{'_features'} = [];
    }

    if ($value) {
        $value->isa("Bio::EnsEMBL::FeaturePair") || $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair");
	push(@{$self->{'_features'}},$value);
    }
}


=head2 get_all_FeaturesbyId

    Title   :   get_all_FeaturesById
    Usage   :   $hash = $self->get_all_FeaturesById;
    Function:   Returns a ref to a hash of features.
                The keys to the hash are distinct feature ids
    Returns :   ref to hash of Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut

sub get_all_FeaturesById {
    my( $self) = @_;
    
    my  %idhash;

    FEAT: foreach my $f ($self->get_all_Features) {
#	print STDERR ("Feature is $f " . $f->seqname . "\t" . $f->hseqname ."\n");
    if (!(defined($f->hseqname))) {
	$self->warn("No hit name for " . $f->seqname . "\n");
	    next FEAT;
	} 
	if (defined($idhash{$f->hseqname})) {
	    push(@{$idhash{$f->hseqname}},$f);
	} else {
	    $idhash{$f->hseqname} = [];
	    push(@{$idhash{$f->hseqname}},$f);
	}

    }

    return (\%idhash);
}


=head2 get_all_Features

    Title   :   get_all_Features
    Usage   :   @f = $self->get_all_Features;
    Function:   Returns the array of features
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   none

=cut


sub get_all_Features {
    my( $self, $value ) = @_;
    
    return (@{$self->{'_features'}});
}


=head2 get_all_FeatureIds

  Title   : get_all_FeatureIds
  Usage   : my @ids = get_all_FeatureIds
  Function: Returns an array of all distinct feature hids 
  Returns : @string
  Args    : none

=cut

sub get_all_FeatureIds {
    my ($self) = @_;

    my %idhash;

    foreach my $f ($self->get_all_Features) {
	if (defined($f->hseqname)) {
	    $idhash{$f->hseqname} = 1;
	} else {
	    $self->warn("No sequence name defined for feature. " . $f->seqname . "\n");
	}
    }

    return keys %idhash;
}

=head2 make_miniseq

  Title   : make_miniseq
  Usage   : 
  Function: makes a mini genomic from the genomic sequence and features list
  Returns : 
  Args    : 

=cut

sub make_miniseq {
    my ($self,@features) = @_;

    my $seqname = $features[0]->seqname;

#    print STDERR "seqname $seqname\n";

    @features = sort {$a->start <=> $b->start} @features;
    my $count  = 0;
    my $mingap = $self->minimum_intron;

    my $pairaln = new Bio::EnsEMBL::Analysis::PairAlign;

    my @genomic_features;

    my $prevend     = 0;
    my $prevcdnaend = 0;
    
  FEAT: foreach my $f (@features) {
#      print STDERR "Found feature - " . $f->hseqname . "\t" . $f->start . "\t" . $f->end . "\t" . $f->strand . "\n"; 

      my $start = $f->start;
      my $end   = $f->end;
      
      $start = $f->start - $self->exon_padding;
      $end   = $f->end   + $self->exon_padding;

      if ($start < 1) { $start = 1;}
      if ($end   > $self->genomic_sequence->length) {$end = $self->genomic_sequence->length;}

      my $gap     =    ($start - $prevend);

#      print STDERR "Feature hstart is " . $f->hstart . "\t" . $prevcdnaend . "\n";
#      print STDERR "Padding feature - new start end are $start $end\n";

#      print STDERR "Count is $count : $mingap " . $gap  . "\n";

      if ($count > 0 && ($gap < $mingap)) {
	# STRANDS!!!!!
	  if ($end < $prevend) { $end = $prevend;}
#	  print(STDERR "Merging exons in " . $f->hseqname . " - resetting end to $end\n");
	    
	  $genomic_features[$#genomic_features]->end($end);
	  $prevend     = $end;
	  $prevcdnaend = $f->hend;
#	  print STDERR "Merged start end are " . $genomic_features[$#genomic_features]->start . "\t" .  $genomic_features[$#genomic_features]->end . "\n";
      } else {
	
	    my $newfeature = new Bio::EnsEMBL::SeqFeature;

        $newfeature->seqname ($f->hseqname);
        $newfeature->start     ($start);
	    $newfeature->end       ($end);
	    $newfeature->strand    (1);
# ???	    $newfeature->strand    ($strand);
	    $newfeature->attach_seq($self->genomic_sequence);

	    push(@genomic_features,$newfeature);
	    
#	    print(STDERR "Added feature $count: " . $newfeature->start  . "\t"  . 
#		  $newfeature->end    . "\t " . 
#		  $newfeature->strand . "\n");

	    $prevend = $end;
	    $prevcdnaend = $f->hend; 
#	    print STDERR "New end is " . $f->hend . "\n";

	}
	$count++;
    }

    # Now we make the cDNA features
    # but presumably only if we actually HAVE any ... 
    return unless scalar(@genomic_features);

    my $current_coord = 1;

    # make a forward strand sequence, est2genome runs -reverse 
    @genomic_features = sort {$a->start <=> $b->start } @genomic_features;

    foreach my $f (@genomic_features) {
	$f->strand(1);
	my $cdna_start = $current_coord;
	my $cdna_end   = $current_coord + ($f->end - $f->start);
	
	my $tmp = new Bio::EnsEMBL::SeqFeature(
					       -seqname => $f->seqname.'.cDNA',
					       -start => $cdna_start,
					       -end   => $cdna_end,
					       -strand => 1);
	
	my $fp  = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
						-feature2 => $tmp);
	
	$pairaln->addFeaturePair($fp);
	
	$self->print_FeaturePair($fp);

	$current_coord = $cdna_end+1;
    }
	
    #changed id from 'Genomic' to seqname
    my $miniseq = new Bio::EnsEMBL::Pipeline::MiniSeq(-id        => $seqname,
						      -pairalign => $pairaln);

    my $newgenomic = $miniseq->get_cDNA_sequence->seq;
    $newgenomic =~ s/(.{72})/$1\n/g;
#    print ("New genomic sequence is " . $newgenomic. "\n");
    return $miniseq;

}


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
    return $self->{'_padding'} || 1000;

}

=head2 print_FeaturePair

  Title   : print_FeaturePair
  Usage   : 
  Function: for debugging
  Returns : 
  Args    : 

=cut

sub print_FeaturePair {
    my ($self,$nf) = @_;
    #changed $nf->id to $nf->seqname
    print(STDERR "FeaturePair is " . $nf->seqname    . "\t" . 
	  $nf->start . "\t" . 
	  $nf->end   . "\t(" . 
	  $nf->strand . ")\t" .
	  $nf->hseqname  . "\t" . 
	  $nf->hstart   . "\t" . 
	  $nf->hend     . "\t(" .
	  $nf->hstrand  . ")\n");
}


=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id)
  Function: Fetches sequences with id $id
  Returns : Bio::PrimarySeq
  Args    : none

=cut

sub get_Sequence {
    my ($self,$id) = @_;
    my $seqfetcher = $self->seqfetcher;

    if (defined($self->{'_seq_cache'}{$id})) {
      return $self->{'_seq_cache'}{$id};
    } 
    
    my $seq;
    eval {
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };

    if ($@) {
      $self->throw("Problem fetching sequence for [$id]: [$@]\n");
    }

    if (!defined($seq)) {
      $self->throw("Couldn't find sequence for [$id]");
    }
    
#    print (STDERR "Found sequence for $id [" . $seq->length() . "]\n");
    
    return $seq;

}

=head2 get_all_Sequences

  Title   : get_all_Sequences
  Usage   : my $seq = get_all_Sequences(@id)
  Function: Fetches sequences with ids in @id
  Returns : nothing, but $self->{'_seq_cache'}{$id} has a Bio::PrimarySeq for each $id in @id
  Args    : array of ids

=cut

sub get_all_Sequences {
  my ($self,@id) = @_;
  
 SEQ: foreach my $id (@id) {
    my $seq = $self->get_Sequence($id);
    if(defined $seq) {
      $self->{'_seq_cache'}{$id} = $seq;
    }
  }
}

=head2 run

  Title   : run
  Usage   : $self->run()
  Function: Runs est2genome on MiniSeq representation of genomic sequence for each EST
  Returns : none
  Args    : 

=cut

sub run {
  my ($self) = @_;
  
  my ($esthash) = $self->get_all_FeaturesById;
  
  my @ests    = keys %$esthash;
  
  $self->get_all_Sequences(@ests);

  my $analysis_obj    = new Bio::EnsEMBL::Analysis
    (-db              => undef,
     -db_version      => undef,
     -program         => "est2genome",
     -program_version => 1,
     -gff_source      => 'est2genome',
     -gff_feature     => 'similarity',);
  
 ID: foreach my $est (@ests) {
    
    my $features = $esthash->{$est};
    my @exons;
    
#    print(STDERR "Processing $est\n");
    next ID unless (ref($features) eq "ARRAY");
    
#    print(STDERR "Features = " . scalar(@$features) . "\n");

    # why > not >= 1?
    next ID unless (scalar(@$features) >= 1);
    
    eval {
      $self->run_blaste2g($est, $features, $analysis_obj);
    };

    if ($@) {
      print STDERR "Error running blaste2g on" . $features->[0]->hseqname . " [$@]\n";
    }

  }
  
}

=head2 run_blaste2g

  Title   : run_blaste2g
  Usage   : $self->run_blaste2g()
  Function: Runs est2genome on a MiniSeq
  Returns : none
  Args    : 

=cut

sub run_blaste2g {
  my ($self,$est,$features,$analysis_obj) = @_;
  
  #?? never did fully understand this.
  my @extras  = $self->find_extras (@$features);
#  print STDERR "Number of extra features = " . scalar(@extras) . "\n";
  return unless (scalar(@extras) >= 1);
  
  my $miniseq = $self->make_miniseq(@$features);
  my $hseq    = $self->get_Sequence($est);
  
  if (!defined($hseq)) {
    $self->throw("Can't fetch sequence for id [$est]\n");
  }
  
  my $eg = new Bio::EnsEMBL::Pipeline::Runnable::Est2Genome(  -genomic => $miniseq->get_cDNA_sequence,
							      -est     => $hseq);
  
  $eg->run;
  
  # output is a list of Features, one per predicted gene. Exons are added as 
  # subseqfeatures of gene features, and supporting evidence featurepairs as 
  # subseqfeatures of exons.
  my @genes = $eg->output;
  
#  print STDERR "number of genes: " . scalar(@genes)  . "\n";
  
  my @newf;
  
  if ( scalar(@genes) >1 ) {
    $self->throw("more than one gene predicted - I'm outta here!\n");
  }
  
  my @genomic_exons;
  my $ec = 0;
  my $strand;

  foreach my $gene(@genes) {
    my @exons = $gene->feature1->sub_SeqFeature;
    print "*numexons: " . scalar(@exons) . "\n";

    $strand = $exons[0]->strand;
  FEAT:    foreach my $ex(@exons){
      $ec++;
      $self->throw("mismatched exon strands\n") unless $ex->strand == $strand;
      
      # exonerate has no concept of phase, but remapping will fail if this is unset
      $ex->phase(0);
      
      # convert back to genomic coords, but leave the EST coordinates alone
      my @converted = $miniseq->convert_FeaturePair($ex);
      if ($#converted > 0) {
	# all hell will break loose as the sub alignments will probably not map cheerfully 
	# for now, ignore this feature.
	print STDERR "Warning : feature converts into > 1 features " . scalar(@converted) . " ignoring exon $ec\n";
	next FEAT;
      }
      
      foreach my $nf(@converted) {
	push(@genomic_exons, $nf->feature1); 
      }
      
      # now sort out sub seqfeatures - details of sub segments making up an exon.
#      print STDERR "minie2g: " . $ex->feature1->sub_SeqFeature .  "\n";
#      print STDERR "minie2g num subf: " . scalar($ex->feature1->sub_SeqFeature) . "\n";

      foreach my $aln($ex->feature1->sub_SeqFeature){
	# strands 
	if($aln->strand != $aln->hstrand){
	  $aln->strand(-1);
	  $aln->hstrand(1);
	  $ex->strand(-1);
	}
	
	# convert to genomic coords
	my @alns = $miniseq->convert_FeaturePair($aln);
	if ($#alns > 0) {
	  # we're in for fun
	  print STDERR "Warning : sub_align feature converts into > 1 features " . scalar(@alns) . "\n";
	}
	
	foreach my $a(@alns) {
	  my $added = 0;
	  $a->strand($aln->strand); # genomic
	  $a->hstrand(1);      # protein
	  
	  $a->seqname($aln->seqname);
	  $a->hseqname($aln->hseqname);
	  # shouldn't need to expand ... as long as we choose the right parent feature to add to!
	  foreach my $g(@genomic_exons){
	    if($a->start >= $g->start && $a->end <=$g->end && !$added){
	      $g->add_sub_SeqFeature($a,'');
	      $added = 1;
	    }
	  }
	  $self->warn("Sub align feature could not be added ...\n") unless $added;
	}
      }
    }
    
    
    foreach my $gex (@genomic_exons) {
      $gex->{_phase} = 0; # e2g doesn;t give us phase
      $gex->strand($strand);
      #BUGFIX: This should probably be fixed in Bio::EnsEMBL::Analysis
      $gex->seqname($gene->seqname); # urrmmmm?
      $gex->score   (100); # probably not
      $gex->analysis($analysis_obj);
      #end BUGFIX
    }
  }   
  
  
  # $fset holds a list of (genomic) SeqFeatures (one fset per gene) plus their constituent exons and
  # sub_SeqFeatures representing ungapped alignments making up the exon:EST alignment
  my $fset = new Bio::EnsEMBL::SeqFeature();
  
  foreach my $nf (@genomic_exons) {
    $fset->add_sub_SeqFeature($nf,'EXPAND');
    $fset->seqname($nf->seqname);
    $fset->analysis($analysis_obj);
    #   $nf->strand($nf->hstrand);
    #    print(STDERR "Realigned output is " . $nf->seqname    . "\t" . 
    #	  $nf->start     . "\t" . 
    #	  $nf->end       . "\t(" . 
    #	  $nf->strand    . ")\t" .
    #	  $nf->hseqname  . "\t" . 
    #	  $nf->hstart    . "\t" . 
    #	  $nf->hend      . "\t(" .
    #	  $nf->hstrand   . ")\t:" .
    #	  $nf->feature1->{_phase} . ":\t:" . 
    #	  $nf->feature2->{_phase} . ":\n");
  }
  
  push(@{$self->{'_output'}},$fset);
  
}

sub find_extras {
  my ($self,@features) = @_;
  
  my @output = $self->output;
  my @new;
  
 FEAT: foreach my $f (@features) {
    my $found = 0;
    if (($f->end - $f->start) < 50) {
      next FEAT;
    }
    #	print ("New feature\n");
    
    #$self->print_FeaturePair($f);
    foreach my $out (@output) {
      foreach my $sf ($out->sub_SeqFeature) {
	
	if (!($f->end < $out->start || $f->start >$out->end)) {
	  $found = 1;
	}
      }
    }
    
    if ($found == 0) {
      push(@new,$f);
    }
  }
  return @new;
}

=head2 output

  Title   : output
  Usage   : $self->output
  Function: Returns results of est2genome as array of FeaturePair
  Returns : An array of Bio::EnsEMBL::FeaturePair
  Args    : none

=cut

sub output {
    my ($self) = @_;
    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }
    return @{$self->{'_output'}};
}

sub _createfeatures {
    my ($self, $f1score, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id,
        $f1source, $f2source, $f1strand, $f2strand, $f1primary, $f2primary) = @_;
    
    #create analysis object
    my $analysis_obj    = new Bio::EnsEMBL::Analysis
                                (-db              => 'est2genome',
                                 -db_version      => 1,
                                 -program         => "est2genome",
                                 -program_version => 1,
                                 -gff_source      => $f1source,
                                 -gff_feature     => $f1primary,);
    
    #create features
    my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start =>  $f1start,
                                              -end =>     $f1end,
                                              -seqname =>      $f1id,
                                              -strand =>  $f1strand,
                                              -score =>   $f1score,
                                              -source =>  $f1source,
                                              -primary => $f1primary,
                                              -analysis => $analysis_obj );
 
     my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start =>  $f2start,
                                                -end =>    $f2end,
                                                -seqname =>$f2id,
                                                -strand => $f2strand,
                                                -score =>  undef,
                                                -source => $f2source,
                                                -primary =>$f2primary,
                                                -analysis => $analysis_obj );
    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
                                             -feature2 => $feat2) ;
 
    $self->_growfplist($fp); 
}

sub _growfplist {
    my ($self, $fp) =@_;
    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
}

sub _createfiles {
    my ($self, $genfile, $estfile, $dirname)= @_;
    
    #check for diskspace
    my $spacelimit = 0.1; # 0.1Gb or about 100 MB
    my $dir ="./";
    unless ($self->_diskspace($dir, $spacelimit)) 
    {
        $self->throw("Not enough disk space ($spacelimit Gb required)");
    }
            
    #if names not provided create unique names based on process ID    
    $genfile = $self->_getname("genfile") unless ($genfile);
    $estfile = $self->_getname("estfile") unless ($estfile);    
    #create tmp directory    
    mkdir ($dirname, 0777) or $self->throw ("Cannot make directory '$dirname' ($?)");
    chdir ($dirname) or $self->throw ("Cannot change to directory '$dirname' ($?)"); 
    return ($genfile, $estfile);
}
    

sub _getname {
    my ($self, $typename) = @_;
    return  $typename."_".$$.".fn"; 
}

sub _diskspace {
    my ($self, $dir, $limit) =@_;
    my $block_size; #could be used where block size != 512 ?
    my $Gb = 1024 ** 3;
    
    open DF, "df $dir |" or $self->throw ("Can't open 'du' pipe");
    while (<DF>) 
    {
        if ($block_size) 
        {
            my @L = split;
            my $space_in_Gb = $L[3] * 512 / $Gb;
            return 0 if ($space_in_Gb < $limit);
            return 1;
        } 
        else 
        {
            ($block_size) = /(\d+).+blocks/i
                || $self->throw ("Can't determine block size from:\n$_");
        }
    }
    close DF || $self->throw("Error from 'df' : $!");
}


sub _deletefiles {
    my ($self, $genfile, $estfile, $dirname) = @_;
    unlink ("$genfile") or $self->throw("Cannot remove $genfile ($?)\n");
    unlink ("$estfile") or $self->throw("Cannot remove $estfile ($?)\n");
    chdir ("../");
    rmdir ($dirname) or $self->throw("Cannot remove $dirname \n");
}

1;
