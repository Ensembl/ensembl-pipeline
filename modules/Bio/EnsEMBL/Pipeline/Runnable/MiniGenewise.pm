#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise->new(-genomic  => $genseq,
								  -features => $features)

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

package Bio::EnsEMBL::Pipeline::Runnable::MiniGenewise;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::Genewise;
use Bio::EnsEMBL::Pipeline::MiniSeq;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::DB::RandomAccessI;

#compile time check for executable
use Bio::PrimarySeqI;
use Bio::SeqIO;
use Bio::Root::RootI;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI );

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $genomic, $features, 
	$seqfetcher, $endbias, $forder) = $self->_rearrange([qw(GENOMIC
								FEATURES
								SEQFETCHER
								ENDBIAS
								FORDER)],
							    @args);
    
    $self->throw("No genomic sequence input")                     unless defined($genomic);
    $self->throw("No seqfetcher provided")                        unless defined($seqfetcher);

    $self->throw("[$genomic] is not a Bio::PrimarySeqI")          unless $genomic->isa("Bio::PrimarySeqI");
    $self->throw("[$seqfetcher] is not a Bio::DB::RandomAccessI") unless $seqfetcher->isa("Bio::DB::RandomAccessI");

    $self->genomic_sequence($genomic) if defined($genomic);
    $self->seqfetcher($seqfetcher)    if defined($seqfetcher);
    $self->endbias($endbias)          if defined($endbias);

    $self->{'_forder'} = $forder      if defined($forder);

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
    
    return $self;
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

=head2 endbias

    Title   :   endbias
    Usage   :   $self->endbias($endbias)
    Function:   Get/set method for endbias
    Returns :   
    Args    :   

=cut

sub endbias {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    $self->{'_endbias'} = $arg;
  }
  if (!defined($self->{'_endbias'})) {
    $self->{'_endbias'} = 0;
  }
  return $self->{'_endbias'};
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
  if (defined($value)) {
    #need to check if we are being passed a Bio::DB::RandomAccessI object
    $self->throw("[$value] is not a Bio::DB::RandomAccessI") unless $value->isa("Bio::DB::RandomAccessI");
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

sub make_miniseq {
    my ($self,@features) = @_;

    my $seqname = $features[0]->seqname;
    @features = sort {$a->start <=> $b->start} @features;
    my $count  = 0;
    my $mingap = $self->minimum_intron;
    
    my $pairaln  = new Bio::EnsEMBL::Analysis::PairAlign;

    my @genomic_features;

    my $prevend     = 0;
    my $prevcdnaend = 0;
    
  FEAT: foreach my $f (@features) {

      my $start = $f->start;
      my $end   = $f->end;
      
      $start = $f->start - $self->exon_padding;
      $end   = $f->end   + $self->exon_padding;

      if ($start < 1) { $start = 1;}
      if ($end   > $self->genomic_sequence->length) {$end = $self->genomic_sequence->length;}

      my $gap     =    ($start - $prevend);

      if ($count > 0 && ($gap < $mingap)) {
	# STRANDS!!!!!
	  if ($end < $prevend) { $end = $prevend;}
	    
	  $genomic_features[$#genomic_features]->end($end);
	  $prevend     = $end;
	  $prevcdnaend = $f->hend;
      } else {
	
	    my $newfeature = new Bio::EnsEMBL::SeqFeature;

        $newfeature->seqname ($f->hseqname);
        $newfeature->start     ($start);
	    $newfeature->end       ($end);
	    $newfeature->strand    (1);
# ???	    $newfeature->strand    ($strand);
	    $newfeature->attach_seq($self->genomic_sequence);

	    push(@genomic_features,$newfeature);
	    
	    $prevend = $end;
	    $prevcdnaend = $f->hend; 
	}
	$count++;
    }

    # Now we make the cDNA features
    # but presumably only if we actually HAVE any ... 
    return unless scalar(@genomic_features);

    my $current_coord = 1;
    
    # make a forward strand sequence, but tell genewise to run reversed if the 
    # features are on the reverse strand - handled by is_reversed
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
	
#	$self->print_FeaturePair($fp);

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

sub minimum_intron {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_minimum_intron'} = $arg;
    }
    return $self->{'_minimum_intron'} || 1000;
}

    
sub exon_padding {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_padding'} = $arg;
    }

    return $self->{'_padding'} || 100;

}

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
  
  if (defined($self->{'_seq_cache'}{$id})) {
    return $self->{'_seq_cache'}{$id};
  } 
  
  my $seqfetcher = $self->seqfetcher;    
  my $seq;
  eval {
    $seq = $seqfetcher->get_Seq_by_acc($id);
  };
  
  if ($@) {
    $self->throw("Problem fetching sequence for [$id]: [$@]\n");
  }
  
  if (!defined $seq) {
    $self->throw("Couldn't find sequence for [$id]");
  }
  
  return $seq;
}

=head2 get_all_Sequences

  Title   : get_all_Sequences
  Usage   : my $seq = get_all_Sequences(@id)
  Function: Fetches sequences with ids in @id
  Returns : nothing, but $self->{'_seq_cache'}{$id} 
            has a Bio::PrimarySeq for each $id in @id
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
  Function: Runs genewises
  Returns : none
  Args    : 

=cut

sub run {
    my ($self) = @_;
    

    my @ids = $self->get_all_FeatureIds;

#    $self->get_all_Sequences(@ids);
    my $analysis_obj    = new Bio::EnsEMBL::Analysis
	(-db              => 'genewise',
	 -db_version      => 1,
	 -program         => "genewise",
	 -program_version => 1,
	 -gff_source      => 'genewise',
	 -gff_feature     => 'exon',
	 -logic_name      => 'minigenewise');

    foreach my $id (@ids) {
	my $hseq = $self->get_Sequence(($id));

	if (!defined($hseq)) {
	    $self->throw("Can't fetch sequence for id [$id]\n");
	}

	
	my $eg = new Bio::EnsEMBL::Pipeline::Runnable::Genewise(-genomic => $self->genomic_sequence,
								-protein => $hseq,
								-memory  => 400000);

	$eg->run;

	my @f = $eg->output;

	foreach my $f (@f) {
	    #print("Aligned output is " . $id . "\t" . $f->start . "\t" . $f->end . "\t" . $f->score . "\n");
	    print $f;
	}

	push(@{$self->{'_output'}},@f);

    }
}

=head2 minirun

  Title   : minirun
  Usage   : $self->minirun()
  Function: Runs genewise on MiniSeq representation of genomic sequence
  Returns : none
  Args    : 

=cut

sub minirun {
  my ($self) = @_;
  
  my ($idhash) = $self->get_all_FeaturesById;
  
  my @ids    = keys %$idhash;
  
  if (defined($self->{'_forder'})) {
    @ids = @{$self->{'_forder'}};
  }
 
  $self->get_all_Sequences(@ids);

  my $analysis_obj    = new Bio::EnsEMBL::Analysis
    (-db              => undef,
     -db_version      => undef,
     -program         => "genewise",
     -program_version => 1,
     -gff_source      => 'genewise',
     -gff_feature     => 'similarity');
  
 ID: foreach my $id (@ids) {
    
    my $features = $idhash->{$id};
    my @exons;
    
    next ID unless (ref($features) eq "ARRAY");
    
    # why > not >= 1?
    next ID unless (scalar(@$features) >= 1);
    
    # forward and reverse split.
    my @forward;
    my @reverse;

    # forward/reverse strand split
    foreach my $feat(@$features) {
      if($feat->strand == 1) { push(@forward,$feat); }
      elsif($feat->strand == -1) { push(@reverse,$feat); }
      else { $self->throw("unstranded feature not much use for gene building\n") }
    }
    
    # run on each strand
    eval {
      $self->run_blastwise($id, \@forward, $analysis_obj);
    };
    if ($@) {
      print STDERR "Error running blastwise for forward strand on " . $features->[0]->hseqname . " [$@]\n";
    }

    eval {
      $self->run_blastwise($id, \@reverse, $analysis_obj);
    };
    if ($@) {
      print STDERR "Error running blastwise for reverse strand on " . $features->[0]->hseqname . " [$@]\n";
    }

  }
  
}

=head2 run_blastwise

  Title   : run_blastwise
  Usage   : $self->run_blastwise()
  Function: Runs genewise on a MiniSeq
  Returns : none
  Args    : 

=cut

sub run_blastwise {
  my ($self,$id,$features,$analysis_obj) = @_;

  my @extras  = $self->find_extras (@$features);
  
#  print STDERR "Number of extra features = " . scalar(@extras) . "\n";
  
  return unless (scalar(@extras) >= 1);
  
  my $miniseq = $self->make_miniseq(@$features);
  my $hseq    = $self->get_Sequence($id);
  
  my $reverse = $self->is_reversed(@$features);
  
  my $endbias = $self->endbias;
     
  if (!defined($hseq)) {
    $self->throw("Can't fetch sequence for id [$id]\n");
  }
  
  my $eg = new Bio::EnsEMBL::Pipeline::Runnable::Genewise(  -genomic => $miniseq->get_cDNA_sequence,
							    -protein => $hseq,
							    -memory  => 400000,
							    "-reverse" => $reverse,
							    "-endbias"  => $endbias);
  
  $eg->run;

  # output is a list of Features (one per exon) with subseqfeatures 
  # representing the ungapped sub alignments for each exon
  my @f = $eg->output;

#  print STDERR "number of exons: " . scalar(@f)  . "\n";
  
  my @newf;

  my $strand = 1;
  if ($reverse == 1) {
    $strand = -1;
  }

  my $ec = 0;
 FEAT: foreach my $f (@f) {
    $ec++;
    $f->strand($strand); # genomic
    my $phase = $f->phase();
    
    # need to convert whole exon back to genomic coordinates
    my @genomics = $miniseq->convert_SeqFeature($f);         
    
    if ($#genomics > 0) {
      # all hell will break loose as the sub alignments will probably not map cheerfully 
      # and we may start introducing in frame stops ...
      # for now, ignore this feature.
      print STDERR "Warning : feature converts into > 1 features " . scalar(@genomics) . " ignoring exon $ec\n";
      next FEAT;
    }
    
    
    # also need to convert each of the sub alignments back to genomic coordinates
    foreach my $aln($f->sub_SeqFeature) {
      my @alns = $miniseq->convert_PepFeaturePair($aln);
      
      if ($#alns > 0) {
	# we're in for fun
	print STDERR "Warning : sub_align feature converts into > 1 features " . scalar(@alns) . "\n";
      }
      foreach my $a(@alns) {
	my $added = 0;
	$a->strand($strand); # genomic
	$a->hstrand(1);      # protein
	
	$a->seqname($aln->seqname);
	$a->hseqname($aln->hseqname);
	# shouldn't need to expand ... as long as we choose the right parent feature to add to!
	foreach my $g(@genomics){
	  if($a->start >= $g->start && $a->end <=$g->end && !$added){
	    $g->add_sub_SeqFeature($a,'');
	    $added = 1;
	  }
	}
	$self->warn("Sub align feature could not be added ...\n") unless $added;
      }
    }

    push(@newf,@genomics);
    
    foreach my $nf (@genomics) {
      $nf->phase($phase);
      $nf->strand($strand);
      #BUGFIX: This should probably be fixed in Bio::EnsEMBL::Analysis
      $nf->seqname($f->seqname);
      $nf->score   (100);
      $nf->analysis($analysis_obj);
      #end BUGFIX
    }
    
  }
 
  # $fset holds a list of (genomic) SeqFeatures (one fset per gene) plus their constituent exons and
  # sub_SeqFeatures representing ungapped alignments making up the exon:protein alignment
  my $fset = new Bio::EnsEMBL::SeqFeature();
  
  foreach my $nf (@newf) {
    $fset->add_sub_SeqFeature($nf,'EXPAND');
    $fset->seqname($nf->seqname);
    $fset->analysis($analysis_obj);
  }
  
  push(@{$self->{'_output'}},$fset);
  
}

sub is_reversed {
    my ($self,@features) = @_;

    my $strand = 0;

    my $fcount = 0;
    my $rcount = 0;

    foreach my $f (@features) {
	if ($f->strand == 1) {
	    $fcount++;
	} elsif ($f->strand == -1) {
	    $rcount++;
	}
    }

    if ($fcount > $rcount) {
	return 0;
    } else {
	return 1;
    }
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
    my ($self, $f1score, $f1start, $f1end, $f1id, 
	$f2start, $f2end, $f2id,
        $f1source, $f2source, $f1strand, $f2strand, 
	$f1primary, $f2primary) = @_;
    
    #create analysis object
    my $analysis_obj    = new Bio::EnsEMBL::Analysis
                                (-db              => 'genewise',
                                 -db_version      => 1,
                                 -program         => "genewise",
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

