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

Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep

=head1 SYNOPSIS

  my $db       = Bio::EnsEMBL::DBLoader->new($locator);
  my $genscan  = Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanDNA->new(
                -dbobj      => $db,
                -input_id   => $input_id
                -analysis   => $analysis);
  $genscan->fetch_input();
  $genscan->run();
  $genscan->output();
  $genscan->write_output(); #writes to DB

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides
constructed from assembling genscan predicted features to peptide
sequence. The resulting blast hits are written back as FeaturePairs.
The appropriate Bio::EnsEMBL::Analysis object must be passed
for extraction of appropriate parameters. A
Bio::EnsEMBL::Pipeline::DBSQL::Obj is required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanDNA;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-genomic     => $genomic
                           -peptide     => $peptide
                           -program     => 'tblastn',
			   -database    => 'dbest',
			   -threshold   => 1e-6);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanDNA object
    Returns :   A Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanDNA object
    Args    :   -genomic   : Bio::Seq containing genomic dna
                -pep       : Bio::EnsEMBL::Transcript
                -program   : the flavour of blast to run
                -database  : the database to run against
                -threshold : the probability to filter the features by



=cut


sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_featurepairs'}= [];
  
  $self->{'_peptide'}     = undef;
  $self->{'_genomic'}     = undef;
  $self->{'_program'}     = undef;
  $self->{'_database'}    = undef;
  $self->{'_threshold'}   = undef;
  $self->{'_threshold_type'}   = undef;
  $self->{'_options'}     = undef;
    
  
  # Read the input parameters and set them in the object

  my ( $genomic,$peptide,$program,$database,$threshold,$threshold_type,$options) = 
    $self->_rearrange ([qw(GENOMIC 
			   PEPTIDE 
			   PROGRAM 
			   DATABASE 
			   THRESHOLD 
			   THRESHOLD_TYPE
			   OPTIONS)], @args);
  
  
  if (defined($genomic) && $genomic->isa("Bio::PrimarySeqI")) {
      $self->genomic($genomic);
  } elsif (defined($genomic)) {
      $self->throw("[$genomic] is not a Bio::PrimarySeqI");
  } else {
      $self->throw("No genomic sequence input");
  }
  
  if (defined($peptide) && $peptide->isa("Bio::EnsEMBL::PredictionTranscript")) {
      $self->peptide($peptide);
  } elsif (defined($peptide)) {
      $self->throw("[$peptide] is not a Bio::EnsEMBL::PredictionTranscript");
  } else {
    $self->throw("No peptide input");
  }
  
  if (defined($program)) {
    $self->program($self->find_executable($program));
  } else {
    $self->throw ("No program input");
  }
  
  if (defined($database)) {
    $self->database($database);
  } else {
    $self->throw("No database defined");
  }

  if (defined($threshold)) {
    $self->threshold($threshold);
  } else {
    $self->threshold(0);
  }

  if (defined($threshold_type)) {
    $self->threshold_type($threshold_type);
  } else {
    $self->threshold_type('PVALUE');
  }
  
  if (defined($options)) {
    $self->options($options);
  } 
  
  return $self;
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;

    my $transcript = $self->peptide;
    
    if (!defined($transcript)) {
      $self->throw("No peptide input");
    }

    print STDERR "Creating BioPrimarySeq ". " " . $transcript->translate() . "\n";

    my $peptide = Bio::PrimarySeq->new(-id         => 'Genscan_prediction',
				       -seq        => $transcript->translate(),
				       -moltype    => 'protein' );

    print STDERR "Peptide length: ", $peptide->length, "\n";
    if ($peptide->length < 3) {
      print "Peptide too short (min length is 3); skipping transcript\n";
      return;
    }

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Blast  (-query     => $peptide,
								 -program   => $self->program,
								 -database  => $self->database,
								 -threshold => $self->threshold,
					
								 -threshold_type => $self->threshold_type,
								 -ungapped  => 0,
								 -options   => $self->options,
                                                                 -filter    => 1);

   $runnable->run();
  
    $self->align_hits_to_contig2($runnable->output);
    #$self->check_features($transcript->translate->seq,$self->featurepairs);
  }

sub check_features {
  my ($self,$pep,@f) = @_;

  print STDERR "Peptide is " . $pep . "\n";

  my %seqhash;

  foreach my $f (@f) {
      eval {
	  if (!defined($seqhash{$f->hseqname})) {
	      my $seq = $self->get_Sequence($f->hseqname);
	      $seqhash{$f->hseqname} = $seq;
	  }
	  my $seq = $seqhash{$f->hseqname};

	  my $fdna = $f->seq;
	  my $rdna = $f->seq->revcom;

	  $rdna = $rdna->seq;
	  $fdna = $fdna->seq;
	  
	  $fdna =~ tr/a-z/A-Z/;
	  $rdna =~ tr/a-z/A-Z/;

	  my $hdna = substr($seq->seq,$f->hstart-1,($f->hend - $f->hstart + 1));
	  
	  
	  $hdna =~ tr/a-z/A-Z/;
	  
	  print "\tFeature " . $f->start . "\t" . $f->end . "\t" . $f->strand . "\t" . $f->phase . "\t" . $f->hstart . "\t" . $f->hend . " "  . "\n" ;
	  print $fdna . "\n$rdna\n" . $hdna . "\n";
      };
      if ($@) {
	  print STDERR "Couldn't fetch sequence for " . $f->hseqname . " No alignment printed [$@]\n";
      }
  }
}

sub get_Sequence {
    my ($self,$id) = @_;

    if ($id =~ /\|/) {
        # Take unigene id or accession
        my ($ug) = $id =~ m{/ug=(.*?)\ };
        if ($id =~ /\|UG\|/) {
           if (length $ug > 0) {
              $id = $ug;
           }
           else {
              $id =~ s/.*\|.*\|(.*)/$1/;
           }
        } else {
	   $id =~ s/.*\|(.*)\|.*/$1/;
        }
    }
    $id =~ s/\..*//;

    next ID unless defined($id);

    print(STDERR "Sequence id :  is [$id]\n");

    open(IN,"efetch -q $id |") || $self->throw("Error fetching sequence for id [$id]");
	
    my $seqstr;
	
    while (<IN>) {
	chomp;
	$seqstr .= $_;
    }
    
    

    if (!defined($seqstr) || $seqstr eq "no match") {
	print("Couldn't find sequence for [$id]");
	return;
    }

    my $seq = new Bio::Seq(-id  => $id,
			   -seq => $seqstr);
    

    print (STDERR "Found sequence for $id [" . $seq->length() . "]\n");

    return $seq;
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    return $self->featurepairs();  
}

sub align_hits_to_contig2 {
  my ( $self, @features )  = @_;
  
  # for each feature
  
  for my $feature ( @features ) {
    print ::LOG join
      ( "\n", 
	( "\n", "Blast result:",
	  "Start ".$feature->start." End ".$feature->end,
	  "hstart ".$feature->hstart." hend ".$feature->hend,
	  "qury: ".$feature->{'qseq'},
	  "subj: ".$feature->{'sseq'},
	  "\n" ));


    my %exon_hash = ();
  # for each ungapped piece in it
    for my $ugFeature ( $feature->ungapped_features() ) {
      
      my $f = $ugFeature;
      print ::LOG "ugfeature: ",join( " ", ( $f->start(), $f->end(), $f->strand(), "-", $f->hstart(), $f->hend(), $f->hstrand() )),"\n";      # ask the $self->peptide to do cdna2genomic
      #   for f->start*3-2 f->end*3
      my @split = $self->peptide->cdna2genomic
	(( $ugFeature->start() * 3 -2 ), 
	 ( $ugFeature->end() * 3 ));
      
      for my $gcoord ( @split ) {
	
	my ( $gstart, $gend, $gstrand, $contig, $exon, $cdna_start, $cdna_end ) =
	  @$gcoord;

	# Take the pieces and make a list of features from it
	# hash them by exon

	# first, eat away non complete codons from start
	while(( $cdna_start - 1 ) % 3 != 0 ) {
	  $cdna_start++;
	  if( $gstrand == 1 ) {
	    $gstart++;
	  } else {
	    $gend--;
	  }
	}
	  
	# and from end
	while( $cdna_end  % 3 != 0 ) {
	  $cdna_end--;
	  if( $gstrand == 1 ) {
	    $gend--;
	  } else {
	    $gstart++;
	  }
	}

	if( $cdna_end <= $cdna_start ) {
	  next;
	}

	my $dna_feat1 = Bio::EnsEMBL::SeqFeature->new 
	  ( -seqname    =>  $contig->name,
	    -start      =>  $gstart, 
	    -end        =>  $gend,
	    -strand     =>  $gstrand,
	    -score      =>  $feature->score,
	    -p_value    =>  $feature->p_value,
	    -percent_id =>  $feature->percent_id,
	    -analysis   =>  $feature->analysis );

	
	my $dna_feat2 = Bio::EnsEMBL::SeqFeature->new 
	  ( -seqname    =>  $feature->hseqname,
	    -start      =>  $cdna_start - ($ugFeature->start()*3-2) + $ugFeature->hstart(),
	    -end        =>  $cdna_end - ($ugFeature->start()*3-2) + $ugFeature->hstart(),
	    -score      =>  $feature->score,
	    -p_value    =>  $feature->p_value,
	    -percent_id =>  $feature->percent_id,
	    -analysis   =>  $feature->analysis );

      
	my $featurepair = Bio::EnsEMBL::FeaturePair->new (-feature1   => $dna_feat1,
							  -feature2   => $dna_feat2 );
	
	push( @{$exon_hash{$exon}}, $featurepair );
      }
    }

    # Take the pieces for each exon and make gapped feature
    foreach my $ex ( keys %exon_hash ) {
      my $dna_align_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new
	(-features => $exon_hash{$ex});
      $self->featurepairs($dna_align_feature);
    }
  }
}



# This function creates a hash which is used to map between the exon genomic position
# and a position within the genscan predicted peptide. The hash is then matched
# against blast peptide hits to return a set of featurepairs of exons and blast
# peptides
sub align_hits_to_contig {
    my ($self, @features) = @_;

    my (%dna_align, @exon_aligns, @featurepairs); #structure for tracking alignment variables
    
    $dna_align {'exons'} = [];

    my $trans = $self->peptide;

    my $pep = $trans->translate();
    my $stop_codon_present = 0;

    #calculate boundaries and map exons to translated peptide
    #Note: Exons have an extra 3 bases for the stop codon. Peptides lack this
    foreach my $exon ($trans->get_all_Exons) {

        my %ex_align;

	my $strand = "+";
	if ($exon->strand == -1) {
	  $strand = "-";
	}

        my $expep = $exon->translate->seq();
	if( $expep =~ /\*$/ ) {
	  $expep =~ s/\*$//g;
	  $stop_codon_present = 1;
	}
	
#        my ($expep) = ($exon->translate->seq =~ /[^\*]+/g);
	if ($expep =~ s/x$//i) {
	    print STDERR "Removed terminal 'X' from exon\n";
	}

        $self->throw("Exon translation not found in peptide") 
                    unless ($pep =~ /$expep/);

        $ex_align {'name'}      = $self->genomic->id;

	if ($exon->strand == 1) {
	  $ex_align {'gen_start'} = $exon->start + (3 - $exon->phase)%3;
	  $ex_align {'gen_end'}   = $exon->end   - $exon->end_phase;  
	  if( $stop_codon_present ) { $ex_align{'gen_end'} -= 3 } 
	} else {
	  $ex_align {'gen_start'} = $exon->start + $exon->end_phase;
	  $ex_align {'gen_end'}   = $exon->end   - (3 - $exon->phase) %3;  
	  if( $stop_codon_present ) { $ex_align{'gen_start'} += 3 } 
	}

        $ex_align {'strand'}    = $exon->strand;
        $ex_align {'phase'}     = $exon->phase;
        $ex_align {'end_phase'} = $exon->end_phase;
        $ex_align {'pep_start'} = index($pep, $expep)+1;
        $ex_align {'pep_end'}   = ($ex_align {'pep_start'} + length($expep))-1;
      
        push (@exon_aligns, \%ex_align);
        
        $dna_align {'exon_dna_limit'} += $exon->length;   
        
    }
    
    $dna_align {'pep_limit'} = $dna_align {'exon_dna_limit'}/3;      
    
    #map each feature to 1 or more exons
    foreach my $gapped_feature ( @features ) {
      my @split_features;
      print STDERR "DEBUG: processing $gapped_feature Cigar:",$gapped_feature->cigar_string(),"\n";

      foreach my $fp ( $gapped_feature->ungapped_features() ) {   
	unless (($fp->end - $fp->start)+1 <= $dna_align{'pep_limit'})
	  {
	    #$self->throw("Feature length (".$fp->start."-".$fp->end. 
	    #   ") is larger than peptide (".$dna_align{'pep_limit'}.")\n");
	  }
	#find each matching exon
	my (@aligned_exons);
	foreach my $ex_align (@exon_aligns) {
	  print STDERR "\tExon " . $ex_align->{gen_start} . " " . $ex_align->{gen_end} . " " . $ex_align->{pep_start} . " " . $ex_align->{pep_end} . "\n";
	  if (!($fp->end < $ex_align->{pep_start} || $fp->start > $ex_align->{pep_end})) {
	    push (@aligned_exons, $ex_align);
	  }
	}
	#create sets of featurepairs mapping peptide features to exons
	push(@split_features, $self->create_peptide_featurepairs($fp, @aligned_exons));
	
      }
      # hash on exon align to reconstruct exon based alignments
      # for each exon align we expect 1 or more seqfeatures, hence the
      # array
      my %exon_hash;
      foreach my $split ( @split_features ) {
	if( !defined $exon_hash{$split->{'_exon_align'}} ) {
	  $exon_hash{$split->{'_exon_align'}} = [];
	}
	push(@{$exon_hash{$split->{'_exon_align'}}},$split);
      }

      foreach my $exon_id ( keys %exon_hash ) {
	foreach my $sf ( @{$exon_hash{$exon_id}} ) {
	  print STDERR "DEBUG features pair ",$sf->start," ",$sf->end," ",$sf->hstart," ",$sf->hend,"\n";
	}

	my $dna_align_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new
	  (-features => $exon_hash{$exon_id});
	$self->featurepairs($dna_align_feature);
      }
    }
}


# This function takes a blast peptide feature hit and a set of matching exons and
# creates a set of featurepairs aligned to genomic coordinates. It will split
# features if they cross exon boundaries
sub create_peptide_featurepairs {
    my ($self, $fp, @aligned_exons) = @_;
    #create featurepairs
    
    my @output_features;

    print STDERR "\nConverting featurepair : PEP " . $fp->start . "\t" . $fp->end . " HIT " . $fp->hstart . "\t" . $fp->hend . "\t" . $fp->hstrand . "\n";

    foreach my $ex_align (@aligned_exons) {
      print "\nFound aligned exon " . $ex_align->{pep_start} . "\t" . $ex_align->{pep_end} . "\t" . $ex_align->{gen_start} . "\t" . $ex_align->{gen_end} . "\n";
      my ($ex_start, $ex_end, $dna_start, $dna_end, $start_phase, $end_phase);
      #This splits features across multiple exons and records phases
      if ($ex_align->{'pep_start'}  < $fp->start) {

#	print "\tFound start of feature\n";
	#feature starts inside current exon
	
	if ($ex_align->{strand} == 1) {
	  $ex_start   = $ex_align->{'gen_start'}
	  + (($fp->start - $ex_align->{'pep_start'})*3);
	} else {
	  $ex_end     = $ex_align->{'gen_end'}
	  - (($fp->start - $ex_align->{'pep_start'})*3);
	}

	$start_phase= 0;
	
	if ($fp->hstrand == 1) {
	  $dna_start  = $fp->hstart;
	} else {
	  $dna_end    = $fp->hend;
	}
	
      } else {
#	print "\tStart of feature = exon start\n";
	#feature starts in a previous exon or absolute start of current exon
	if ($ex_align->{strand} == 1) {
	  $ex_start   = $ex_align->{'gen_start'};
	} else {
	  $ex_end     = $ex_align->{'gen_end'};
	}
	$start_phase= $ex_align->{'phase'};
	
	# check for strand
	if ($fp->hstrand == 1) {
	  $dna_start  = $fp->hstart + ($ex_align->{'pep_start'} - $fp->start)*3;
	} else {
	  $dna_end    = $fp->hend   - ($ex_align->{pep_start}   - $fp->start)*3;
	}
      }
        
      if ($$ex_align{'pep_end'}    > $fp->end) {
#	print "\tFeature ends in exon\n";
	#feature ends in current exon
	if ($ex_align->{strand} == 1) {
	  $ex_end     = $ex_align->{'gen_start'}
	  + (($fp->end -  $ex_align->{'pep_start'})*3)+2;
	} else {
	  $ex_start   = $ex_align->{'gen_end'}
	  - (($fp->end -  $ex_align->{'pep_start'})*3)-2;
	}
	$end_phase  = 0;
	
	#Check for strand
	if ($fp->hstrand == 1) {
	  $dna_end    = $fp->hend;
	} else {
	  $dna_start  = $fp->hstart;
	}
      } else {
#	print "\tFeature end is exon end\n";
	#feature ends in a later exon or absolute end of current exon
	if ($ex_align->{strand} == 1) {
	  $ex_end     = $ex_align->{'gen_end'};
	} else {
	  $ex_start   = $ex_align->{'gen_start'};
	}
	$end_phase  = $ex_align->{'end_phase'};
	
	# Check for strand
	if ($fp->hstrand == 1) {
#	  print "Setting dna end\n";
	  $dna_end    = $fp->hstart + ($ex_align->{'pep_end'} - $fp->start)*3 +2;
	} else {
	  $dna_start  = $fp->hend   - ($ex_align->{pep_end}   - $fp->start)*3 -2;
        }
      }

      # Need to sort out strand - do we need strand or hstrand here?
      my $strand = 1;
      if ($ex_align->{strand} == 1 ) {
	$strand  = $fp->hstrand;
      } else {
	$strand  = $fp->hstrand * -1;
      }

      my $start_frac = $ex_align->{'phase'} + 1;
      my $end_frac   = (( 3 - $$ex_align{'end_phase'})%3) + 1;
      my $dna_feat = Bio::EnsEMBL::SeqFeature->new (
						    -seqname    =>  $ex_align->{'name'},
						    -start      =>  $ex_start, 
						    -end        =>  $ex_end,
						    -strand     =>  $strand,
						    -score      =>  $fp->score,
						    -p_value    =>  $fp->p_value,
						    -percent_id =>  $fp->percent_id,
						    -analysis   =>  $fp->analysis,
						    -primary_tag=>  $fp->primary_tag,
						    -source_tag =>  $fp->source_tag, 
						    -phase      =>  $start_phase,  
						    -end_phase  =>  $end_phase );
        
      my $pep_feat = Bio::EnsEMBL::SeqFeature->new (
							-seqname    =>  $fp->hseqname,
							-start      =>  $dna_start,
							-end        =>  $dna_end,
							-strand     =>  $strand,
							-score      =>  $fp->score,
							-p_value    =>  $fp->p_value,
							-percent_id =>  $fp->percent_id,
							-analysis   =>  $fp->analysis,
							-primary_tag=>  $fp->primary_tag,
							-source_tag =>  $fp->source_tag );
                                    
      my $featurepair = Bio::EnsEMBL::FeaturePair->new (
							-feature1   => $dna_feat,
							-feature2   => $pep_feat );

      $featurepair->{'_exon_align'} = $ex_align;

      $featurepair->attach_seq($self->genomic);
      print STDERR "After conversion ",$featurepair->start," ",$featurepair->end," ",$featurepair->hstart," ",$featurepair->hend,"\n";

      push(@output_features,$featurepair);


      
      print "\n" . $featurepair->gffstring .  " " . ($featurepair->feature1->end-$featurepair->feature1->start) . " " .( $featurepair->feature2->end-$featurepair->feature2->start) ."\n";

    }   

    return @output_features;
}


# Get/set functions for the data follow here

=head2 genomic

    Title   :   genomic
    Usage   :   $self->genomic
    Function:   Get/set for the genomic dna sequence
    Returns :   Bio::PrimarySeqI
    Args    :   Bio::PrimarySeqI

=cut

sub genomic {
    my($self,$seq) = @_;
    
    if (defined($seq)) {
      if (!($seq->isa("Bio::PrimarySeqI"))) {
	$self->throw("[$seq] is not a Bio::PrimarySeqI");
      }

      $self->{'_genomic'} = $seq;
    }

    return $self->{'_genomic'};
}

=head2 peptide

    Title   :   peptide
    Usage   :   $self->peptide($pep)
    Function:   Get/set for the peptide transcript
    Returns :   Bio::EnsEMBL::Transcript
    Args    :   Bio::EnsEMBL::Transcript

=cut

sub peptide {
    my($self,$seq) = @_;
    
    if (defined($seq)) {
      if (!($seq->isa("Bio::EnsEMBL::PredictionTranscript"))) {
	$self->throw("[$seq] is not a Bio::EnsEMBL::PredictionTranscript");
      }
      $self->{'_peptide'} = $seq;
    }

    return $self->{'_peptide'};
}

=head2 program

    Title   :   program
    Usage   :   $self->program('blastp');
    Function:   Get/set for the flavour of blast to run
    Returns :   String
    Args    :   String

=cut

sub program {
    my($self,$arg) = @_;
    
    if (defined($arg)) {
      $self->{'_program'} = $arg;
    }

    return $self->{'_program'};
}

=head2 database

    Title   :   database
    Usage   :   $self->database('swir');
    Function:   Get/set for the database to search against
    Returns :   String
    Args    :   String

=cut

sub database {
    my($self,$arg) = @_;
    
    if (defined($arg)) {
      $self->{'_database'} = $arg;
    }

    return $self->{'_database'};
}

=head2 threshold

    Title   :   threshold
    Usage   :   $self->threshold(100)
    Function:   Get/set for the score threshold to filter with
    Returns :   int
    Args    :   int

=cut


=head2 options

    Title   :   options
    Usage   :   $self->options('V=1000000');
    Function:   Get/set for the options to pass to the blast command
    Returns :   string
    Args    :   string

=cut


sub featurepairs {
    my ($self, $fp) = @_;
    if ($fp)
    {
        $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair") 
                unless $fp->isa("Bio::EnsEMBL::FeaturePairI");
        push (@{$self->{'_featurepairs'}}, $fp);
    }
    print STDERR   "FEATURES: ".(@{$self->{'_featurepairs'}})."\n";
    return @{$self->{'_featurepairs'}};
}


