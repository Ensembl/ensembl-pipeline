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

  my $blastgenscan = Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep->new ( 
                                                    -genomic    => $genomic,
                                                    -peptide    => $pep,
                                                    -database   => 'swir',
                                                    -threshold  => 1e-6,
                                                    -options    => 'B=1000'
                                                    );
  $blastgenscan->run();
  $blastgenscan->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides
constructed from assembling genscan predicted features to peptide
sequence. The resulting blast hits are written back as FeaturePairs.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::SearchFilter;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Pep_SeqFeature;
use Bio::EnsEMBL::Pipeline::Runnable::SearchFilter;
#use Data::Dumper;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-genomic     => $genomic
                           -peptide     => $peptide
                           -program     => 'blastp',
			   -database    => 'swir',
			   -threshold   => 1e-6);
			   
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep object
    Args    :   -genomic : Bio::Seq containing genomic dna
                -pep     : Bio::EnsEMBL::Transcript
                -parameters : ref. to a hash of parameters

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
    $self->_rearrange ([qw(GENOMIC PEPTIDE PROGRAM DATABASE THRESHOLD THRESHOLD_TYPE
			   OPTIONS)], @args);
  
  if (defined($genomic) && $genomic->isa("Bio::PrimarySeqI")) {
    $self->genomic($genomic);
  } else {
    $self->throw("No genomic sequence input");
  }
  
  if (defined($peptide) && $peptide->isa("Bio::EnsEMBL::Transcript")) {
    $self->peptide($peptide);
  } else {
    $self->throw("No peptide input");
  }
  
  if (defined($program)) {
    $self->program($program);
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
    Usage   :   $self->run;
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep->run
    Returns :   none
    Args    :   none

=cut

sub run {
  my ($self) = @_;
  
  my $transcript = $self->peptide;

  if (!defined($transcript)) {
    $self->throw("No peptide input");
  }
  
  print STDERR "Creating BioPrimarySeq for peptide ".$transcript->id."\n";

  my $peptide = Bio::PrimarySeq->new(-id         => $transcript->id,
				     -seq        => $transcript->translate->seq(),
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
							       -options   => $self->options,
                                                               -filter    => 1);
  $runnable->run();
  
  $self->align_hits_to_contig($runnable->output);
#  $self->check_features($transcript->translate->seq,$self->featurepairs);
  
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    return $self->featurepairs();  
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
    $trans->sort;

#    print STDERR "Peptide translation is " . $trans->translate->seq . "\n";


    # Calculate boundaries and map exons to translated peptide
    # Each exon is trimmed at either end so it has a whole number
    # of residues

    my $pep = $trans->translate->seq;

    foreach my $exon ($trans->each_Exon) {

        my %ex_align;

#	print STDERR "\tGenscan exon " . $exon->start . "\t" . $exon->end . "\t" . $exon->strand . "\t" . $exon->phase . "\n";
	my $strand = "+";
	if ($exon->strand == -1) {
	  $strand = "-";
	}
#	print STDERR "exon\tgenscan\tsimilarity\t" . $exon->start . "\t" . $exon->end . "\t100\t" . $strand ."\t" . $exon->phase . "\t" . $trans->id . "\n";

        my ($expep) = $exon->translate->seq =~ /[^\*]+/g;
	if ($expep =~ s/x$//i) {
	    print STDERR "Removed terminal 'X' from exon @{[$exon->id]}\n";
	}

        $self->throw("Exon translation not found in peptide") 
                    unless ($pep =~ /$expep/);

        $ex_align {'name'}      = $exon->id;

	# Trim the start and end of the exon
	if ($exon->strand == 1) {
	  $ex_align {'gen_start'} = $exon->start + (3 - $exon->phase)%3;
	  $ex_align {'gen_end'}   = $exon->end   - $exon->end_phase;  
	} else {
	  $ex_align {'gen_start'} = $exon->start + $exon->end_phase;
	  $ex_align {'gen_end'}   = $exon->end   - (3 - $exon->phase) %3;  
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
    foreach my $fp (@features) {
      unless (($fp->end - $fp->start)+1 <= $dna_align{'pep_limit'}) {
#	$self->throw("Feature length (".$fp->start."-".$fp->end. 
#		     ") is larger than peptide (".$dna_align{'pep_limit'}.")\n");
      }

      #find each matching exon
      my @aligned_exons;

      foreach my $ex_align (@exon_aligns) {
#	print STDERR "\texon " . $fp->start . "\t" . $fp->end . "\t" . $ex_align->{pep_start} . "\t" . $ex_align->{pep_end} . "\n";
	if (!($fp->end < $ex_align->{pep_start} || $fp->start > $ex_align->{pep_end})) {
	  push (@aligned_exons, $ex_align);
	}
      }
      #create sets of featurepairs mapping peptide features to exons
      $self->create_peptide_featurepairs($fp, @aligned_exons);
      
    }
  }

# This function takes a blast peptide feature hit and a set of matching exons and
# creates a set of featurepairs aligned to genomic coordinates. It will split
# features if they cross exon boundaries
sub create_peptide_featurepairs {
    my ($self, $fp, @aligned_exons) = @_;
    #create featurepairs
    
    #print "Converting featurepair : PEP " . $fp->start . "\t" . $fp->end . " HIT " . $fp->hstart . "\t" . $fp->hend . "\n";

    foreach my $ex_align (@aligned_exons)
    {
#      print "Found aligned exon " . $ex_align->{pep_start} . "\t" . $ex_align->{pep_end} . "\t" . $ex_align->{gen_start} . "\t" . $ex_align->{gen_end} . "\n";

        my ($ex_start, $ex_end, $pep_start, $pep_end, $start_phase, $end_phase);

      #This splits features across multiple exons and records phases

      
      if ($ex_align->{'pep_start'}  < $fp->start) {
	#feature starts inside current exon

	if ($ex_align->{strand} == 1) {
	  $ex_start   = $ex_align->{'gen_start'}
	  + (($fp->start - $ex_align->{'pep_start'})*3);
	} else {
	  $ex_end     = $ex_align->{'gen_end'}
	  - (($fp->start - $ex_align->{'pep_start'})*3);
	}

	$start_phase = 0;
	$pep_start   = $fp->hstart;
      } else {
	#feature starts in a previous exon or absolute start of current exon

	if ($ex_align->{strand} == 1) {
	  $ex_start   = $ex_align->{'gen_start'};
	} else {
	  $ex_end     = $ex_align->{'gen_end'};
	}
	$start_phase = 0;
	$pep_start  = $fp->hstart + ($ex_align->{'pep_start'} - $fp->start);
      }
        
      if ($$ex_align{'pep_end'}    > $fp->end) {
	#feature ends in current exon
	if ($ex_align->{strand} == 1) {
	  $ex_end     = $ex_align->{'gen_start'}
	  + (($fp->end -  $ex_align->{'pep_start'})*3)+2;
	} else {
	  $ex_start   = $ex_align->{'gen_end'}
	  - (($fp->end -  $ex_align->{'pep_start'})*3)-2;
	}
	$end_phase  = 0;
	$pep_end    = $fp->hend;
      } else {

	#feature ends in a later exon or absolute end of current exon
	if ($ex_align->{strand} == 1) {
	  $ex_end     = $ex_align->{'gen_end'};
	} else {
	  $ex_start   = $ex_align->{'gen_start'};
	}
	$end_phase  = $ex_align->{'end_phase'};
	$pep_end    = $fp->hstart + ($ex_align->{'pep_end'} - $fp->start);
	
      }
        
      my $start_frac = $ex_align->{'phase'} + 1;
      my $end_frac   = (( 3 - $$ex_align{'end_phase'})%3) + 1;
      my $dna_feat = Bio::EnsEMBL::SeqFeature->new (-seqname    =>  $ex_align->{'name'},
						    -start      =>  $ex_start, 
						    -end        =>  $ex_end,
						    -strand     =>  $ex_align->{'strand'},
						    -score      =>  $fp->score,
						    -p_value    =>  $fp->p_value,
						    -percent_id =>  $fp->percent_id,
						    -analysis   =>  $fp->analysis,
						    -primary_tag=>  $fp->primary_tag,
						    -source_tag =>  $fp->source_tag, 
						    -phase      =>  $start_phase,  
						    -end_phase  =>  $end_phase );
      
	my $pep_feat = Bio::EnsEMBL::SeqFeature->new (-seqname    =>  $fp->hseqname,
							-start      =>  $pep_start,
							-end        =>  $pep_end,
							-strand     =>  $ex_align->{'strand'},
							-start_frac =>  $start_frac,
							-end_frac   =>  $end_frac,
							-score      =>  $fp->score,
							-p_value    =>  $fp->p_value,
							-percent_id =>  $fp->percent_id,
							-analysis   =>  $fp->analysis,
							-primary_tag=>  $fp->primary_tag,
							-source_tag =>  $fp->source_tag );
      
      my $featurepair = Bio::EnsEMBL::FeaturePair->new (-feature1   => $dna_feat,
							-feature2   => $pep_feat );

      $featurepair->attach_seq($self->genomic);
      
      $self->featurepairs($featurepair);    

#      print "\n" . $featurepair->gffstring . "\n";

    }
}

sub check_features {
  my ($self,$pep,@f) = @_;

  print STDERR "Peptide is " . $pep . "\n";

  my %seqhash;

  foreach my $f (@f) {
    if (!defined($seqhash{$f->hseqname})) {
      my $seq = $self->get_Sequence($f->hseqname);
      $seqhash{$f->hseqname} = $seq;
    }
    my $seq = $seqhash{$f->hseqname};
    my $fdna = $f->seq;
    my $fpep = $fdna->translate('*','X',$f->phase);
    my $hpep = substr($seq->seq,$f->hstart-1,($f->hend - $f->hstart + 1));
    my $pos  = index($pep,$fpep->seq);

    print "\tFeature " . $f->start . "\t" . $f->end . "\t" . $f->strand . "\t" . $f->phase . "\t" . $f->hstart . "\t" . $f->hend . " " . $pos . "\n" ;
      print $fpep->seq . "\n" . $hpep . "\n";

  }
}

sub get_Sequence {
    my ($self,$id) = @_;


    next ID unless defined($id);

    print(STDERR "Sequence id :  is [$id]\n");

    open(IN,"pfetch -q $id |") || $self->throw("Error fetching sequence for id [$id]");
	
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
      if (!($seq->isa("Bio::EnsEMBL::Transcript"))) {
	$self->throw("[$seq] is not a Bio::EnsEMBL::Transcript");
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

sub threshold {
    my($self,$arg) = @_;
    
    if (defined($arg)) {
      $self->{'_threshold'} = $arg;
    }

    return $self->{'_threshold'};
}

sub threshold_type {
    my($self,$arg) = @_;
    
    if (defined($arg)) {
      $self->{'_threshold_type'} = $arg;
    }

    return $self->{'_threshold_type'};
}

=head2 options

    Title   :   options
    Usage   :   $self->options('V=1000000');
    Function:   Get/set for the options to pass to the blast command
    Returns :   string
    Args    :   string

=cut

sub options {
    my($self,$arg) = @_;
    
    if (defined($arg)) {
      $self->{'_options'} = $arg;
    }

    return $self->{'_options'};
}


=head2 featurepairs

    Title   :   featurepairs
    Usage   :   @out = $obj->featurepairs
    Function:   Returns the output featurepairs after they\'ve
                been converted back to genomic coords
                Also adds feature pairs to the array
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   Bio::EnsEMBL::FeaturePair

=cut


sub featurepairs {
    my ($self, $fp) = @_;
    if ($fp)
    {
        $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair") 
                unless $fp->isa("Bio::EnsEMBL::FeaturePairI");
        push (@{$self->{'_featurepairs'}}, $fp);
    }

    return @{$self->{'_featurepairs'}};
}

1;
