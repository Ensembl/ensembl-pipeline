# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# based on Michelle Clamp's Blast.pm
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::BlastWorm

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $clonefile,
			   -format => 'Fasta');

  my $blast =  Bio::EnsEMBL::Pipeline::Runnable::BlastWorm->new 
    ('-clone'          => $query,
     '-program'        => 'wublastp' or '/usr/local/pubseq/bin/wublastp',
     '-database'       => 'swissprot',
     '-threshold'      => 1e-3,
     '-threshold_type' => 'PVALUE'
     '-filter'         => 1,
     '-options'        => 'V=1000000');

  $blast->run();
  @featurepairs = $blast->output();

  foreach my $fp (@featurepairs) {
      print $fp->gffstring."\n";
  }

  # If you have blast runs lying around that need parsing, you can do
  open(BLAST,"<blast.out");
  my @featurepairs = Bio::EnsEMBL::Pipeline::Runnable::BlastWorm->parse_results(\*BLAST);
  close(BLAST);


=head1 DESCRIPTION

  Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs blast; the
  output is parsed by BPLite and stored as Bio::EnsEMBL::FeaturePairs. 

=head1 CONTACT

  Describe contact details here

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BlastWorm;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter;
use Bio::Tools::BPlite;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

 Title    : new
 Usage    : my $seg =  Bio::EnsEMBL::Pipeline::Runnable::BlastpWorm->new ()
 Function : initialises BlastWorm object
 Returns  : a BlastWorm object
 Args     : 
 Throws   :

=cut

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new (@_);    
  
    $self->{'_flist'} = [];               # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef;        # location of Bio::Seq object
    $self->{'_database'}  = undef;        # name of database
    $self->{'_program'} = undef;          # location of program executable
    $self->{'_threshold'} = undef;        # threshold for hit filterting
    $self->{'_options'}   = undef;        # additional arguments for blast
    $self->{'_filter'}    = 0;            # do we filter features?
    $self->{'_workdir'}   = undef;        # location of tmp directory
    $self->{'_filename'}  = undef;        # file to store Bio::Seq object
    $self->{'_results'}   = undef;        # file to store results of seg run
    $self->{'_protected'} = [];           # a list of files protected from deletion
  
    my ($clone, $program, $database,
        $threshold,  $threshold_type, $options, $filter) = $self->_rearrange([qw(CLONE 
                                                                                 PROGRAM
                                                                                 DATABASE
                                                                                 THRESHOLD
                                                                                 THRESHOLD_TYPE
                                                                                 OPTIONS         
                                                                                 FILTER)], 
                                                                              @args);
  
    $self->clone ($clone) if ($clone);       
    $self->program ($self->find_executable ($program));
  
    if ($database) {
        $self->database($database);
    } else {
        $self->throw("BlastWorm needs a database");
    }
    
    $self->threshold ($threshold) if ($threshold);
    $self->threshold_type ($threshold_type) if ($threshold_type);
    $self->options ($options) if ($options);
    $self->filter ($filter) if ($filter);

    return $self;
}


###################
# get/set methods 
###################

=head2 clone

 Title    : clone
 Usage    : $self->clone ($clone);
 Function : get/set method for the Sequence object; assigns clone,
            seq-filename and result-filename
 Example  :
 Returns  : a Bio::Seq or Bio::PrimarySeq object
 Args     : a Bio::Seq or Bio::PrimarySeq object (optional)
 Throws   :

=cut

sub clone {
    my ($self, $seq) = @_;
    if ($seq) {
        ($seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI"))
            || $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
        $self->{'_sequence'} = $seq ;
        $self->clonename ($self->clone->id);
        $self->filename ($self->clone->id.".$$.seq");
        $self->results ($self->filename.".out");
    }
    return $self->{'_sequence'};
}


=head2 program

 Title    : program
 Usage    : $self->program ('/usr/local/pubseq/bin/wublastp');
 Function : get/set method for the path to the program binaries
 Example  :
 Returns  : File path
 Args     : File path (optional)
 Throws   :

=cut

sub program {
    my ($self, $location) = @_;
    if ($location) {
        unless (-e $location) {
            $self->throw ($self->program." not found at $location");
        }
        $self->{'_program'} = $location ;
    }
    return $self->{'_program'};
}


=head2 database

 Title    : database
 Usage    : $self->database ($database);
 Function : get/set method for the analysisId
 Example  :
 Returns  : analysisId
 Args     : analysisId (optional)
 Throws   :

=cut

sub database {
    my $self = shift;
    if (@_) {
        $self->{'_database'} = shift;
    }
    return $self->{'_database'};
} 


=head2 threshold

 Title    : threshold
 Usage    : $self->threshold ($threshold);
 Function : get/set method for the analysisId
 Example  :
 Returns  : analysisId
 Args     : analysisId (optional)
 Throws   :

=cut

sub threshold {
    my $self = shift;
    if (@_) {
        $self->{'_threshold'} = shift;
    }
    return $self->{'_threshold'};
} 


=head2 threshold_type

 Title    : threshold_type
 Usage    : $self->threshold_type ($threshold_type);
 Function : get/set method for the analysisId
 Example  :
 Returns  : analysisId
 Args     : analysisId (optional)
 Throws   :

=cut

sub threshold_type {
    my $self = shift;
    if (@_) {
        $self->{'_threshold_type'} = shift;
    }
    return $self->{'_threshold_type'};
}


=head2 options

 Title    : options
 Usage    : $self->options ($options);
 Function : get/set method for the analysisId
 Example  :
 Returns  : analysisId
 Args     : analysisId (optional)
 Throws   :

=cut

sub options {
    my $self = shift;
    if (@_) {
        $self->{'_options'} = shift;
    }
    return $self->{'_options'};
} 


=head2 filter

 Title    : filter
 Usage    : $self->filter ($filter);
 Function : get/set method for the analysisId
 Example  :
 Returns  : analysisId
 Args     : analysisId (optional)
 Throws   :

=cut

sub filter {
    my $self = shift;
    if (@_) {
        $self->{'_filter'} = shift;
    }
    return $self->{'_filter'};
} 


###################
# analysis methods
###################

=head2 run

 Title    : run
 Usage    : $self->run ($workdir, $args)
 Function : runs program and populates @{$self->{'_flist'}} (array of features)
 Example  :
 Returns  :   
 Args     : workdir, args (optional)
 Throws   :

=cut

sub run {
    my ($self, $dir, $args) = @_;

    # nothing to be done with $args

    # check clone
    my $seq = $self->clone || $self->throw("Clone required for ".$self->program."\n");

    # set directory if provided
    $self->workdir ('/tmp') unless ($self->workdir($dir));
    $self->checkdir;

    # reset filename and results as necessary (adding the directory path)
    my $tmp = $self->workdir;
    my $input = $tmp."/".$self->filename;
    $self->filename ($input);
    $tmp .= "/".$self->results;
    $self->results ($tmp);

    # write sequence to file
    $self->writefile;        

    # run program
    $self->run_program;

    # parse output
    $self->parse_results;
    $self->deletefiles;
}


=head2 run_program

 Title    : run_program
 Usage    : $self->program
 Function : makes the system call to program
 Example  :
 Returns  : 
 Args     :
 Throws   :

=cut

sub run_program {
    my ($self) = @_;

    # this routine expands the database name into $db-1 etc.
    # for split databases
    my @databases = $self->fetch_databases;

    # run program
    foreach my $database (@databases) {
        print STDERR "running ".$self->program." against ".$self->database."\n";

        my $cmd = $self->program .' '. 
	   	  $database      .' '.
		  $self->filename.' '.
		  $self->options . ' >> '. 
		  $self->results;
        print STDERR "$cmd\n";   

        $self->throw ("Error running ".$self->program." on ".$self->filename." against ".$self->database) 
            unless ((system ($cmd)) == 0);
    }
}


=head2 fetch_databases

 Title    : fetch_databases
 Usage    : $self->fetch_databases
 Function : checks whether the database in $self->database
            exists.  If it doesn\'t it checks whether the
            database has been split into $database-1, $database-2 etc.
 Example  : 
 Returns  : array of database names
 Args     :
 Throws   :

=cut

sub fetch_databases {
    my ($self) = @_;
    
    my @databases;
    my $fulldbname;

    # If we have passed a full path name don't append the $BLASTDB
    # environment variable.
    if ($self->database =~ /\//) {
	$fulldbname = $self->database;
    } else {
	$fulldbname = $ENV{BLASTDB}."/".$self->database;
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
	    push(@databases,$fulldbname."-$count");
	    $count++;
	}
    }
    if (scalar(@databases) == 0) {
	$self->throw("No databases exist for ".$self->database);
    }

    return @databases;
}


=head2 parse_results

 Title    : parse_results
 Usage    : $self->parse_results ($filename)
 Function : Parses the Blast results and stores the output in
            an array of feature pairs.
            Giving a filehandle as argument makes this method 
            read from this file rather than the filename 
            stored in the object.
 Example  :
 Returns  : 
 Args     : filehandle (optional)
 Throws   :

=cut

sub parse_results {
    my ($self, $fh) = @_;

    my $filehandle;

    if ($fh) {
        $filehandle = $fh;
    }
    else {
        my $resfile = $self->results;
    
        if (-e $resfile) {
            # it's a filename
            if (-z $self->results) {  
                print STDERR $self->program." didn't find anything\n";
                return;
            }       
            else {
                open (BLASTOUT, "<$resfile") or $self->throw ("Error opening $resfile");
                $filehandle = \*BLASTOUT;
            }
        }
        else {
            # it'a a filehandle
            $filehandle = $resfile;
        }
    }

    # BPlite does most of the work
    my $parser = Bio::Tools::BPlite->new (-fh => $filehandle);

    # Loop over each blast hit
    my %ids;

    if ($self->filter) {
        print STDERR "filtering hits: ";
        %ids = $self->filter_hits($parser);

        close $filehandle;

        if ($fh) {
            $filehandle = $fh;
        }
        else {
            my $resfile = $self->results;
    
            if (-e $resfile) {
                # it's a filename
                if (-z $self->results) {  
                    print STDERR $self->program." didn't find anything\n";
                    return;
                }       
                else {
                    open (BLASTOUT, "<$resfile") or $self->throw ("Error opening $resfile");
                    $filehandle = \*BLASTOUT;
                }
            }
            else {
                # it'a a filehandle
                $filehandle = $resfile;
            }
        }

        # BPlite does most of the work
        $parser = Bio::Tools::BPlite->new (-fh => $filehandle);
    }

    SBJCT:while (my $sbjct = $parser->nextSbjct)  {
  
        # check if the match has been filtered out
        if (($self->filter == 1) && !exists($ids{$sbjct->name})) {
	    next SBJCT;
        }
	  
        # loop over all hsp's
        # (hsp's are Featurepairs (pairs of Bio::SeqFeatureI's), and
        # $hsp->query, $hsp->subject is equal to $hsp->feature1, $hsp->feature2)
        HSP: while (my $hsp = $sbjct->nextHSP) {
            if ($self->threshold_type eq "PID") {
	        next HSP if ($hsp->percent < $self->threshold);
	    } elsif ($self->threshold_type eq "PVALUE") {
	        next HSP if ($hsp->P > $self->threshold);
	    }

            my ($strand,$hstrand) = $self->_findStrands ($hsp);

            # Each HSP is a gapped alignment:
	    # split the gapped alignment into ungapped subalignments
  	    my @align_coordinates = $self->split_HSP($hsp);

            my %feature;
            ($feature{name}) = $parser->query =~ /^(\S+)/;
            $feature{score} = $hsp->score;
            $feature{p_value} = sprintf ("%.3e", $hsp->P);
            $feature{percent_id} = $hsp->percent;
            $feature{start} = $hsp->query->start;
            $feature{end} = $hsp->query->end;
            $feature{strand} = $strand;
            ($feature{hname}) = $sbjct->name =~ /^(\S+)/;
            $feature{hstart} = $hsp->subject->start;
            $feature{hend} = $hsp->subject->end;
            $feature{hstrand} = $hstrand;
            ($feature{source}) = $self->program =~ /([^\/]+)$/;
            $feature{primary} = 'similarity';
            ($feature{program}) = $self->program =~ /([^\/]+)$/;
            ($feature{db}) = $self->database =~ /([^\/]+)$/;
            ($feature{logic_name}) = $self->program =~ /([^\/]+)$/;
            $feature {align_coor} = \@align_coordinates;
            $self->create_feature (\%feature);
        }
    }
    close $filehandle;
}


=head2 filter_hits

 Title    : filter_hits
 Usage    : $self->filter_hits ($parser);
 Function : 
 Example  :
 Returns  : hash of id\'s that passed the filtering step
 Args     : a Bio::Tools::BPlite object
 Throws   :

=cut

sub filter_hits {
    my ($self,$parser) = @_;

    my %ids;

    my @features;
    NAME: while(my $sbjct = $parser->nextSbjct)  {
      
        my $name = $sbjct->name ;
        
        HSP: while (my $hsp = $sbjct->nextHSP) {
            if ($self->threshold_type eq "PID") {
	        next HSP if ($hsp->percent < $self->threshold);
            } elsif ($self->threshold_type eq "PVALUE") {
	        next HSP if ($hsp->P > $self->threshold);
	    }

  	    my $qstart = $hsp->query->start();
	    my $hstart = $hsp->subject->start();
	
	    my $qend   = $hsp->query->end();
	    my $hend   = $hsp->subject->end();      

	    my ($qstrand,$hstrand) = $self->_findStrands   ($hsp);

	    my $score  = $hsp->score;

	    my $feature1 = new Bio::EnsEMBL::SeqFeature();
	    $feature1->start($qstart);
	    $feature1->end  ($qend);
	    $feature1->strand($qstrand);
	    $feature1->score($score);
	    $feature1->source_tag('tmp');
	    $feature1->primary_tag('similarity');

	    my $feature2 = new Bio::EnsEMBL::SeqFeature();
	    $feature2->start  ($hstart);
	    $feature2->end    ($hend);
	    $feature2->strand ($hstrand);
	    $feature2->score  ($score);
	    $feature2->seqname($name);
	    $feature2->source_tag('tmp');
	    $feature2->primary_tag('similarity');

	    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
		  			           -feature2 => $feature2);

	    push(@features,$fp);
        }
    }

    # don't give the e/p-value to the FeatureFilter for the moment and
    # set the score threshold very low => filter only for coverage
    my $search = Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter->new ( -coverage  => 10,
                                                                        -minscore  => -1000,
                                                                        -maxevalue => 0.1,
                                                                      );
    my @newfeatures = $search->run(@features);

    print STDERR "there were ".@features." HSP's before, and ".@newfeatures." after\n\n";

    foreach my $f (@newfeatures) {
        my $id = $f->hseqname;
        $ids{$id} = 1;
    }

  return %ids;
}
    

=head2 split_HSP

 Title    : split_HSP
 Usage    : $self->split_HSP
 Function : takes a gapped blast HSP,
            and turns it into an array of ungapped subalignments
 Example  : 
 Returns  : 
 Args     : a BPlite::HSP
 Throws   :

=cut

sub split_HSP {
    my ($self, $hsp) = @_;

    my @align_coordinates;

    # First of all some jiggery pokery to find out what sort of alignment
    # we have - (dna-dna, dna-pep, pep-dna etc).
    # We also work out which strand each sequence is on.
    # 
    # Both these variables are needed to work out the increment
    # values for the query and the hit sequences.  As we split
    # up the gapped alignment we need to increment the coordinates.
    #
    # We have to take care with dna-pep alignments to increment the dna
    # sequences by 3 as there are 3 bases in a codon.
    #
    # For dna-dna alignments this will be 1 or -1 depending on
    # the strand.  
    #
    # For pep-dna alignments and vice-versa the increments will be +-3 
    # for the dna sequence and +- 1 for the peptide sequence. 

    my ($qtype,  $htype)   = $self->_findTypes     ($hsp);
    my ($qstrand,$hstrand) = $self->_findStrands   ($hsp);
    my ($qinc,   $hinc)    = $self->_findIncrements($hsp,$qstrand,$hstrand,$qtype,$htype);

#    print STDERR "\t$qtype ($qinc); $htype ($hinc); qstrand $qstrand, hstrand $hstrand\n";

    # We split the alignment strings into arrays of one char each.  
    # We then loop over this array and when we come to a gap
    # in either the query sequence or the hit sequence we make 
    # a new feature pair. 
    #

    # Before making a new feature pair we check whether we have
    # something to make a feature pair out of - i.e. we aren't just on
    # the 2nd or greater base in a gapped region.  
    #
    # As we loop over the array we need to keep track of both the
    # query coordinates and the hit coordinates.  We track the last
    # start of a feature pair and also the current end of a feature
    # pair.  The feature pair start is reset every time we hit a
    # gapped position and the feature pair end is reset every time we
    # hit an aligned region.

    my @gap;

    my @qchars = split(//,$hsp->querySeq);        # split query alignment into array of chars
    my @hchars = split(//,$hsp->sbjctSeq);        # ditto for hit
    
    my $qstart = ($qstrand == 1 ) ? $hsp->query->start : $hsp->query->end;             # set the query start
    my $hstart = ($hstrand == 1 ) ? $hsp->subject->start : $hsp->subject->end;         # ditto for hit
    
    my $qend   = $qstart;                         # set the ungapped subalignment end
    my $hend   = $hstart;                         # ditto for hit
    
    my $count = 0;                                # counter for the bases in the alignment
    my $found = 0;                                # flag saying whether we have an ungapped subalignment
    
    while ($count <= $#qchars) {
	# We have hit an ungapped region. Increase the query and hit counters.
	# and set the subalignment flag
	if ($qchars[$count] ne '-' && $hchars[$count] ne '-') {
	    $qend += $qinc;
	    $hend += $hinc;
	    $found = 1;
	} else {
	    # We have hit a gapped region. If the subalignment flag is set,
	    # push the coordinates on the array, and reset the start and end variables
	    if ($found == 1) {
                push (@align_coordinates, [$qstart, $hstart]);
	    }
	
	    # We're in a gapped region. We need to increment the sequence that
	    # doesn't have the gap in it to keep the coordinates correct.
	    # We also need to reset the current end coordinates.

	    if ($qchars[$count] ne '-') {
		$qstart = $qend+$qinc;
	    } else {
		$qstart = $qend;
	    }
	    if ($hchars[$count] ne '-') {
		$hstart = $hend+$hinc;
	    } else {
		$hstart = $hend;
	    }
	    $qend = $qstart;
	    $hend = $hstart;
	    $found = 0;
	}
	$count++;
    }			     
    # Remember the last feature
    if ($found == 1) {
        push (@align_coordinates, [$qstart, $hstart]);
    }

    # take the HDP start coordinates off the array
    shift @align_coordinates;
    return @align_coordinates;
}


sub _findIncrements {
    my ($self,$hsp,$qstrand,$hstrand,$qtype,$htype) = @_;

    my $qinc   = 1 * $qstrand;
    my $hinc   = 1 * $hstrand;
    
    if ($qtype eq 'dna' && $htype eq 'pep') {
	$qinc = 3 * $qinc;
    } 
    if ($qtype eq 'pep' && $htype eq 'dna') {
	$hinc = 3 * $hinc;
    }

    return ($qinc,$hinc);
}

sub _findStrands {
    my ($self,$hsp) = @_;

    return ( $hsp->query->strand(),
	     $hsp->subject->strand());
}

sub _findTypes {
    my ($self,$hsp) = @_;

    my $type1;
    my $type2;
    my $len1 = $hsp->query->length();
    my $len2 = $hsp->subject->length();

    if ($len1/$len2 > 2) {
	$type1 = 'dna';
	$type2 = 'pep';
    } elsif ($len2/$len1 > 2) {
	$type1 = 'pep';
	$type2 = 'dna';
    } else {
	$type1 = 'dna';
	$type2 = 'dna';
    }

    return ($type1,$type2);
}


=head2 create_feature

 Title    : create_feature
 Usage    : $self->create_feature ($feature)
 Function : creates a Bio::EnsEMBL::FeaturePair object from %feature,
            and pushes it onto @{$self->{'_flist'}}
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub create_feature {
    my ($self, $feat) = @_;

    # create analysis object (will end up in the analysis table)
    my $analysis = Bio::EnsEMBL::Analysis->new ();
    $analysis->db ($feat->{db});
    $analysis->program ($feat->{program});
    $analysis->gff_source ($feat->{source});
    $analysis->gff_feature ($feat->{primary});
    $analysis->logic_name ($feat->{logic_name});

    # create featurepair object
    my $feature1 = Bio::EnsEMBL::SeqFeature->new ();
    $feature1->seqname ($feat->{name});
    $feature1->start ($feat->{start});
    $feature1->end ($feat->{end});
    $feature1->strand ($feat->{strand});
    $feature1->score ($feat->{score});
    $feature1->p_value ($feat->{p_value});
    $feature1->percent_id ($feat->{percent_id});
    $feature1->source_tag ($feat->{source});
    $feature1->primary_tag ($feat->{primary});
    $feature1->analysis ($analysis);

    $feature1->add_tag_value ('align_coor', $feat->{align_coor});

    my $feature2 = Bio::EnsEMBL::SeqFeature->new ();
    $feature2->seqname ($feat->{hname});
    $feature2->start ($feat->{hstart});
    $feature2->end ($feat->{hend});
    $feature2->strand ($feat->{hstrand});
    $feature2->score ($feat->{score});
    $feature2->p_value ($feat->{p_value});
    $feature2->percent_id ($feat->{percent_id});
    $feature2->source_tag ($feat->{source});
    $feature2->primary_tag ($feat->{primary});
    $feature2->analysis ($analysis);

    my $featurepair = Bio::EnsEMBL::FeaturePair->new ();
    $featurepair->feature1 ($feature1);
    $featurepair->feature2 ($feature2);

    if ($featurepair) {
        $featurepair->feature1->validate_prot_feature (1);
        $featurepair->feature2->validate_prot_feature (1);
        # add to _flist
        push (@{$self->{'_flist'}}, $featurepair);
    }
}


=head2 output

 Title    : output
 Usage    : $self->output
 Function : returns an array of featurepair objects
 Example  :
 Returns  : an array of Bio::EnsEMBL::FeaturePair objects
 Args     :
 Throws   :

=cut

sub output {
    my ($self) = @_;
    my @list = @{$self->{'_flist'}};
    return @{$self->{'_flist'}};
}

1;
