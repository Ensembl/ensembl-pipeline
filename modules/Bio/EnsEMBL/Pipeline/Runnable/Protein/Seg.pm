# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg

=head1 SYNOPSIS

  my $seqstream = Bio::SeqIO->new ( -file => $queryfile,
                                    -fmt => 'Fasta',
                                  );
  $seq = $seqstream->next_seq;

  my $seg = Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg->new ( -QUERY => $seq);
  $seg->workdir ($workdir);
  $seg->run;
  my @results = $seg->output;

=head1 DESCRIPTION

  Seg takes a Bio::Seq (or Bio::PrimarySeq) object
  and runs seg on it (detecting low complexity sequences). 
  The resulting output file is parsed to produce a set of features.

=head1 CONTACT
  
  Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

  The rest of the documentation details each of the object methods. 
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

 Title    : new
 Usage    : my $seg =  Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg->new
                       ( -program    => '/usr/local/pubseq/bin/seg',
                         -query      => $query,
                         -analysis   => $analysis,
                       );
 Function : initialises Seg object
 Returns  : a Seg object
 Args     : a Bio::Seq object, the path to the program binaries and an analysisId
 Throws   :

=cut

sub new {
    my ($class, @args) = @_;
  
    my $self = $class->SUPER::new (@_);    
  
    $self->{'_flist'}     = [];           # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef;        # location of Bio::Seq object
    $self->{'_program'}   = undef;        # location of executable
    $self->{'_workdir'}   = undef;        # location of tmp directory
    $self->{'_filename'}  = undef;        # file to store Bio::Seq object
    $self->{'_results'}   = undef;        # file to store results of program run
    $self->{'_protected'} = [];           # a list of files protected from deletion
  
    my ($query, $analysis, $program) = $self->_rearrange([qw(QUERY 
						             ANALYSIS
                                                             PROGRAM)], 
					                 @args);
    $self->query ($query) if ($query);

    if ($analysis) {
        $self->analysis ($analysis);
    }
    else {
        $self->throw("Seg needs an analysis");
    }

    $self->program ($self->find_executable ($self->analysis->program_file));
  
    return $self;

}

###################
# get/set methods 
###################

=head2 query

 Title    : query
 Usage    : $self->query ($query);
 Function : get/set method for the Sequence object; assigns query, filename and results
 Example  :
 Returns  : a Bio::Seq or Bio::PrimarySeq object
 Args     : a Bio::Seq or Bio::PrimarySeq object (optional)
 Throws   :

=cut

sub query {
    my ($self, $seq) = @_;
    if ($seq) {
	eval {
	    $seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")
	};

	if (!$@) {
	    $self->{'_sequence'} = $seq ;
	    $self->filename ("query.$$.seq");
	    $self->results ($self->filename.".out");
	}
	else {
	    print STDERR "SEG WARNING: The input_id is not a Seq object but if its a peptide fasta file, it should go fine\n";
	    $self->{'_sequence'} = $seq ;
	    $self->filename ("$$.tmp.seq");
	    
	    $self->results ("seg.$$.out");
	    
	}
    }
    return $self->{'_sequence'};
}


=head2 analysis

 Title    : analysis
 Usage    : $self->analysis ($analysis);
 Function : get/set method for the analysis
 Example  :
 Returns  : analysis
 Args     : analysis (optional)
 Throws   :

=cut

sub analysis {
    my $self = shift;
    if (@_) {
        $self->{'_analysis'} = shift;
    }
    return $self->{'_analysis'};
} 


=head2 program

 Title    : program
 Usage    : $self->program ('/usr/local/pubseq/bin/hmmpfam');
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

####################
# analysis methods
####################

=head2 run

 Title    : run
 Usage    : $self->run ($workdir, $args)
 Function : runs program and populates @{$self->{'_flist'}} (array of features)
 Example  :
 Returns  :   
 Args     : workdir (optional)
 Throws   :

=cut

sub run {
    my ($self, $dir) = @_;

    # check query
    my $seq = $self->query || $self->throw("Query required for Program");

    # set directory if provided
    $self->workdir ('/tmp') unless ($self->workdir($dir));
    $self->checkdir;

    # reset filename and results as necessary (adding the directory path)
    my $tmp = $self->workdir;
    my $input = $tmp."/".$self->filename;
    $self->filename ($input);
    $tmp .= "/".$self->results;
    $self->results ($tmp);


    eval {
	$seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI");
    };        
    if (!$@) {
	#The inputId is a sequence file...got the normal way...

	# write sequence to file
	$self->writefile;        

	# run program
	$self->run_program;

	# parse output
	$self->parse_results;
	$self->deletefiles;
    }
    else {
	#The query object is not a seq object but a file.
	#Perhaps should check here or before if this file is fasta format...if not die
	#Here the file does not need to be created or deleted. Its already written and may be used by other runnables.

	$self->filename($self->query);

	# run program
	$self->run_program;

	# parse output
	$self->parse_results;
    }
    
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
    # run program
    # print STDERR "Running ".$self->program." ".$self->filename." -l > ".$self->results."\n";
    $self->throw ("Error running ".$self->program." on ".$self->filename) 
        unless ((system ($self->program." ".$self->filename." -l > ".$self->results)) == 0); 
}


=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses program output to give a set of features
 Example  :
 Returns  : 
 Args     : filename (optional, can be filename, filehandle or pipe, not implemented)
 Throws   :

=cut

sub parse_results {
    my ($self) = @_;
    my $filehandle;
    my $resfile = $self->results;
    
    if (-e $resfile) {
        # it's a filename
        if (-z $self->results) {  
	    # print STDERR $self->program." didn't find anything\n";
	    return;
        }       
        else {
            open (OUT, "<$resfile") or $self->throw ("Error opening $resfile");
            $filehandle = \*OUT;
      }
    }
    else {
        # it'a a filehandle
        $filehandle = $resfile;
    }
    
    # parse
    while (<$filehandle>) {
        chomp;
        next if /^$/;
        if (/^\>/) {
          #>(110-122) complexity=2.20 (12/2.20/2.50)
          /^\>(\S?)\((\d+)\-(\d+)\)\s*complexity=(\S+)/;

	  my $tid = $1;
          my $start = $2;
          my $end = $3;
          my $score = $4;

          my (%feature);
	  $feature{name} = $tid;
          $feature{score} = $score;
          $feature{start} = $start;
          $feature{end} = $end;
          ($feature{source}) = $self->program =~ /([^\/]+)$/;
          $feature{primary} = 'low_complexity';
          ($feature{program}) = $self->program =~ /([^\/]+)$/;
          $feature{logic_name} = 'low_complexity';
          $self->create_feature (\%feature);
          # print STDERR "have start ".$start." end ".$end."\n";
	}
    }
    close $filehandle;   
}


=head2 create_feature

 Title    : create_feature
 Usage    : $self->create_feature ($feature)
 Function : creates a Bio::EnsEMBL::SeqFeature object from %feature,
            and pushes it onto @{$self->{'_flist'}}
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub create_feature {
    my ($self, $feat) = @_;

    my $analysis = $self->analysis;


    # create feature object
    my $feat1 = Bio::EnsEMBL::SeqFeature->new ( -seqname     => $feat->{name},
						-start       => $feat->{start},
						-end         => $feat->{end},
						-score       => 0,
						-analysis    => $analysis,
						-percent_id => 0,
						-p_value => 0,
						); 



    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => 0,
					      -end => 0,
					      -analysis => $analysis,
					      -seqname => 'Seg');
    
    
    my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
                                                -feature2 => $feat2);

    if ($feature) {
      push (@{$self->{'_flist'}}, $feature);
    }
}


=head2 output

 Title    : output
 Usage    : $self->output
 Function : returns an array of feature objects
 Example  :
 Returns  : an array of Bio::EnsEMBL::SeqFeature objects
 Args     :
 Throws   :

=cut

sub output {
    my ($self) = @_;
    my @list = @{$self->{'_flist'}};
    return @{$self->{'_flist'}};
}

sub get_low_complexity_length {
    my ($self) = @_;

    if ($self->query->length > 0) {    
	my $lc_length = 0;
	
	foreach my $feat ($self->output) {
	    #print STDERR "Start ".$feat->start." End ".$feat->end."\n";
	    $lc_length += abs($feat->end - $feat->start) + 1;
	}
	
	my $low_complexity = ($lc_length)/($self->query->length);
	
	$low_complexity *= 100;
	
	return $low_complexity;
    }
    else {
	return 0;
    }
}

		
1;
