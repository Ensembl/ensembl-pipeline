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

Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam

=head1 SYNOPSIS

  # something like this
  my $query = new Bio::Seq(-file   => $queryfile,
			   -format => 'Fasta');

  my $hmm =  Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam->new 
    ('-query'          => $query,
     '-program'        => 'hmmpfam' or '/usr/local/pubseq/bin/hmmpfam',
     '-database'       => 'Pfam');

  $hmm->workdir ($workdir);
  $hmm->run;
  my @results = $hmm->output;

=head1 DESCRIPTION

  Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs hmmpfam.
  The resulting output file is parsed to produce a set of Bio::EnsEMBL::FeaturePairs.

=head1 CONTACT

   Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

 Title    : new
 Usage    : my $hmm =  Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam->new ()
 Function : initialises Hmmpfam object
 Returns  : a Hmmpfam object
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
    $self->{'_options'}   = undef;        # additional arguments for blast
    $self->{'_workdir'}   = undef;        # location of tmp directory
    $self->{'_filename'}  = undef;        # file to store Bio::Seq object
    $self->{'_results'}   = undef;        # file to store results of seg run
    $self->{'_protected'} = [];           # a list of files protected from deletion
  
    my ($query, $analysis, $options) = $self->_rearrange([qw(QUERY 
                                                             ANALYSIS
                                                             OPTIONS)],        
                                                          @args);


    
    $self->query ($query) if ($query);       
    $self->analysis ($analysis) if ($analysis);
	
    $self->program ($self->analysis->program_file);
  
    if ($self->analysis->db_file) {
        $self->database($self->analysis->db_file);
    } else {
        $self->throw("Hmmpfam needs a database");
    }
    
    $self->options ($options) if ($options);

    return $self;
}


###################
# get/set methods 
###################

=head2 query

 Title    : query
 Usage    : $self->query ($query);
 Function : get/set method for the Sequence object; assigns query,
            seq-filename and result-filename
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
	    $self->queryname ($self->query);
	    $self->filename ($self->query.".$$.seq");
	    $self->lsresults ($self->filename.".lsout");
	    $self->fsresults ($self->filename.".fsout");
	}
	else {
	    print STDERR "WARNING: The input_id is not a Seq object but if its a peptide fasta file, it should go fine\n";
	    $self->{'_sequence'} = $seq ;
	    $self->filename ("$$.tmp.seq");
	    
	    $self->lsresults ("hmm.$$.lsout");
	     $self->fsresults ("hmm.$$.fsout");
	}
    }
    return $self->{'_sequence'};
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
        ($self->{'_analysis'}->isa ("Bio::EnsEMBL::Analysis") || $self->{'_analysis'}->isa ("Bio::EnsEMBL::Analysis"))             || $self->throw("Input isn't a Bio::EnsEMBL::Analysis");
    }
    return $self->{'_analysis'};
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
 my ($self, $dir) = @_;

    # check query
    my $seq = $self->query || $self->throw("Query required for Program");

    # set directory if provided
    $self->workdir ('/tmp') unless ($self->workdir($dir));
    $self->checkdir;

    # reset filename and results as necessary (adding the directory path)
    my $tmp1 = $self->workdir;
    my $tmp2 = $self->workdir;
 my $input = $tmp1."/".$self->filename;
    $self->filename ($input);
    $tmp1 .= "/".$self->lsresults;
    $self->lsresults ($tmp1);
    $tmp2 .= "/".$self->fsresults;
    $self->fsresults ($tmp2);


 eval {
     $seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")
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
    print STDERR "running ".$self->program." against ".$self->database."\n";

    # some of these options require HMMER 2.2g (August 2001)
    
    my @dbfiles = split(/;/,$self->analysis->db_file);

    if ($dbfiles[0] =~ /ls/) {
	 

	my $cmd = $self->program .' '.
	        '--acc --cut_ga --cpu 1 '.
	        $self->options .' '.
	        $dbfiles[0]      .' '.
	        $self->filename.' > '.
		$self->lsresults;
    print STDERR "$cmd\n";   
    $self->throw ("Error running ".$self->program." on ".$self->filename." against ".$dbfiles[0]) 
        unless ((system ($cmd)) == 0);
    }
    else {
	die || "ls pfam file has not been provided";
    }

   if ($dbfiles[1] =~ /fs/) { 
       
	my $cmd = $self->program .' '.
	        '--acc --cut_ga --cpu 1 '.
	        $self->options .' '.
	        $dbfiles[1]      .' '.
	        $self->filename.' > '.
		$self->fsresults;
	print STDERR "$cmd\n";   
	$self->throw ("Error running ".$self->program." on ".$self->filename." against ".$dbfiles[1]) 
	    unless ((system ($cmd)) == 0);
    }
    else {
	die || "fs pfam file has not been provided";
    }
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
    my $fshandle;
    my $id;
    my $resfile = $self->lsresults;
    my $fsfile  = $self->fsresults;

    #$resfile = "/tmp/hmm.4063364.lsout";
    #$fsfile = "/tmp/hmm.4063364.fsout";
		

    if (-e $resfile) {
        # it's a filename
        if ((-z $self->lsresults) && (-z $self->fsresults)) {  
            print STDERR $self->program." didn't find anything\n";
            return;
        }       
        else {
	    open (OUT, "<$resfile") or $self->throw ("Error opening $resfile");
            $filehandle = \*OUT;
	

	    open (OUT1, "<$fsfile")  or $self->throw ("Error opening $fsfile");
	    $fshandle = \*OUT1;
	}
    }
    else {
        # it'a a filehandle
        $filehandle = $resfile;
	$fshandle = $fsfile;
    }
   


#First parse what comes from the ls mode matches. Every match in that case is taken
    while (<$filehandle>) {
	chomp;
        last if /^Alignments of top-scoring domains/;
        next if (/^Model/ || /^\-/ || /^$/);
        if (/^Query sequence:\s+(\S+)/) {
            $id = $1;
	}
	
        if (my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
	
            my %feature;
            ($feature{name}) = $id;
            $feature{score} = $score;
            $feature{p_value} = sprintf ("%.3e", $evalue);
            $feature{start} = $start;
            $feature{end} = $end;
            $feature{hname} = $hid;
            $feature{hstart} = $hstart;
            $feature{hend} = $hend;
            ($feature{source}) = $self->program =~ /([^\/]+)$/;
            $feature{primary} = 'similarity';
            ($feature{program}) = $self->program =~ /([^\/]+)$/;
            ($feature{db}) = $self->database =~ /([^\/]+)$/;
            ($feature{logic_name}) = $self->program =~ /([^\/]+)$/;
            $self->create_feature (\%feature);
	}
    }
    close FILEHANDLE; 

#Then read all of the fs mode matches. If a match does not overlap with any ls match thus its taken
    while (<$fshandle>) {
	my ($hid, $start, $end, $hstart, $hend, $score, $evalue);
	my $overlap = undef;
        
	chomp;
	
        last if /^Alignments of top-scoring domains/;
        next if (/^Model/ || /^\-/ || /^$/);
        if (/^Query sequence:\s+(\S+)/) {
            $id = $1;
        }
	if (($hid, $start, $end, $hstart, $hend, $score, $evalue) = /^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
	   
	
	    foreach my $featpair(@{$self->{'_flist'}}) {
		my $lsstart = $featpair->feature1->start;
		my $lsend = $featpair->feature1->end;
		if ((($start >= $lsstart) && ($start <= $lsend)) || (($end >= $lsstart) && ($end <= $lsend))) {
		    $overlap = 1;
		}
	    }

	    if (!defined $overlap) {
		#There is no ovelap thus create a feature
		my %feature;
		($feature{name}) = $id;
		$feature{score} = $score;
		$feature{p_value} = sprintf ("%.3e", $evalue);
		$feature{start} = $start;
		$feature{end} = $end;
		$feature{hname} = $hid;
		$feature{hstart} = $hstart;
		$feature{hend} = $hend;
		($feature{source}) = $self->program =~ /([^\/]+)$/;
		$feature{primary} = 'similarity';
		($feature{program}) = $self->program =~ /([^\/]+)$/;
		($feature{db}) = $self->database =~ /([^\/]+)$/;
		($feature{logic_name}) = $self->program =~ /([^\/]+)$/;
		$self->create_feature (\%feature);
	    }
	}
    }
    close (FS);
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

    # create featurepair object
    my $feature1 = Bio::EnsEMBL::SeqFeature->new;
    $feature1->seqname ($feat->{name});
    $feature1->start ($feat->{start});
    $feature1->end ($feat->{end});
    $feature1->score ($feat->{score});
    $feature1->p_value ($feat->{p_value});
    $feature1->source_tag ($feat->{source});
    $feature1->primary_tag ($feat->{primary});
    $feature1->analysis ($self->analysis);

    my $feature2 = Bio::EnsEMBL::SeqFeature->new;
    $feature2->seqname ($feat->{hname});
    $feature2->start ($feat->{hstart});
    $feature2->end ($feat->{hend});
    $feature2->score ($feat->{score});
    $feature2->p_value ($feat->{p_value});
    $feature2->source_tag ($feat->{source});
    $feature2->primary_tag ($feat->{primary});
    $feature2->analysis ($self->analysis);

    my $featurepair = Bio::EnsEMBL::FeaturePair->new;
    $featurepair->feature1 ($feature1);
    $featurepair->feature2 ($feature2);

    if ($featurepair) {
        $featurepair->feature1->validate_prot_feature;
        $featurepair->feature2->validate_prot_feature;
        # add to _flist
        push (@{$self->{'_flist'}}, $featurepair);
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


sub lsresults {
    my ($self, $results) = @_;
    $self->{_lsresults} = $results if ($results);
    return $self->{_lsresults};
}


sub fsresults {
    my ($self, $results) = @_;
    $self->{_fsresults} = $results if ($results);
    return $self->{_fsresults};
}

1;
