# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

 Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp

=head1 SYNOPSIS

 my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                   -fmt => 'Fasta',
                                 );
 $seq = $seqstream->next_seq;

 my $signalp = Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp->new ( -CLONE => $seq);
 $signalp->workdir ($workdir);
 $signalp->run;
 my @results = $signalp->output;

=head1 DESCRIPTION

 Signalp takes a Bio::Seq (or Bio::PrimarySeq) object
 and runs signalp on it (detecting signal peptides). 
 The resulting output file is parsed to produce a set of features.

=head1 CONTACT

 Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX
 
 The rest of the documentation details each of the object methods. 
 Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp;

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

 Title    : new
 Usage    : my $signalp =  Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp->new
                           ( -program    => '/usr/local/pubseq/bin/signalp',
                             -clone      => $clone,
                             -analysisid => 4,
                           );
 Function : initialises Signalp object
 Returns  : a Signalp object
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
  
    my ($clone, $program, $analysisid) = $self->_rearrange([qw(CLONE 
					                       PROGRAM
                                                               ANALYSISID)], 
					                    @args);

    print STDERR "PROGR: $program\n";
  
    $self->clone ($clone) if ($clone);       
    $self->analysisid ($analysisid) if ($analysisid);

    if (!defined $program) {
	$self->program ($self->find_executable ($program));
    }
    else {
	$self->program($program);
    }

    return $self;
}

###################
# get/set methods 
###################

=head2 clone

 Title    : clone
 Usage    : $self->clone ($clone);
 Function : get/set method for the Sequence object; assigns clone, filename
iprscan/bin/scanregexpf.pl       | /analysis/iprscan/data/confirm.patterns | NULL                                                                                                                and results
 Example  :
 Returns  : a Bio::Seq or Bio::PrimarySeq object
 Args     : a Bio::Seq or Bio::PrimarySeq object (optional)
 Throws   :

=cut

sub clone {
    my ($self, $seq) = @_;
    if ($seq) {
	if ($seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")) {
	    $self->{'_sequence'} = $seq ;
	    $self->clonename ($self->clone->id);
	    $self->filename ($self->clone->id.".$$.seq");
	    $self->results ($self->filename.".out");
	}
	else {
	    print STDERR "WARNING: The input_id is not a Seq object but if its a peptide fasta file, it should go fine\n";
	    $self->{'_sequence'} = $seq ;
	    $self->filename ("$$.tmp.seq");

	    print STDERR "FILENAMEN: ".$self->filename."\n";
	    $self->results ("sigp.$$.out");
	    
	}
    }
    return $self->{'_sequence'};
}


=head2 analysisid

 Title    : analysisid
 Usage    : $self->analysisid ($analysisid);
 Function : get/set method for the analysisId
 Example  :
 Returns  : analysisId
 Args     : analysisId (optional)
 Throws   :

=cut

sub analysisid {
    my $self = shift;
    if (@_) {
        $self->{'_analysisid'} = shift;
    }
    return $self->{'_analysisid'};
} 


=head2 program

 Title    : program
 Usage    : $self->program ('/usr/local/pubseq/bin/signalp');
 Function : get/set method for the path to the executable
 Example  :
 Returns  : File path
 Args     : File path (optional)
 Throws   :

=cut

sub program {
    my $self = shift;
    if (@_) {
        $self->{'_program'} = shift;
    }
    return $self->{'_program'};
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
 Args     : workdir (optional)
 Throws   :

=cut

sub run {
    my ($self, $dir) = @_;

    # check clone
    my $seq = $self->clone || $self->throw("Clone required for Program\n");

    # set directory if provided
    $self->workdir ('/tmp') unless ($self->workdir($dir));
    $self->checkdir;

    # reset filename and results as necessary (adding the directory path)
    my $tmp = $self->workdir;
    my $input = $tmp."/".$self->filename;
    print STDERR "INTPUT: $input\n";
    $self->filename ($input);
    $tmp .= "/".$self->results;
    $self->results ($tmp);

    if ($seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")) {
	
	# truncate the sequence (the first 50 aa are enough for signalp)
	my $sub_seq = substr ($seq->seq, 0, 50);
	$seq->seq ($sub_seq);
	
	# write sequence to file
	$self->writefile;        

	# run program
	$self->run_program;
	
	# parse output
	$self->parse_results;
	$self->deletefiles;
    }
    else {
	my $in  = Bio::SeqIO->new(-file => $seq, '-format' =>'Fasta');
	my $out = Bio::SeqIO->new(-file => ">".$self->filename.".cutted", '-format' =>'Fasta');

	open (OUT,">".$self->filename.".cutted");

	print STDERR "SEQ: ".$self->filename.".cutted\n";
	while ( my $tmpseq = $in->next_seq() ) {
	    print STDERR "AC: ".$tmpseq->display_id."\n";
	    #print STDERR "AC: ".$tmpseq->seq."\n";
	    
	    my $sub_seq = substr ($tmpseq->seq, 0, 50);
	    
	    $tmpseq->seq($sub_seq);

	    print OUT ">".$tmpseq->display_id."\n".$tmpseq->seq."\n";

	    #$out->write_seq($tmpseq);
	    
	    
	}
	close(OUT);
	$self->filename($self->filename.".cutted");
	$self->run_program;
	$self->parse_results;
	    #$self->clone ($tmpseq);
	    
	    # write sequence to file
	 #   $self->writefile;        
	    
	    # run program
	  #  $self->run_program;
	    
	    # parse output
	   # $self->parse_results;
	   # $self->deletefiles;
	    
	#}
	
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
    
    print STDERR "RUNNING: ".$self->program." -t euk ".$self->filename." > ".$self->results."\n";
    
    my $filename = $self->filename;

    $filename = "/tmp/1980455.tmp.seq.cutted";
    

    # run program
    print STDERR "running ".$self->program."\n";
    
    $self->throw ("Error running ".$self->program." on ".$self->filename) 
        unless ((system ($self->program." -t euk ".$self->filename. " > ".$self->results)) == 0); 
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
	    print STDERR $self->program." didn't find anything\n";
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
    my ($id, $fact1, $fact2, $end);
    while (<$filehandle>) {
        chomp;
        if (/^\>(\S+)/) {
            $id = $1;
	}
        elsif (/max\.\s+Y\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
            $fact1 = $2;
        }
        elsif (/mean\s+S\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
            $fact2 = $2;
            if ($fact1 eq "YES" && $fact2 eq "YES") {
                my $line = <$filehandle>;
                if ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
                    $end = $1;
		}
                else {
                    $self->throw ("parsing problem in ".$self->program);
		}
                my (%feature);
	        $feature{name} = $id;
       	        $feature{start} = 1;
	        $feature{end} = $end;
                ($feature{source}) = $self->program =~ /([^\/]+)$/;
	        $feature{primary}= 'signal_peptide';
	        ($feature{program}) = $self->program =~ /([^\/]+)$/;
                $feature{logic_name} = 'signal_peptide';
  	        $self->create_feature (\%feature);
	    }
        }
    }
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

    # create analysis object (will end up in the analysis table)
    my $analysis = Bio::EnsEMBL::Analysis->new ( -program         => $feat->{program},
                                                 -gff_source      => $feat->{source},
                                                 -gff_feature     => $feat->{primary},
                                                 -logic_name      => $feat->{logic_name},
                                               );

    # create feature object
    my $feat1 = Bio::EnsEMBL::SeqFeature->new ( -seqname     => $feat->{name},
						-start       => $feat->{start},
						-end         => $feat->{end},
						-source_tag  => $feat->{source},
						-primary_tag => $feat->{primary},
						-analysis    => $analysis,
						-percent_id => 'NULL',
						-p_value => 'NULL',
                                                 ); 



    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => 0,
					      -end => 0,
					      -analysis => $analysis,
					      -seqname => 'Sigp');
    
    
    my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
						-feature2 => $feat2);

    if ($feature) {
	#$feature->validate_prot_feature;
	# add to _flist
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

1;
