# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Signalp

=head1 SYNOPSIS

my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                  -fmt => 'Fasta',
                                );
$seq = $seqstream->next_seq;

my $signalp = Bio::EnsEMBL::Pipeline::Runnable::Signalp->new ( -CLONE => $seq);
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

package Bio::EnsEMBL::Pipeline::Runnable::Signalp;

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

 Title    : new
 Usage    : my $signalp =  Bio::EnsEMBL::Pipeline::Runnable::Signalp->new
                           ( -program    => '/usr/local/pubseq/bin/signalp',
                             -clone      => $clone,
                             -analysisid => 4,
                           );
 Function : initialises Signalp object
 Returns  : a Signalp object
 Args     : a Bio::Seq object, the path to the signalp binaries and an analysisId (all optional)
 Throws   :

=cut

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new (@_);    
  
    $self->{'_flist'} = [];               # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef;        # location of Bio::Seq object
    $self->{'_signalp'} = undef;          # location of signalp executable
    $self->{'_workdir'}   = undef;        # location of tmp directory
    $self->{'_filename'}  = undef;        # file to store Bio::Seq object
    $self->{'_results'}   = undef;        # file to store results of signalp run
    $self->{'_protected'} = [];           # a list of files protected from deletion
  
    my ($clone, $signalp, $analysisid) = $self->_rearrange([qw(CLONE 
					                       PROGRAM
                                                               ANALYSISID)], 
					                    @args);
  
    $self->clone ($clone) if ($clone);       
    $self->analysisid ($analysisid) if ($analysisid);
  
    my $bindir = $::pipeConf{'bindir'} || undef;

    if (-x $signalp) { 
        # passed from RunnableDB (full path assumed)  
        $self->signalp ($signalp); 
    }
    elsif (defined $bindir && -x ($signalp = "$bindir/signalp")) {
        $self->signalp ($signalp);
    }
    else {   
        # search shell $PATH
        $self->signalp ($self->locate_executable('signalp'));
    }
  
    return $self;
}

###################
# get/set methods 
###################

=head2 clone

 Title    : clone
 Usage    : $self->clone ($clone);
 Function : get/set method for the Sequence object; assigns clone, filename and results
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


=head2 signalp

 Title    : signalp
 Usage    : $self->signalp ('/usr/local/pubseq/bin/signalp');
 Function : get/set method for the path to the signalp binaries
 Example  :
 Returns  : File path
 Args     : File path (optional)
 Throws   :

=cut

sub signalp {
    my ($self, $location) = @_;
    if ($location) {
        unless (-x $location) {
            $self->throw ("signalp not found at $location");
	}
        $self->{'_signalp'} = $location ;
    }
    return $self->{'_signalp'};
}

###################
# analysis methods
###################

=head2 run

 Title    : run
 Usage    : $self->run ($workdir, $args)
 Function : runs signalp and populates @{$self->{'_flist'}} (array of features)
 Example  :
 Returns  :   
 Args     : workdir, args (optional)
 Throws   :

=cut

sub run {
    my ($self, $dir, $args) = @_;

    # nothing to be done with $args

    # check clone
    my $seq = $self->clone || $self->throw("Clone required for Signalp\n");

    # set directory if provided
    $self->workdir ('/tmp') unless ($self->workdir($dir));
    $self->checkdir;

    # reset filename and results as necessary (adding the directory path)
    my $tmp = $self->workdir;
    my $input = $tmp."/".$self->filename;
    $self->filename ($input);
    $tmp .= "/".$self->results;
    $self->results ($tmp);

    # truncate the sequence (the first 50 aa are enough for signalp)
    my $sub_seq = substr ($seq->seq, 0, 50);
    $seq->seq ($sub_seq);

    # write sequence to file
    $self->writefile;        

    # run signalp
    $self->run_signalp;

    # parse output
    $self->parse_results;
    $self->deletefiles;
}


=head2 run_signalp

 Title    : run_signalp
 Usage    : $self->signalp
 Function : makes the system call to signalp
 Example  :
 Returns  : 
 Args     :
 Throws   :

=cut

sub run_signalp {
    my ($self) = @_;
    # run signalp
    print STDERR "running signalp\n";
    $self->throw ("Error running signalp on ".$self->filename) 
        unless ((system ($self->signalp." -t euk ".$self->filename." > ".$self->results)) == 0); 
}


=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses signalp output to give a set of features
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
	    print STDERR "signalp didn't find anything\n";
	    return;
        }       
        else {
            open (SIGNALPOUT, "<$resfile") or $self->throw ("Error opening $resfile");
            $filehandle = \*SIGNALPOUT;
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
                    $self->throw ("problem in signalp");
		}
                my (%feature);
	        $feature {name} = $id;
       	        $feature {start} = 1;
	        $feature {end} = $end;
 	        $feature {source}= 'signalp';
	        $feature {primary}= 'signal_peptide';
	        $feature {program} = 'signalp';
  	        $self->create_feature (\%feature);
	    }
        }
    }
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
    my $analysis = Bio::EnsEMBL::Analysis->new ( -db              => $feat->{db},               # optional
                                                 -db_version      => $feat->{db_version},       # optional
                                                 -program         => $feat->{program},
                                                 -program_version => $feat->{program_version},
                                                 -gff_source      => $feat->{source},
                                                 -gff_feature     => $feat->{primary},
                                               );

    # create featurepair object
    my $feature1 = Bio::EnsEMBL::SeqFeature->new ( -seqname     => $feat->{name},
                                                   -start       => $feat->{start},
                                                   -end         => $feat->{end},
                                                   -score       => $feat->{score},
                                                   -p_value     => $feat->{p_value},            # optional
                                                   -percent_id  => $feat->{percent_id},
                                                   -source_tag  => $feat->{source},
                                                   -primary_tag => $feat->{primary},
                                                   -analysis    => $analysis,
                                                 ); 
    $feature1->add_tag_value ('analysisid', $self->analysisid);

    my $feature2 = Bio::EnsEMBL::SeqFeature->new;

    my $featurepair = Bio::EnsEMBL::FeaturePair->new ( -feature1 => $feature1,
                                                       -feature2 => $feature2,
						     );

    if ($featurepair) {
	$featurepair->feature1->validate_prot_feature;
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
