# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Seg

=head1 SYNOPSIS

my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                  -fmt => 'Fasta',
                                );
$seq = $seqstream->next_seq;

my $seg = Bio::EnsEMBL::Pipeline::Runnable::Seg->new ( -CLONE => $seq);
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

package Bio::EnsEMBL::Pipeline::Runnable::Seg;

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI );


=head2 new

 Title    : new
 Usage    : my $seg =  Bio::Worm::Pipeline::Runnable::Seg->new ( -program    => '/usr/local/pubseq/bin/seg',
                                                                 -clone      => $clone,
                                                                 -analysisid => 4,
                                                               );
 Function : initialises Seg object
 Returns  : a Seg object
 Args     : a Bio::Seq object, the path to the seg binaries and an analysisId (all optional)
 Throws   :

=cut

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new (@_);    
  
  $self->{'_flist'} = [];               # an array of Bio::SeqFeatures
  $self->{'_sequence'}  = undef;        # location of Bio::Seq object
  $self->{'_seg'} = undef;              # location of seg executable
  $self->{'_workdir'}   = undef;        # location of tmp directory
  $self->{'_filename'}  = undef;        # file to store Bio::Seq object
  $self->{'_results'}   = undef;        # file to store results of seg run
  $self->{'_protected'} = [];           # a list of files protected from deletion
  
  my ($clone, $seg, $analysisid) = $self->_rearrange([qw(CLONE 
					                 PROGRAM
                                                         ANALYSISID)], 
					              @args);
  
  $self->clone ($clone) if ($clone);       
  $self->analysisid ($analysisid) if ($analysisid);
  
  if ($seg) {   
      $self->seg ($seg); 
  }
  else {   
      $self->seg ($self->locate_executable('seg'));
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


=head2 seg

 Title    : seg
 Usage    : $self->seg ('/usr/local/pubseq/bin/seg');
 Function : get/set method for the path to the seg binaries
 Example  :
 Returns  : File path
 Args     : File path (optional)
 Throws   :

=cut

sub seg {
    my ($self, $location) = @_;
    if ($location) {
        unless (-e $location) {
            $self->throw ("seg not found at $location");
	}
        $self->{'_seg'} = $location ;
    }
    return $self->{'_seg'};
}

# analysis methods

=head2 run

 Title    : run
 Usage    : $self->run ($workdir, $args)
 Function : runs seg and populates @{$self->{'_flist'}} (array of features)
 Example  :
 Returns  :   
 Args     : workdir, args (optional)
 Throws   :

=cut

sub run {
    my ($self, $dir, $args) = @_;

    # nothing to be done with $args

    # check clone
    my $seq = $self->clone || $self->throw("Clone required for Seg\n");

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

    # run seg
    $self->run_seg;

    # parse output
    $self->parse_results;
    $self->deletefiles;
}


=head2 run_seg

 Title    : run_seg
 Usage    : $self->seg
 Function : makes the system call to seg
 Example  :
 Returns  : 
 Args     :
 Throws   :

=cut

sub run_seg {
    my ($self) = @_;
    # run seg
    print "running seg\n";
    $self->throw ("Error running seg on ".$self->filename) 
        if (system ($self->seg.' '.$self->filename." -l > ".$self->results)); 
}


=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses seg output to give a set of features
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
	    print "seg didn't find anything\n";
	    return;
        }       
        else {
            open (SEGOUT, "<$resfile") or $self->throw ("Error opening $resfile");
            $filehandle = *SEGOUT;
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
            /^\>\s*(\S+)\s*\((\d+)\-(\d+)\)\s*complexity=(\S+)/;
            my $id = $1;
            my $start = $2;
            my $end = $3;
            my $score = $4;
            my (%feature);
	    $feature {name} = $id;
     	    $feature {score} = $score;
	    $feature {start} = $start;
	    $feature {end} = $end;
            $feature {strand} = '.';
	    $feature {source}= 'seg';
	    $feature {primary}= 'low_complexity';
	    $feature {program} = 'seg';
  	    $self->create_feature(\%feature);
	}
    }
    close $filehandle;   
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
    my $analysis = Bio::EnsEMBL::Analysis->new ( -db              => $feat->{'db'},               # optional
                                                 -db_version      => $feat->{'db_version'},       # optional
                                                 -program         => $feat->{'program'},
                                                 -program_version => $feat->{'program_version'},
                                                 -gff_source      => $feat->{'source'},
                                                 -gff_feature     => $feat->{'primary'},
                                               );

    # create featurepair object
    my $feature1 = Bio::EnsEMBL::SeqFeature->new ( -seqname     => $feat->{'name'},
                                                   -start       => $feat->{'start'},
                                                   -end         => $feat->{'end'},
                                                   -score       => $feat->{'score'},
                                                   -p_value     => $feat->{'p_value'},            # optional
                                                   -percent_id  => $feat->{'percent_id'},
                                                   -source_tag  => $feat->{'source'},
                                                   -primary_tag => $feat->{'primary'},
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
