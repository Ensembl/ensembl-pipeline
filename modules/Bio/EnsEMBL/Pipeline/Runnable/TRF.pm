#
#
# Written by Simon Potter
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::TRF

=head1 SYNOPSIS

  # create and fill Bio::Seq object
  my $seqfile = '/path/to/seq.fa';
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(
    -file => $seqfile,
    -fmt => 'Fasta'
  );
  $seq = $seqstream->next_seq();

  # create Bio::EnsEMBL::Pipeline::Runnable::TRF object
  my $trf = Bio::EnsEMBL::Pipeline::Runnable::TRF->new(
    -CLONE => $seq
  );
  $trf->workdir($workdir);
  $trf->run();

  # get results
  my @results = $trf->output();

=head1 DESCRIPTION

TRF takes a Bio::Seq (or Bio::PrimarySeq) object and runs TRF on it. The
resulting .dat file is parsed to produce a set of feature pairs.

=head2 Methods:

=over4

=item new($seq_obj)

=item trf($path_to_TRF)

=item workdir($directory_name)

=item run()

=item output()

=back

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::TRF;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use FileHandle;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    $self->{'_fplist'}    = [];     # an array of feature pairs
    $self->{'_clone'}     = undef;  # location of Bio::Seq object
    $self->{'_trf'}       = undef;  # location of TRF executable
    $self->{'_workdir'}   = undef;  # location of temp directory
    $self->{'_filename'}  = undef;  # file to store Bio::Seq object
    $self->{'_results'}   = undef;  # file to store results of TRF
    $self->{'_protected'} = [];     # a list of files protected from deletion

    # TRF takes seven parameters - defaults ripped from a humpub script
    # see accessor methods for description
    $self->{'_match'}     = 2;
    $self->{'_mismatch'}  = 7;
    $self->{'_delta'}     = 7;
    $self->{'_pm'}        = 80;
    $self->{'_pi'}        = 10;
    $self->{'_minscore'}  = 20;
    $self->{'_maxperiod'} = 500; 

    my ($clone, $trf, $match, $mismatch, $delta, $pm, $pi, $minscore, $maxperiod) = $self->_rearrange([qw(
	CLONE
	TRF
        MATCH
        MISMATCH
        DELTA
        PM
        PI
        MINSCORE
        MAXPERIOD
    )], @args);
    

    $trf = 'trf301' unless defined($trf);
    $trf = $self->find_executable($trf);
    $self->trf($trf);

    $self->match     ($match)     if defined $match;
    $self->mismatch  ($mismatch)  if defined $mismatch;
    $self->delta     ($delta)     if defined $delta;
    $self->pm        ($pm)        if defined $pm;
    $self->pi        ($pi)        if defined $pi;
    $self->minscore  ($minscore)  if defined $minscore;
    $self->maxperiod ($maxperiod) if defined $maxperiod;

    $self->clone($clone) if ($clone);       

    return $self;
}

#################
# get/set methods 
#################

sub clone {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
        }
        $self->{'_clone'} = $seq ;

	my @arg_list = (
            $self->match,
            $self->mismatch,
            $self->delta,
            $self->pm,
            $self->pi,
            $self->minscore,
            $self->maxperiod
	);
	$self->options(join(' ', @arg_list));
        
        $self->clonename($self->clone->id);

	# TRF seems to want to truncate the filename to 13 chars ...
        $self->filename(substr($self->clone->id.".$$", 0, 13));

        # output has (equally annoyingly) parameters in the filename
        $self->results($self->filename . "." . join('.', @arg_list) . '.dat');
    }
    return $self->{'_clone'};
}

=head2 protect

=cut

=head2 trf

=cut

sub trf {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("TRF not found at $location: $!\n") 
         unless -e $location;
        $self->{_trf} = $location ;
    }
    return $self->{_trf};
}

=head2 workdir

=cut

=head2 clonename

=cut

###########
# Analysis methods
##########

=head2 run

=cut

sub run {
    my ($self) = @_;

    # check seq
    my $seq = $self->clone() || $self->throw("Seq required for TRF\n");

    # set directory if provided
    $self->workdir('/tmp') unless $self->workdir();

    $self->checkdir();

    # write sequence to file
    $self->writefile();        
    $self->run_trf();

    # parse output of trf
    $self->parse_results();
    $self->deletefiles();
}

=head2 parsefile

=cut

sub run_trf {
    my ($self) = @_;

    # options() set in clone() at same time as results()
    # both are a catenation of the parameters

    print "Running TRF; ";
    my $cmd = join(" ", $self->trf, $self->filename, $self->options, '-d');
    print "command: $cmd\n";

    # Don't test return value of system() as TRF exits with
    # non-zero status (at time of writing)
    system($cmd);
}

#New and improved! takes filenames and handles, therefore pipe compliant!
sub parse_results {
    my ($self) = @_;
    my $filehandle;

    my $results = $self->results;
    if (ref ($results) !~ /GLOB/)
    {
        $filehandle = new FileHandle;
        $filehandle->open("< $results")
	 or $self->throw("Error opening $results \n");
    }
    else
    {
        $filehandle = $results;
    } 
    
    my $found = 0;
    my $seqname;
    while (<$filehandle>)
    {  
	if (/Sequence: (\w+)/) {
	    $seqname = $1;
	}
	next unless (/^\d/); # ignore introductory lines
	$found++;

        my (
            $start,          $end,            $period_size,
            $copy_number,    $consensus_size, $percent_matches,
            $percent_indels, $score,          $A,
            $T,              $G,              $C,
            $entropy,        $mer
        ) = split;

        my (%feat1, %feat2);

        $feat1 {name} = $mer;
        $feat1 {score} = $score;
        $feat1 {start} = $start;
        $feat1 {end} = $end;
        $feat1 {strand} = 0;
        $feat1 {source}= 'trf';
        $feat1 {primary}= 'repeat';

        $feat2 {name} = $mer;
        $feat2 {score} = $score;
        $feat2 {start} = $start;
        $feat2 {end} = $end;
        $feat2 {strand} = 0;

        $feat2 {db} = undef;
        $feat2 {db_version} = undef;
        $feat2 {program} = 'trf';
        $feat2 {p_version}='1';
        $feat2 {source}= 'trf';
        $feat2 {primary}= 'repeat';

        $self->create_repeat(\%feat1, \%feat2);

    }
    print STDERR "No tandem repeats found\n" unless $found;
    close $filehandle;   
}


##############
# input/output methods
#############

=head2 output

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_fplist'}};
}

# Accessor methods

# TODO: proper doc for these methods
# The following lines ripped straight from trf

# Please use: trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]
# Where: (all scores are positive)
#   File = sequences input file
#   Match  = matching weight
#   Mismatch  = mismatching penalty
#   Delta = indel penalty
#   PM = match probability (whole number)
#   PI = indel probability (whole number)
#   Minscore = minimum alignment score to report
#   MaxPeriod = maximum period size to report
#   [options] = one or more of the following :
#                -m    masked sequence file
#                -f    flanking sequence

sub match {
    my ($self, $arg) = @_;

    $self->{'_match'} = $arg if defined $arg;

    return $self->{'_match'};
}

sub mismatch {
    my ($self, $arg) = @_;

    $self->{'_mismatch'} = $arg if defined $arg;

    return $self->{'_mismatch'};
}

sub delta {
    my ($self, $arg) = @_;

    $self->{'_delta'} = $arg if defined $arg;

    return $self->{'_delta'};
}

sub pm {
    my ($self, $arg) = @_;

    $self->{'_pm'} = $arg if defined $arg;

    return $self->{'_pm'};
}

sub pi {
    my ($self, $arg) = @_;

    $self->{'_pi'} = $arg if defined $arg;

    return $self->{'_pi'};
}

sub minscore {
    my ($self, $arg) = @_;

    $self->{'_minscore'} = $arg if defined $arg;

    return $self->{'_minscore'};
}

sub maxperiod {
    my ($self, $arg) = @_;

    $self->{'_maxperiod'} = $arg if defined $arg;

    return $self->{'_maxperiod'};
}

sub options {
    my ($self, $arg) = @_;

    $self->{'_options'} = $arg if defined $arg;

    return $self->{'_options'};
}

1;
