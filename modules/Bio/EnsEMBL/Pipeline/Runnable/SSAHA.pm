#
# Written by Simon Potter
#
# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::SSAHA

=head1 SYNOPSIS

#create and fill Bio::Seq object

my $clonefile = '/path/to/seq.fa';
my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
$seq = $seqstream->next_seq();

#create Bio::EnsEMBL::Pipeline::Runnable::SSAHA object

my $ssaha = Bio::EnsEMBL::Pipeline::Runnable::SSAHA->new(
    -CLONE => $seq,
    -DB    => 'mydb.fa'
);

$ssaha->workdir($workdir);
$ssaha->run();
my @results = $ssaha->output();

=head1 DESCRIPTION

SSAHA takes a Bio::Seq (or Bio::PrimarySeq) object and runs SSAHA
against a set of sequences.  The resulting output file is parsed
to produce a set of features.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::SSAHA;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_fplist'}    = [];    # an array of Bio::SeqFeatures
  $self->{'_sequence'}  = undef; # location of Bio::Seq object
  $self->{'_ssaha'}     = undef; # location of ssaha executable
  $self->{'_workdir'}   = undef; # location of temp directory
  $self->{'_filename'}  = undef; # file to store Bio::Seq object
  $self->{'_results'}   = undef; # file to store results of ssaha
  $self->{'_protected'} = [];    # a list of files protected from deletion ???
  $self->{'_min_length'}= 0;     # minimum length required for a
                                 # predicted match to be reported

  my ($seq, $len, $ssaha, $db) = $self->_rearrange([
    qw(CLONE LENGTH SSAHA DB)
  ], @args);

  $ssaha = 'ssaha' unless ($ssaha);
  $self->ssaha($self->find_executable($ssaha));

  $self->query($seq) if ($seq);
  $self->min_length($len) if ($len);
  $self->database($db) if ($db);

  return $self; # success - we hope!
}

#################
# get/set methods
#################
# really ough to be renamed "sequence" but this involves rewriting RunnableI::writefile and also any other modules that inherit from it.
# to do!
sub query {
  my ($self, $seq) = @_;
  if ($seq)
  {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI"))
      {
	  $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
      }
      $self->{'_sequence'} = $seq ;

      $self->filename($self->query->id.".$$.seq");
      $self->results($self->filename.".out");
  }
  return $self->{'_sequence'};
}


sub ssaha {
    my ($self, $location) = @_;
    if ($location) {
        $self->throw("SSAHA not found at $location: $!\n")
         unless (-e $location);
        $self->{'_ssaha'} = $location ;
    }
    return $self->{'_ssaha'};
}


sub min_length {
    my ($self, $len) = @_;
    if (defined $len && $len > 0) {
      $self->{'_min_length'} = $len ;
    }
    return $self->{'_min_length'};
}


sub database {
    my ($self, $db) = @_;
    if ($db) {
      $self->{'_database'} = $db ;
    }
    return $self->{'_database'};
}


###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run($workdir, $args)
    Function:   Runs ssaha script and creates array of features
    Returns :   none
    Args    :   optional $workdir and $args (e.g. '-ace' for ace file output)

=cut

sub run {
    my ($self) = @_;
    #check clone
    my $seq = $self->query() || $self->throw("Clone required for ssaha\n");
    #set directory if provided
    $self->workdir('/tmp') unless ($self->workdir());
    $self->checkdir();

    # reset filename and results as necessary
    $self->filename($self->workdir()."/".$self->filename());
    $self->results($self->workdir()."/".$self->results());

    #write sequence to file
    $self->writefile();
    $self->run_ssaha();
    #parse output of ssaha
    $self->parse_results();
    $self->deletefiles();

    1;
}

sub run_ssaha {
    my ($self) = @_;
    #run ssaha

    # TODO - need to deal with databases with pre-computed
    # hash tables
    my @args = (
	"-pf",
	"-rq",
        "-qf fasta",
	"-sf fasta",
	"-lm null"
    );

    push @args, "-mp " . $self->min_length if $self->min_length;
    my $cmd = join " ", $self->ssaha, $self->filename, $self->database, @args;
#    print STDERR "Running ssaha: $cmd\n";

    $self->throw("Error running ssaha on ".$self->filename."\n")
        if (system ($cmd . " > " . $self->results));
}

=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses ssaha output to give a set of features
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut
sub parse_results {
    my ($self, $filehandle) = @_;
    my $resfile = $self->results();

    if (-e $resfile) {
        if (-z $self->results) {
	    print STDERR "SSAHA didn't find any matches\n";
	    return; 
        } else {
            open (OUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");
            $filehandle = \*OUT;
        }
    }
    else { #it'a a filehandle
        $filehandle = $resfile;
    }

    #extract values
    while (<$filehandle>)
    {
        if (/^[FR]/) # ignore lines before matches
        {
	    my (
	       $strand, $query, $qstart, $qend,
	       $sbjct, $sstart, $send, $match, $pc
	    ) = split;

            my (%feat1, %feat2);

            $feat1 {name} = $query;
            $feat2 {name} = $sbjct;

            $feat1 {start} = $qstart;
            $feat1 {end}   = $qend;

            $feat2 {start} = $sstart;
            $feat2 {end}   = $send;

            $feat1 {score} = $match;
            $feat2 {score} = $feat1 {score};

            $feat1 {percent} = $pc;
            $feat2 {percent} = $feat1 {percent};

            if ($strand eq 'F') {
                $feat2 {strand} = 1;
            } elsif ($strand eq 'R') {
                $feat2 {strand} = -1 ;
            }
            $feat1 {strand} = $feat2 {strand};

            $feat2 {db} = undef;
            $feat2 {db_version} = undef;
            $feat2 {program} = 'ssaha';
            $feat2 {p_version}='1';
            $feat2 {source}= 'ssaha';
            $feat2 {primary}= 'similarity';
            $feat1 {source}= 'ssaha';
            $feat1 {primary}= 'similarity';

            $self->create_FeaturePair(\%feat1, \%feat2);

        }
    }
    close $filehandle;
}


##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of features
    Returns :   Returns an array of features
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_fplist'}};
}

1;

