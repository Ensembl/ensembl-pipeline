#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Dust

=head1 SYNOPSIS

  # create and fill Bio::Seq object

  my $seqfile = '/my/dir/clone.seq';

  my $seq = Bio::SeqIO->new(
      -file => $seqfile,
      -fmt  => 'Fasta'
  )->next_seq;

  # create Bio::EnsEMBL::Pipeline::Runnable::Dust object

  my $dust = Bio::EnsEMBL::Pipeline::Runnable::Dust->new(
      -query => $seq
  );

  $dust->workdir($workdir);

  $dust->run;
  my @results = $dust->output;

=head1 DESCRIPTION

Dust takes a Bio::Seq (or Bio::PrimarySeq) object and runs dust
on it. The resulting output file is parsed to produce a set of features.

=head1 CONTACT

Mail to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Dust;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::RepeatConsensus;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Root;
use FileHandle;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

  Args       : various
  Description: Runnable constructor
  Returntype : Bio::EnsEMBL::Pipeline::Runnable::Dust
  Caller     : general

=cut


sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    $self->{'_rflist'}    = [];    # an array of Bio::EnsEMBL::RepeatFeatures
    $self->{'_query'}     = undef; # location of Bio::Seq object
    $self->{'_dust'}      = undef; # location of dust executable
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of dust

    my( $sequence, $dust, $level) = $self->_rearrange([qw(
        QUERY
        DUST
        LEVEL
    )], @args);

    $level ||= 20;
    $dust  ||= 'dust';

    $self->dust($self->find_executable($dust));

    $self->query($sequence) if ($sequence);
    $self->level($level)    if ($level);

    return $self;
}


=head2 query

  Arg [1]    : Bio::SeqI $query
  Description: accessor for query sequence
  Returntype : Bio::SeqI
  Exceptions : query not a Bio::PrimarySeqI or Bio::SeqI
  Caller     : general

=cut

sub query {
    my ($self, $seq) = @_;
    if ($seq) {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) {
            $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
        }
        $self->{'_query'} = $seq ;

        $self->filename($self->query->id.".$$.seq");
        $self->results($self->filename.".out");
    }
    return $self->{'_query'};
}


=head2 dust

  Arg [1]    : string $dust
  Description: accessor for name of dust binary
  Returntype : string
  Exceptions : binary not found, or not execuable
  Caller     : general

=cut

sub dust {
    my ($self, $file) = @_;

    if ($file) {
        $self->throw("dust not found at $file: $!\n") unless (-e $file);
        $self->{'_dust'} = $file;
        $self->throw("dust not executable: $!\n") unless (-x $file);
    }
    return $self->{'_dust'};
}


=head2 level

  Arg [1]    : int $dust_level
  Description: accessor for level of dusting (default 20)
  Returntype : int
  Exceptions : value not an integer
  Caller     : general

=cut

sub level {
    my ($self, $arg) = @_;

    if (defined $arg) {
        $self->{'_level'} = $arg;
        $self->throw("Invalid dust level [$arg]")
         unless $arg =~ /(^\d+$)/;
    }
    return $self->{'_level'};
}


=head2 run

  Arg [none] :
  Description: runs the Dust runnable
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run {
    my ($self) = @_;

    # check clone
    my $seq = $self->query || $self->throw("Clone required for dust\n");
    #set directory if provided
    $self->workdir('/tmp') unless $self->workdir;
    $self->checkdir;

    # reset filename and results as necessary
    $self->filename($self->workdir()."/".$self->filename);
    $self->results($self->workdir()."/".$self->results);

    # write sequence to file
    $self->writefile;
    $self->run_dust;

    # get output
    $self->parse_results;
    $self->deletefiles;
}


sub run_dust {
    my ($self) = @_;

    # dust takes one compulsory argument (fasta file)
    # default output (stdout) is a dusted fasta file
    # optional second argument is dust level (integer)
    # adding another argument (of any value) returns
    # results as a list of start/end pairs, rather than fasta
    # (N.B. you can't get this output without specifying the
    # dust level as well!)

    my $cmd = join(" ", $self->dust, $self->filename, $self->level, "1 >", $self->results);

    print "Running dust\n";
    $self->throw("Error running dust on " . $self->filename . "\n")
     if (system $cmd);
}


=head2 parse_results

  Arg [none] :
  Description: parses output of dust
               produces a set of Bio::EnsEMBL::RepeatFeature
               these can be retrieved using the "output" method (qv)
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub parse_results {
    my ($self) = @_;

    my $fh;
    my $results = $self->results();

    if (-f $results) {
        if (-z $results) {
            print "No dust found\n";
            return;
        }
        open $fh, "< $results";
    }
    else {
        $fh = $results;
    }

    while (<$fh>)
    {
        chomp;
        if (/(\d+)\.\.(\d+)/) {

            my ($start, $end) = ($1, $2);
            $start++;
            $end++;
            my $length = $end - $start + 1;

            my $rc = Bio::EnsEMBL::RepeatConsensus->new;
            $rc->name             ("dust");
            $rc->repeat_class     ("dust");
            $rc->repeat_consensus ("N");

            my $rf = Bio::EnsEMBL::RepeatFeature->new;
            $rf->score            (0);
            $rf->start            ($start);
            $rf->end              ($end);
            $rf->strand           (0);
            $rf->hstart           (1);
            $rf->hend             ($end - $start + 1);
            $rf->repeat_consensus ($rc);

            push @{$self->{'_rflist'}}, $rf;
        }
        else {
            $self->warn("Unexpected line in dust output [$_]")
        }
    }
    close $fh;
}


=head2 output

  Arg [none] :
  Description: returns output of running dust
  Returntype : @{Bio::EnsEMBL::RepeatFeature}
  Exceptions : none
  Caller     : general

=cut

sub output {
    my ($self) = @_;

    return @{$self->{'_rflist'}};
}

1;
