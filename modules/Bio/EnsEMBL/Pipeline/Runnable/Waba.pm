# Author: Marc Sohrmann (ms2@sanger.ac.uk)
#
# You may distribute this module under the same terms as perl itself

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::Runnable::Waba

=head1 SYNOPSIS

  my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                    -fmt => 'Fasta',
                                  );
  $seq = $seqstream->next_seq;

  my $waba = Bio::EnsEMBL::Pipeline::Runnable::Waba->new ( -QUERY => $seq);
  $waba->workdir ($workdir);
  $waba->run;
  my @results = $waba->output;

=head1 DESCRIPTION

  Waba takes a Bio::Seq (or Bio::PrimarySeq) object
  and runs waba against a fasta database (Jim Kents aligner,
  initially written for elegans-briggsae comparison)
  The resulting output file is parsed to produce a set of features.

=head1 CONTACT

  B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Waba;

use vars qw(@ISA);
use strict;
$| = 1;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::SearchIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

    Arg [1]    : Bio::Seq $query
    Arg [2]    : char $waba_binary
    Arg [3]    : Bio::EnsEMBL::Analysis $analysis
    Description: initialises Waba object
    Returns    : a Waba object
    Exceptions : none
    Caller     : general

=cut

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new (@_);

    $self->{'_fplist'} = [];              # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef;        # location of Bio::Seq object
    $self->{'_waba'}      = undef;        # location of waba executable
    $self->{'_database'}  = undef;        # name of database
    $self->{'_workdir'}   = undef;        # location of tmp directory
    $self->{'_filename'}  = undef;        # file to store Bio::Seq object
    $self->{'_results'}   = undef;        # file to store results of waba run
    $self->{'_protected'} = [];           # a list of files protected from deletion

    my ($query, $waba, $database) = $self->_rearrange([qw(
	QUERY
        WABA
	DATABASE
    )], @args);

    $waba ||= 'waba';
    $self->query   ($query)    if ($query);
    $self->database($database) if ($database);
    $self->waba($self->find_executable($waba));

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
	($seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI"))
	    || $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
	$self->{'_sequence'} = $seq;
	$self->filename ($self->query->id.".$$.seq");
	$self->results ($self->filename.".out");
    }
    return $self->{'_sequence'};
}


=head2 waba

    Arg [1]    : string $waba
    Description: accessor for waba binary
    Returntype : string
    Exceptions : binary not found and executable
    Caller     : general

=cut

sub waba {
    my ($self, $location) = @_;

    if ($location) {
        unless (-x $location) {
            $self->throw ("program not found at $location");
	}
        $self->{'_waba'} = $location;
    }
    return $self->{'_waba'};
}


=head2 database

    Arg [1]    : string $db
    Description: accessor for target database
    Returntype : string
    Exceptions : none
    Caller     : general

=cut

sub database {
    my $self = shift;
    if (@_) {
        $self->{'_database'} = shift;
    }
    return $self->{'_database'};
}


=head2 run

    Args       : none
    Description: runs the waba program and creates set of features
    Returntype : none
    Exceptions : none
    Caller     : general

=cut

sub run {
    my ($self, $dir) = @_;

    # check clone
    my $seq = $self->query || $self->throw("Clone required for Program\n");

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

    Args       : none
    Description: makes the system call to program
    Returntype : none
    Exceptions : system call unsuccessful
    Caller     : general

=cut

sub run_program {
    my ($self) = @_;

    # write the files containing the subject file lists

    my $path = $self->filename;
    open (QP, ">$path.query") || $self->throw("CANNOT create small program input file");
    print QP "$path\n";
    close QP;

    # write the files containing the query file lists
    # check database is a fasta file (first char '>') or a list of files

    my $db = $self->database;
    open DB, "< $db";
    my $char = getc DB;
    close DB;

    if ($char eq '>') {
	$db = $self->workdir . "/waba.$$.db";
        open DB, "> $db" || $self->throw("cannot create large program input file");
	print DB $self->database, "\n";
        close DB;
    }

    my $cmd = $self->waba . " all $path.query $db " . $self->results .
     "> /dev/null";

    $self->throw("Error running ".$self->waba." on ".$self->filename)
        unless ((system($cmd)) == 0);

    unlink "$path.query";
    unlink "$db" if $db ne $self->database;
}


=head2 parse_results

    Arg [1]    : filehandle containing results file (optional)
    Description: parses output to give a set of features
    Returntype : none
    Exceptions : none
    Caller     : general

=cut

sub parse_results {
    my ($self) = @_;
    my $results = $self->results;
    my $waba;

    my $hname = $self->database;
    $hname =~ s{.*/}{};

    if (-e $results) {
        # it's a filename
        if (-z $self->results) {
	    print STDERR $self->waba." didn't find anything\n";
	    return;
        }
	$waba = Bio::SearchIO->new(
	    -format => 'waba',
	    -file   => $results
	);
    }
    else {
        # it'a a filehandle
	$waba = Bio::SearchIO->new(
	    -format => 'waba',
            -fh     => $results
	);
    }

    while (my $result = $waba->next_result) {
       while (my $hit = $result->next_hit) {
          while (my $hsp = $hit->next_hsp) {
	      $self->split_HSP($hsp, $hit->name);
	  }
       }
    }
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

    return @{$self->{'_fplist'}};
}


=head2 split_HSP

    Arg [1]    : Bio::Search::HSP::WABAHSP
    Arg [2]    : feature name
    Description: splits an HSP into ungapped FeaturePairs
    Returntype : none
    Exceptions : none
    Caller     : general

=cut

sub split_HSP {
    my ($self, $hsp, $name) = @_;

    my $qstrand = $hsp->query->strand;
    my $hstrand = $hsp->subject->strand;

    # We split the alignment strings into arrays of one char each.
    # We then loop over this array and when we come to a gap
    # in either the query sequence or the hit sequence we make
    # a new feature pair.

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

    my @qchars = split(//,$hsp->query_string);  # split alignment into array of char
    my @hchars = split(//,$hsp->hit_string);    # ditto for hit sequence

    my ($qstart, $hstart, $qend, $hend);

    if ($qstrand == 1) {
        $qstart = $hsp->query->start();   # Start off the feature pair start
    }
    else {
        $qstart = $hsp->query->end;
    }

    if ($hstrand == 1) {
        $hstart = $hsp->subject->start();   # ditto
    }
    else {
        $hstart = $hsp->subject->end;
    }

    $qend = $qstart;
    $hend = $hstart;

    my $count = 0;   # counter for the bases in the alignment
    my $found = 0;   # flag saying whether we have a feature pair

    my $source = $self->waba;
    $source =~ s/\/.*\/(.*)/$1/;

    my @tmpf;

    my $analysis = new Bio::EnsEMBL::Analysis(-db              => $self->database,
                                              -db_version      => 1,
                                              -program         => $source,
                                              -program_version => 1,
                                              -gff_source      => $source,
                                              -gff_feature     => 'similarity',
                                              -logic_name      => 'blast');

    # Here goes...

    while ($count <= $#qchars) {

        # We have hit an ungapped region.  Increase the query and hit counters.
        # and flag that we have a feature pair.

        if ($qchars[$count] ne '-' &&
            $hchars[$count] ne '-') {

            $qend += $qstrand;
            $hend += $hstrand;

            $found = 1;

        } else {

            # We have hit a gapped region.  If the feature pair flag is
	    # set ($found) then make a feature pair, store it and reset
	    # the start and end variables.

            if ($found) {

	        my $fp = $self->_convert2FeaturePair($qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$hsp,$name, $analysis);
	        push @tmpf, $fp;
            }

            # We're in a gapped region.  We need to increment the sequence
            # that doesn't have the gap in it to keep the coordinates
	    # correct. We also need to reset the current end coordinates.

            if ($qchars[$count] ne '-') {
                $qstart = $qend + $qstrand;
            } else {
                $qstart = $qend;
            }
            if ($hchars[$count] ne '-') {
                $hstart = $hend + $hstrand;
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
    if ($found) {

        my $fp = $self->_convert2FeaturePair($qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$hsp,$name, $analysis);
        push @tmpf, $fp;
    }

    my $fp = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tmpf);

    # helps debugging subsequent steps
    $fp->{'qseq'} = $hsp->query_string();
    $fp->{'sseq'} = $hsp->hit_string();

    $self->growfplist($fp);
}


=head2 _convert2FeaturePair

    Arg [1]    : int $query_start
    Arg [2]    : int $query_end
    Arg [3]    : int $query_strand
    Arg [4]    : int $subject_start
    Arg [5]    : int $subject_end
    Arg [6]    : int $subject_strand
    Arg [7]    : Bio::Search::HSP::WABAHSP $hsp
    Arg [8]    : char $name
    Arg [9]    : Bio::EnsEMBL::Analysis $analysis
    Description: takes a set of coords and converts them into a FeaturePair
    Returntype : Bio::EnsEMBL::FeaturePair
    Exceptions : none
    Caller     : general

=cut

sub _convert2FeaturePair {
    my ($self, $qstart, $qend, $qstrand, $hstart, $hend, $hstrand, $hsp, $name, $analysis) = @_;

    # The actual end of the alignment is the previous character.

    $qend -= $qstrand;
    $hend -= $hstrand;

    # Make sure start is always < end

    if ($qstart > $qend) {
	($qstart, $qend) = ($qend, $qstart);
    }
    if ($hstart > $hend) {
	($hstart, $hend) = ($hend, $hstart);
    }

    return $self->_makeFeaturePair(
	$qstart, $qend, $qstrand,
	$hstart, $hend, $hstrand,
	$hsp->score,    $hsp->percent_identity,
	$hsp->pvalue,   $name, $analysis
    );
}


=head2 _makeFeaturePair

    Arg  [1]   : int $query_start
    Arg  [2]   : int $query_end
    Arg  [3]   : int $query_strand
    Arg  [4]   : int $subject_start
    Arg  [5]   : int $subject_end
    Arg  [6]   : int $subject_strand
    Arg  [7]   : int $score
    Arg  [8]   : int $percent_id
    Arg  [9]   : float $e_value
    Arg [10]   : char $hit_name
    Arg [11]   : Bio::EnsEMBL::Analysis $analysis
    Description: internal function that makes feature pairs
    Returntype : Bio::EnsEMBL::FeaturePair
    Exceptions : none
    Caller     : general

=cut

sub _makeFeaturePair {
    my ($self, $qstart, $qend, $qstrand, $hstart, $hend, $hstrand, $score, $pid, $evalue, $name, $analysis) = @_;

    my $feature1 = new Bio::EnsEMBL::SeqFeature(
	-seqname     => $self->query->id,
        -start       => $qstart,
        -end         => $qend,
        -strand      => $qstrand,
        -analysis    => $analysis,
        -score       => $score
    );
    $feature1->percent_id($pid);
    $feature1->p_value($evalue);

    my $feature2 = new Bio::EnsEMBL::SeqFeature(
	-seqname => $name,
        -start   => $hstart,
        -end     => $hend,
        -strand  => $hstrand,
        -analysis => $analysis,
        -score    => $score
    );

    my $fp = new Bio::EnsEMBL::FeaturePair(
        -feature1 => $feature1,
        -feature2 => $feature2
    );

    $feature2->percent_id($pid);
    $feature2->p_value($evalue);

    return $fp;
}

1;
