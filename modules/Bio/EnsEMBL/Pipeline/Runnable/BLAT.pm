

# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BLAT

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Blat->new(
								 -database    => $database,
								 -query  => \@sequences,
								 -blat        => $self->blat,
								 -options     => $self->options,
								);

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
 where @results is an array of SeqFeatures, each one representing an aligment (e.g. a transcript), 
 and each feature contains a list of alignment blocks (e.g. exons) as sub_SeqFeatures, which are
 in fact feature pairs.
 
=head1 DESCRIPTION

Blat takes a Bio::Seq (or Bio::PrimarySeq) object and runs Blat
against a set of sequences.  The resulting output file is parsed
to produce a set of features.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BLAT;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use Carp 'cluck';

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'} = [];
    my ( $database, $query, $blat, $options ) = $self->_rearrange(
        [
            qw(
            DATABASE
            QUERY
            BLAT
            OPTIONS
            )
        ],
        @args
    );

    if ($query) {
        $self->clone($query);
    }
    else {
        $self->throw("No query sequence input.");
    }

    # you can pass a sequence object for the target or a database (multiple fasta file);
    if ($database) {
        $self->database($database);
    }
    else {
        $self->throw("Blat needs a target - database: $database");
    }

    # can choose which blat to use
    $self->blat('blat') unless $blat;
    $self->blat( $self->find_executable($blat) );

    # can add extra options as a string
    if ($options) {
        $self->options($options);
    }
    return $self;
}

sub clone {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }

        $self->{'_query'} = $seq;

        $self->filename( $self->clone->id . ".$$.seq" );
        $self->results( $self->filename . ".blat.out" );

    }
    return $self->{'_query'};
}

=head2 run

Usage   :   $obj->run($workdir, $args)
Function:   Runs blat script and puts the results into the file $self->results
            It calls $self->parse_restuls, and results are stored in $self->output
=cut
sub run {
    my ($self) = @_;

    # set the working directory (usually /tmp)
    print STDERR" working directory " . $self->workdir() . "\n";
    $self->workdir('/tmp') unless ( $self->workdir() );
    $self->checkdir();
    $self->writefile();
    my $gen_query = join ( '/', $self->workdir, $self->filename );
    my $options = $self->options;
    my $blat    = $self->blat;

    my @db_list = $self->fetch_databases;
    #my $command =
#      "sed 's/>\\([:A-Z0-9:]*\\) \\([:A-Z0-9.:]*\\)/>\\2 \\1/' $db_list[0] | $blat stdin $gen_query $options"
#      . $self->results;

    foreach my $db (@db_list) {
        my $command     = "/acari/work7a/keenan/tools/reformat_headers $db | $blat stdin $gen_query $options " . $self->results;
        print $command, "\n";
        $self->throw("$command\nFailed during BLAT run $!") unless ( system($command) == 0 );
        $self->parse_results();
    }

    #close BLAT or die "Error running blat pipe '$cmd' : exit($?)";
    #$self->deletefiles;
}


############################################################

=head2 parse_results

 Usage   :   $obj->parse_results
 Function:   reads the Blat output (in PSL format or psLayout ) which has been written to
             a local file $self->results. can accept filenames, filehandles or pipes (\*STDIN)
 Returns :   a list of Features (each alignment), each feature with a list of sub_Seqfeatures (each block of the alignment)
 Args    :   optional filename

=cut

sub parse_results {
    my ( $self, $filehandle ) = @_;

    my $filehandle;
    my $resfile = $self->results();

    if ( -e $resfile ) {
        if ( -z $self->results ) {
            print STDERR "Blat didn't find any matches\n";
            return;
        }
        else {
            open( OUT, "<$resfile" ) or $self->throw( "Error opening ", $resfile, " \n" );
            $filehandle = \*OUT;
        }
    }
    else {    #it'a a filehandle
        $filehandle = $resfile;
    }

    my @feats;

    #extract values
    while (<OUT>) {

        ############################################################
        #  PSL lines represent alignments and are typically taken from files generated 
        # by BLAT or psLayout. See the BLAT documentation for more details. 
        #
        # 1.matches - Number of bases that match that aren't repeats 
        # 2.misMatches - Number of bases that don't match 
        # 3.repMatches - Number of bases that match but are part of repeats 
        # 4.nCount - Number of 'N' bases 
        # 5.qNumInsert - Number of inserts in query 
        # 6.qBaseInsert - Number of bases inserted in query 
        # 7.tNumInsert - Number of inserts in target 
        # 8.tBaseInsert - Number of bases inserted in target 
        # 9.strand - '+' or '-' for query strand. In mouse, second '+'or '-' is for genomic strand 
        #10.qName - Query sequence name 
        #11.qSize - Query sequence size 
        #12.qStart - Alignment start position in query 
        #13.qEnd - Alignment end position in query 
        #14.tName - Target sequence name 
        #15.tSize - Target sequence size 
        #16.tStart - Alignment start position in target 
        #17.tEnd - Alignment end position in target 
        #18.blockCount - Number of blocks in the alignment 
        #19.blockSizes - Comma-separated list of sizes of each block 
        #20.qStarts - Comma-separated list of starting positions of each block in query 
        #21.tStarts - Comma-separated list of starting positions of each block in target 
        ############################################################

        # first split on spaces:
        chomp;
        
        print $_;
       
        my (
            $matches,      $mismatches,    $rep_matches, $n_count, $q_num_insert, $q_base_insert,
            $t_num_insert, $t_base_insert, $strand,      $q_name,  $q_length,     $q_start,
            $q_end,        $t_name,        $t_length,    $t_start, $t_end,        $block_count,
            $block_sizes,  $q_starts,      $t_starts
          )
          = split;

        my $superfeature = Bio::EnsEMBL::SeqFeature->new();

        # ignore any preceeding text
        unless ( $matches =~ /^\d+$/ ) {
            next;
        }

        #print $_."n";

        # create as many features as blocks there are in each output line
        my ( %feat1, %feat2 );

        # all the block sizes add up to $matches + $mismatches + $rep_matches
        # percentage identity =  ( matches not in repeats + matches in repeats ) / ( alignment length )
        #print STDERR "calculating percent_id and score:\n";
        #print STDERR "matches: $matches, rep_matches: $rep_matches, mismatches: $mismatches, q_length: $q_length\n";
        #print STDERR "percent_id = 100x".($matches + $rep_matches)."/".( $matches + $mismatches + $rep_matches )."\n";
        my $percent_id = sprintf "%.2f",
          ( 100 * ( $matches + $rep_matches ) / ( $matches + $mismatches + $rep_matches ) );

        # or is it ...?
        ## percentage identity =  ( matches not in repeats + matches in repeats ) / query length
        #my $percent_id = sprintf "%.2d", (100 * ($matches + $rep_matches)/$q_length );

        # we put basically score = coverage = ( $matches + $mismatches + $rep_matches ) / $q_length
        #print STDERR "score = 100x".($matches + $mismatches + $rep_matches)."/".( $q_length )."\n";

        my $score = sprintf "%.2f", ( 100 * ( $matches + $mismatches + $rep_matches ) / $q_length );

        # size of each block of alignment (inclusive)
        my @block_sizes = split ",", $block_sizes;

        # start position of each block (you must add 1 as psl output is off by one in the start coordinate)
        my @q_start_positions = split ",", $q_starts;
        my @t_start_positions = split ",", $t_starts;

        $superfeature->seqname($q_name);
        $superfeature->score($score);
        $superfeature->percent_id($percent_id);

        # each line of output represents one possible entire aligment of the query (feat1) and the target(feat2)
        for ( my $i = 0 ; $i < $block_count ; $i++ ) {

            my ( $query_start, $query_end );

            # qStart = qSize - revQEnd
            # qEnd = qSize - revQStart

            if ( $strand eq '+' ) {
                $query_start = $q_start_positions[$i] + 1;
                $query_end   = $query_start + $block_sizes[$i] - 1;
            }
            else {
                $query_end   = $q_length - $q_start_positions[$i];
                $query_start = $query_end - $block_sizes[$i] +1;
            }

            $feat2{db}         = undef;
            $feat2{db_version} = undef;
            $feat2{program}    = 'blat';
            $feat2{p_version}  = '1';
            $feat2{logic_name} = 'Blat';

            $feat1{name}    = $q_name;
            $feat1{start}   = $query_start;
            $feat1{end}     = $query_end;
            $feat1{strand}  = $strand;
            $feat1{score}   = $score;
            $feat1{source}  = 'blat';
            $feat1{primary} = 'similarity';
            $feat1{percent} = $percent_id;
            $feat1{p}       = undef;

            $feat2{name}    = $t_name;
            $feat2{start}   = $t_start_positions[$i] + 1;
            $feat2{end}     = $feat2{start} + $block_sizes[$i] - 1;
            $feat2{strand}  = 1;
            $feat2{score}   = $score;
            $feat2{source}  = 'blat';
            $feat2{primary} = 'similarity';
            $feat2{percent} = $percent_id;
            $feat2{p}       = undef;

            my $fp = $self->createfeaturepair( \%feat1, \%feat2 );

        }

    }

    close $filehandle;

}

sub blat {
    my ( $self, $location ) = @_;
    if ($location) {
        $self->throw("Blat not found at $location: $!\n") unless ( -e $location );
        $self->{_blat} = $location;
    }
    return $self->{_blat};
}

sub options {
    my ( $self, $options ) = @_;
    if ($options) {
        $self->{_options} = $options;
    }
    return $self->{_options};
}

sub output {
    my ( $self, @arg ) = @_;

    if (@arg) {
        @{ $self->{'_fplist'} } = @arg;
    }

    if ( !defined( $self->{'_fplist'} ) ) {
        $self->{'_fplist'} = [];
    }

    return @{ $self->{'_fplist'} };
}

sub database {
    my ( $self, $database ) = @_;
    if ($database) {
        $self->{_database} = $database;
    }
    return $self->{_database};
}

sub fetch_databases {
    my ($self) = @_;

    my $db_string = $self->database;
    my @split_db_string = split /,/, $db_string;
    warn @split_db_string, "\n\n";
    my @databases;

    foreach my $db (@split_db_string) {

        my $fulldbname;

        # If we have passed a full path name don't append the $BLASTDB
        # environment variable.

        if ( $db =~ /\// ) {
            $fulldbname = $db;
        }
        else {
            $fulldbname = $ENV{BLASTDB} . "/" . $db;
        }

        # If the expanded database name exists put this in
        # the database array.
        #
        # If it doesn't exist then see if $database-1,$database-2 exist
        # and put them in the database array

        if ( !( -d $fulldbname ) && -e $fulldbname ) {
            push ( @databases, $db );
        }
        else {
            my $count = 1;

            while ( -e $fulldbname . "-$count" ) {
                push ( @databases, $fulldbname . "-$count" );
                $count++;
            }
        }

    }

    if ( scalar(@databases) == 0 ) {
        $self->throw( "No databases exist for " . $self->database );
    }

    return @databases;

}

1;
