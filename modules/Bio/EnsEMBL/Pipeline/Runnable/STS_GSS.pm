package Bio::EnsEMBL::Pipeline::Runnable::STS_GSS;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;

use FileHandle;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

  Arg [1]   : hash of parameters 
  Function  : to create a new STS_GSS object
  Returntype: an STS_GSS object
  Exceptions: It will throw is it recieves no unmasked sequence, it will throw if it isn''t running a blast an no blast features are passed to it. If it is running a blast it will throw if it receives no blast db location or masked sequence. It will also throw if it gets no seqfetcher 
  Caller    : 
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);

    $self->{'_query'}           = undef;    # location of Bio::Seq object
    $self->{'_unmasked'}        = undef;
    $self->{'_blast_program'}   = undef;    # location of Blast
    $self->{'_est_program'}     = undef;    # location of Blast
    $self->{'_database'}        = undef;    # name of database
    $self->{'_threshold'}       = undef;    # Threshold for hit filterting
    $self->{'_options'}         = undef;    # arguments for blast
    $self->{'_filter'}          = 0;        # Do we filter features?
    $self->{'_fplist'}          = [];       # an array of feature pairs (the output)
    $self->{'_no_blast'}        = 0;
    $self->{'_features'}        = [];
    $self->{'_workdir'}         = undef;    # location of temp directory
    $self->{'_filename'}        = undef;    # file to store Bio::Seq object
    $self->{'_results'}         = undef;    # file to store results of analysis
    $self->{'_seqfetcher'}      = undef;
    $self->{'_prune'}           = 0;        #
    $self->{'_coverage'}        = 10;
    $self->{'_merged_features'} = [];
    $self->{'_percent_id'}      = undef;
    $self->{'_padding'}         = undef;
    $self->{'_percent_filter'}  = undef;
    $self->{'_tandem'}          = undef;

    # Now parse the input options and store them in the object
    #print "@args\n";
    my (
        $query,    $unmasked,   $program, $database, $threshold,  $threshold_type,
        $filter,   $coverage,   $prune,   $options,  $seqfetcher, $features,
        $no_blast, $percent_id, $padding, $tandem
      )
      = $self->_rearrange(
        [
            qw(QUERY
            UNMASKED
            PROGRAM
            DATABASE
            THRESHOLD
            THRESHOLD_TYPE
            FILTER
            COVERAGE
            PRUNE
            OPTIONS
            SEQFETCHER
            FEATURES
            NO_BLAST
            PERCENT_ID
            PADDING
            PERCENT_FILTER
            TANDEM_CHECK)
        ],
        @args
    );

    if ( $no_blast == 1 ) {

        #print "setting no_blast = ".$no_blast."\n";
        $self->{'_no_blast'} = $no_blast;
    }

    #print "1no_blast = ".$no_blast."\n";

    #print "2no_blast = ".$self->no_blast."\n";
    #print "running new\n";
    #print "1unmasked = ".$unmasked."\n";
    if ($unmasked) {
        $self->unmasked($unmasked);
    }
    else {
        $self->throw("No unmasked query sequence input.");
    }

    #print "2unmasked = ".$self->unmasked."\n";
    #print STDERR "find executable output ".$self->find_executable($program)."\n";
    if ($program) {
        $self->program( $self->find_executable($program) );
    }
    if ( $self->no_blast ) {
        if ( defined($features) ) {
            if ( ref($features) eq "ARRAY" ) {
                push ( @{ $self->{'_features'} }, @$features );
            }
            else {
                $self->throw("[$features] is not an array ref.");
            }
        }
        else {
            $self->throw("should pass in features haven't \n");
        }
    }
    if ( $no_blast == 0 ) {

        #print $no_blast."\n";
        if ($query) {
            $self->clone($query);
        }
        else {
            if ( !$no_blast ) {
                $self->throw("No query sequence input.");
            }
        }
        if ($database) {
            $self->database($database);
        }
        else {
            {
                $self->throw("No database input");
            }
        }
    }
    if ($options) {

        #this option varies the number of HSP displayed proportionally to the query contig length
        $self->options($options);
    }
    else {
        $self->options(' -cpus=1 ');
    }
    if ($percent_id) {

        $self->percent_id($percent_id);
    }
    else {
        $self->percent_id(0);
    }
    if ($percent_id) {

        $self->padding($padding);
    }
    else {
        $self->padding(400);
    }

    #print "options = ".$self->options."\n";
    if ( defined($threshold) ) {
        $self->threshold($threshold);
    }
    else {
        $self->threshold(0);
    }

    if ( defined($threshold_type) ) {
        $self->threshold_type($threshold_type);
    }
    else {
        $self->threshold_type('PID');
    }

    if ( defined($filter) ) {
        $self->filter($filter);
    }

    if ( defined($prune) ) {
        $self->prune($prune);
    }
    if ( defined($coverage) ) {
        $self->coverage($coverage);
    }
    if ( defined($tandem) ) {
        $self->tandem_check($tandem);
    }
    $self->throw("No seqfetcher provided") unless defined($seqfetcher);

    $self->seqfetcher($seqfetcher) if defined($seqfetcher);
    return $self;    # success - we hope!
}

#################
#get/set methods#
#################

=head2 accessor methods

  Arg [1]   : variable to be set 
  Function  : if passed set the varible to the argument passed the return the argument
  Returntype: the varible to be set
  Exceptions: some throw exceptions if not passed the correct varible
  Caller    : $self
  Example   : $clone = $self->clone;

=cut

sub clone {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }

        $self->{'_query'} = $seq;

    }
    return $self->{'_query'};
}

sub unmasked {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }

        $self->{'_unmasked'} = $seq;
        $self->filename( $self->unmasked->id . ".$$.seq" );
        $self->results( $self->filename . ".blast.out" );

    }
    return $self->{'_unmasked'};
}

sub program {
    my ( $self, $location ) = @_;

    if ($location) {
        $self->throw("executable not found at $location: $!\n") unless ( -e $location && -x $location );
        $self->{'_blast_program'} = $location;
    }
    return $self->{'_blast_program'};
}

sub database {
    my ( $self, $db ) = @_;

    if ( defined($db) ) {
        $self->{'_database'} = $db;
    }
    return $self->{'_database'};
}

sub no_blast {
    my ($self) = @_;

    return $self->{'_no_blast'};
}

sub options {
    my ( $self, $args ) = @_;

    if ( defined($args) ) {
        $self->{'_options'} = $args;
    }
    return $self->{'_options'};
}

sub percent_filter {
    my ( $self, $args ) = @_;

    if ( defined($args) ) {
        $self->{'_percent_filter'} = $args;
    }
    return $self->{'_percent_filter'};
}

sub percent_id {
    my ( $self, $args ) = @_;

    if ( defined($args) ) {
        $self->{'_percent_id'} = $args;
    }
    return $self->{'_percent_id'};
}

sub padding {
    my ( $self, $args ) = @_;

    if ( defined($args) ) {
        $self->{'_padding'} = $args;
    }
    return $self->{'_padding'};
}

sub prune {
    my ( $self, $arg ) = @_;

    if ( defined($arg) ) {
        $self->{_prune} = $arg;
    }
    return $self->{_prune};
}

sub tandem_check {
    my ( $self, $arg ) = @_;

    if ( defined($arg) ) {
        $self->{_tandem} = $arg;
    }
    return $self->{_tandem};
}

sub coverage {
    my ( $self, $arg ) = @_;

    if ( defined($arg) ) {
        $self->{_coverage} = $arg;
    }
    return $self->{_coverage};
}

sub filter {
    my ( $self, $args ) = @_;

    if ( defined($args) ) {
        if ( $args != 0 && $args != 1 ) {
            $self->throw("Filter option must be 0 or 1");
        }
        $self->{'_filter'} = $args;
    }
    return $self->{'_filter'};
}

sub get_threshold_types {
    my ($self) = @_;

    return ( "PID", "PVALUE" );
}

sub threshold_type {
    my ( $self, $type ) = @_;

    my @types = $self->get_threshold_types;

    if ( defined($type) ) {
        my $found = 0;
        foreach my $allowed_type ( $self->get_threshold_types ) {
            if ( $type eq $allowed_type ) {
                $found = 1;
            }
        }
        if ( $found == 0 ) {

            $self->throw("Type [$type] is not an allowed type.  Allowed types are [@types]");
        }
        else {
            $self->{_threshold_type} = $type;
        }
    }
    return $self->{_threshold_type} || $types[0];
}

sub get_pars {
    my ($self) = @_;

    if ( !defined( $self->{_hits} ) ) {
        $self->{_hits} = [];
    }

    return @{ $self->{_hits} };

}

sub seqfetcher {
    my ( $self, $value ) = @_;
    if ($value) {

        #need to check if passed sequence is Bio::DB::RandomAccessI object
        $value->isa("Bio::DB::RandomAccessI") || $self->throw("Input isn't a Bio::DB::RandomAccessI");
        $self->{'_seqfetcher'} = $value;
    }
    return $self->{'_seqfetcher'};
}

sub add_merged_feature {
    my ( $self, $feature ) = @_;
    if ( !$feature ) {
        $self->throw("no feature passed : $!\n");
    }
    else {
        push ( @{ $self->{'_merged_features'} }, $feature );
    }
}

sub each_merged_feature {

    my ($self) = @_;

    return @{ $self->{'_merged_features'} };

}

sub features {
    my ($self) = @_;
    return @{ $self->{'_features'} };
}

###########
# Analysis methods
##########

=head2 run

  Arg [1]   : directory result to be written too (optional) 
  Function  : runs the blast and est2genome analysis on the given contig
  Returntype: none
  Exceptions: if a blast analysis is to be run it throws if no masked sequence is present
  Caller    : 
  Example   : 

=cut

sub run {
    my ( $self, $dir ) = @_;

    $self->workdir('/tmp') unless ( $self->workdir($dir) );
    $self->checkdir();

    #print STDERR "checked directories\n";
    #write sequence to file
    #$self->writefile('unmasked');
    my @blast_output;
    if ( $self->no_blast == 1 ) {

        @blast_output = $self->features;

    }
    else {
        my $seq = $self->clone || $self->throw("Query seq required for blast\n");

        @blast_output = $self->run_blasts();
    }

    #print STDERR "ran analysis\n";
    #parse output and create features
    print "there are " . scalar(@blast_output) . " features\n";
    $self->expand_and_merge_features( \@blast_output );
    my @features = $self->each_merged_feature;
    my @results;
    foreach my $feature (@features) {
        print STDERR "running on " . $feature->hseqname . " with score " . $feature->score . " and evalue "
          . $feature->p_value . "\n";

        #print "there are ".scalar(@genes)." results\n";
        push ( @results, $self->run_est2genome($feature) );
    }
    $self->deletefiles();

    $self->output( \@results );

}

=head2 run_blasts

  Arg [1]   : none 
  Function  : run blast on given contig and return results 
  Returntype: array of blast results
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub run_blasts {

    my ($self) = @_;

    #print "running blasts\n";

    my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(
        -query          => $self->clone,
        -program        => $self->program,
        -database       => $self->database,
        -threshold      => $self->threshold,
        -threshold_type => $self->threshold_type,
        -filter         => $self->filter,
        -coverage       => $self->coverage,
        -prune          => $self->prune,
        -options        => $self->options,

    );

    $blast->run;

    my @features = $blast->output;

    return @features;

}

=head2 expand_and_merge_features

  Arg [1]   : ref to array of feature pairs 
  Function  : filter features on percent_id, extend features start and end to cover whole hit plus padding, and merge overlapping features
  Returntype: array of merged features
  Exceptions: throws if features don't have 1 or -1 as a strand'
  Caller    : 
  Example   : 

=cut

sub expand_and_merge_features {
    my ( $self, $features ) = @_;

    my @features;
    my @unfiltered = @$features;

    if ( $self->percent_filter ) {

        #print "filtering on percent_id\n";
        foreach my $unf (@unfiltered) {

            #print $unf->percent_id."\n"; 
            if ( $unf->percent_id >= $self->percent_id ) {
                push ( @features, $unf );
            }
        }
    }
    else {
        push ( @features, @unfiltered );
    }

    #return @features;
    #print "there are ".scalar(@features)." unmerged features\n";
    my %unique_hids;

    foreach my $feature (@features) {
        my $hseqname = $feature->hseqname();

        if ( !$unique_hids{$hseqname} ) {
            $unique_hids{$hseqname} = [];
        }
        push ( @{ $unique_hids{$hseqname} }, $feature );
    }

    my @hids = keys(%unique_hids);

    ## For each sequence hit, all the features which correspond to that sequence
    ## have their start and end extended.  The start and end are extended so
    ## there is enough query sequence to cover the entire hit sequence plus a bit of padding
    ## which is set to 200 as default    
    my $count          = 0;
    my $padding        = $self->padding;
    my $genomic_length = $self->unmasked->length;
    foreach my $hid (@hids) {

        #print "feature ".$hid." has ".scalar(@{$unique_hids{$hid}})." hits\n";
        my $feature_array = $unique_hids{$hid};

        #print $self->seqfetcher."\n";
        my $hid_seq = $self->seqfetcher->get_Seq_by_acc($hid);

        my $hid_len = $hid_seq->length;

        foreach my $feature (@$feature_array) {

#print "before seqname = ".$hid." length = ".$hid_len." start ".$feature->start." end ".$feature->end." strand ".$feature->strand." hstart ".$feature->hstart." hend ".$feature->hend."\n";
            my $genomic_start = $feature->start;
            my $genomic_end   = $feature->end;
            my $hstart        = $feature->hstart;
            my $hend          = $feature->hend;
            if ( $feature->strand == 1 ) {
                $feature->start( $genomic_start - ( $hstart + $padding ) );
                if ( $feature->start < 1 ) {
                    $feature->start(1);
                }
                $feature->end( $genomic_end + ( $hid_len - $hend ) + $padding );
                if ( $feature->end > $genomic_length ) {
                    $feature->end($genomic_length);
                }
            }
            elsif ( $feature->strand == -1 ) {
                $feature->start( $genomic_start - ( $hid_len - $hend ) - $padding );
                if ( $feature->start < 1 ) {
                    $feature->start(1);
                }
                $feature->end( $genomic_end + $hstart + $padding );

                #print $count." trying to get unmasked from ".$self."\n";
                $count++;
                if ( $feature->end > $genomic_length ) {
                    $feature->end($genomic_length);
                }
            }
            else {
                $self->throw( $feature->hseqname . " has got " . $feature->strand . " as a stand : $!\n" );
            }

#print "after seqname = ".$hid."start ".$feature->start." end ".$feature->end." hstart ".$feature->hstart." hend ".$feature->hend." strand ".$feature->strand."\n";
        }

## the features are then sorted into array of features on the forward and reverse strand so overlapping features can be merged##
        my @sorted_forward_features;
        my @sorted_reverse_features;
        foreach my $feature (@$feature_array) {
            my $strand = $feature->strand;
            if ( $strand == -1 ) {
                push ( @sorted_reverse_features, $feature );
            }
            elsif ( $strand == 1 ) {
                push ( @sorted_forward_features, $feature );
            }
            else {
                $self->throw( $feature->hid . " has strand '$strand' : $!\n" );
            }
        }
        @sorted_forward_features = sort { $a->start <=> $b->start || $a->end <=> $b->end } @sorted_forward_features;
        @sorted_reverse_features = sort { $a->start <=> $b->start || $a->end <=> $b->end } @sorted_reverse_features;
        $self->fuse_overlapping_features( \@sorted_forward_features );
        $self->fuse_overlapping_features( \@sorted_reverse_features );
    }

    #my @merged = $self->each_merged_feature;

    #print "there are ".scalar(@merged)." merged features\n";
}

=head2 fuse_overlapping_features

  Arg [1]   : refence to an array of features 
  Function  : merge together overapping features abutting features then adds each merged feature to an array
  Returntype: none
  Exceptions: throws if feature start becomes greater than end
  Caller    : 
  Example   : 

=cut

sub fuse_overlapping_features {
    my ( $self, $features ) = @_;

    for ( my $i = 1 ; $i < @$features ; ) {
        my $f_a = $features->[ $i - 1 ];
        my $f_b = $features->[$i];

        # If coordinates overlap or abut
        if ( $f_a->end + 1 >= $f_b->start ) {

            # Transfer the end of b to a if b ends beyond a
            if ( $f_a->end < $f_b->end ) {
                $f_a->end( $f_b->end );
            }

            # Remove b from the list
            splice( @$features, $i, 1 );

            # Critical error check!
            if ( $f_a->start > $f_a->end ) {
                die "Feature start (", $f_a->start, ") greater than end (", $f_a->end, ")";
            }
        }
        else {

            # Haven't spliced, so move pointer
            $i++;
        }
    }

    foreach my $feat (@$features) {
        $self->add_merged_feature($feat);
    }
}

=head2 run_est2genome

  Arg [1]   :a feature 
  Function  : this uses the finshed_est2genome module to run est2genome on a given subseq of the contig and a particular est/gss/sts seq 
  Returntype: returns and array of supporting features for that est/gss hit
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub run_est2genome {

    my ( $self, $feature ) = @_;

    my $est    = $self->seqfetcher->get_Seq_by_acc( $feature->hseqname );
    my $seq    = $self->unmasked->subseq( $feature->start, $feature->end );
    my $start  = $feature->start;
    my $strand = $feature->strand;
    my $end    = $feature->end;

    my $query_seq = Bio::Seq->new(
        -seq       => $seq,
        -id        => $self->unmasked->id,
        -accession => $self->unmasked->id
    );

    my $est_length = length( $est->primary_seq->seq );

    #print "running est2genome on ".$feature->hseqname." strand ".$feature->strand." start = ".$start." end ".$end."\n";
    my $est2genome = Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome->new(
        -genomic => $query_seq,
        -est     => $est,
    );

    $est2genome->run;
    my @features = $est2genome->output;

    #print "est2 genome outputted ".scalar(@features)."\n";
    my @output;
    foreach my $f (@features) {

        my $coverage = ( $f->length / $est_length );
        $f->start( $f->start + $start - 1 );
        $f->end( $f->end + $start - 1 );

        if ( ( $f->hstart <= 5 || $f->hend >= ( $est_length - 4 ) ) && $coverage >= 0.9 && $f->percent_id >= 95 ) {
            push ( @output, $f );
        }
    }

    #print "there are ".scalar(@output)." outputted\n";

    return @output;

}

=head2 output

  Arg [1]   : refence to array to be set for output (optional)
  Function  : pushs any given array ref on to the outout array 
  Returntype: returns the output array
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub output {

    my ( $self, $output ) = @_;
    if ($output) {
        my @outputs = @$output;
        push ( @{ $self->{'_fplist'} }, @outputs );
    }

    return @{ $self->{'_fplist'} };
}

1;
