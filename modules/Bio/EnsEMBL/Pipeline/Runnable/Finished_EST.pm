
### Bio::EnsEMBL::Pipeline::Runnable::Finished_EST


package Bio::EnsEMBL::Pipeline::Runnable::Finished_EST;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_MiniEst2Genome;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast;
use Bio::EnsEMBL::Pipeline::Config::Blast;
use Data::Dumper;

use vars qw(@ISA $verbose);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

$verbose = 0;

sub new {
    my ( $class, @args ) = @_;
    
    #my $self = $class->SUPER::new(@args);
    my $self = bless {}, $class;
    my ( $query, $unmasked, $analysis, $seqfetcher, ) = $self->_rearrange(
        [
        qw{
            QUERY
            UNMASKED
            ANALYSIS
            SEQFETCHER
        }
        ],
        @args
    );
    
    die "No QUERY (masked genomic sequence) given" unless $query;
    die "No UNMASKED (genomic sequence) given"     unless $unmasked;
    
    $self->query($query);
    $self->unmasked($unmasked);
    $self->analysis($analysis);
    # this should make a OBDAIndexSeqFetcher, but will default to Pfetch
    my $indices = [ split(",", $self->analysis->db) ];
    my $seqfetcher = $self->_make_seqfetcher($indices);
    $self->seqfetcher($seqfetcher);
    
    return $self;
}

sub query {
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
            die Dumper($seq);
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_unmasked'} = $seq;
    }
    return $self->{'_unmasked'};
}

sub seqfetcher {
    my ( $self, $seqfetcher ) = @_;
    
    if ($seqfetcher) {
        $self->{'_seqfetcher'} = $seqfetcher;
    }
    return $self->{'_seqfetcher'};
}

sub analysis {
    my ( $self, $analysis ) = @_;
    
    if ($analysis) {
        $self->{'_analysis'} = $analysis;
    }
    return $self->{'_analysis'};
}

sub default_blast_parameters {
    my ( $self, $ana ) = @_;
    
    return (
        '-query'          => $self->query,
        '-database'       => $ana->db_file,
        '-program'        => $ana->program || 'wublastn',
        '-threshold_type' => 'PVALUE',
        '-threshold'      => 1e-4,
        '-options'        => 'Z=500000000 cpus=1',
    );
}

sub _make_blast_parameters {
    my ($self) = @_;
    
    # Set parameters from analysis object, or use defaults
    my $ana = $self->analysis or $self->throw('analysis not set');
    my %param = $self->default_blast_parameters($ana);
    
    my ($arguments);
    
    foreach my $ele ( split /\s*,\s*/, $ana->parameters ) {
        if ( my ( $key, $value ) = split /\s*=>\s*/, $ele ) {
            if ( defined $value ) {
                #                $param{$key} = $value;
                print STDERR "Need to use Blast Config to set module parameters\n";
                print STDERR "The analysis.parameters field [varchar(40)] should only be for program params\n";
            }
            else {
                # remaining arguments not of '=>' form
                # are simple flags (like -p1)
                $arguments .= " $key";
            }
        }
    }
    $param{'-options'} = $arguments if $arguments;
    
    return %param;
}

sub run {
    my ($self) = @_;
    
    my $seq   = $self->query;
    my $blast = Bio::EnsEMBL::Pipeline::Runnable::Finished_Blast->new( $self->_make_blast_parameters );
    # run the blast
    $blast->run();
    # keep the features
    my $features = [ $blast->output ];
    # set the db_version_searched
    my $dbv = $blast->db_version_searched;
    print "using $dbv\n";
    $self->db_version_searched($dbv);
    
    # print STDERR Dumper($features);
    
    print STDERR "\nPlus strand est_genome\n" if $verbose;
    $self->run_est_genome_on_strand( 1, $features );
    
    print STDERR "\nMinus strand est_genome\n" if $verbose;
    $self->run_est_genome_on_strand( -1, $features );
    
}

sub run_est_genome_on_strand {
    my ( $self, $strand, $feat ) = @_;
    
    my $hit_features = {};
    for ( my $i = 0 ; $i < @$feat ; $i++ ) {
        my $f   = $feat->[$i];
        my $hid = $f->hseqname or $self->throw("Missing hid");
        
        next unless $f->strand == $strand;
        $hit_features->{$hid} ||= [];
        push ( @{ $hit_features->{$hid} }, $f );
    }
    
    my ($is_linear);
    if ( $strand == 1 ) {
        $is_linear = sub {
            my ( $x, $y ) = @_;
            
            return $x->hstart <= $y->hstart;
        };
    }
    else {
        $is_linear = sub {
            my ( $x, $y ) = @_;
            return $x->hend >= $y->hend;
        };
    }
    
    while ( my ( $hid, $flist ) = each %$hit_features ) {
        
        print STDERR "$hid\n" if $verbose;
        print STDERR "PRESORTED\n" if $verbose; #
        foreach my $f (@$flist) {
            print STDERR "hit:  " . $f->gffstring . "\n" if $verbose;
        }
        
        # Sort features by start/end in genomic
        @$flist = sort { $a->start <=> $b->start or $a->end <=> $b->end } @$flist;
        
        print STDERR "SORTED\n" if $verbose;
        foreach my $f (@$flist) {
            print STDERR "hit:  " . $f->gffstring . "\n" if $verbose;
        }
        
        # Group into linear matches
        my @sets = ( [ $flist->[0] ] );
        my $curr = $sets[0];
        for ( my $i  = 1 ; $i < @$flist ; $i++ ) {
            
            my $prev = $flist->[ $i - 1 ];
            my $this = $flist->[$i];
            
            if ( &$is_linear( $prev, $this ) ) {
                push ( @$curr, $this );
            }
            else {
                $curr = [$this];
                push ( @sets, $curr );
            }
        }
        
        print STDERR "MADE ",scalar(@sets)," LINEAR MATCHES SET(S)\n" if $verbose;
        foreach my $lin (@sets) {
            $self->do_mini_est_genome($lin);
        }
        #        $self->do_mini_est_genome($flist);
    }
}

sub do_mini_est_genome {
    my ( $self, $linear ) = @_;
    
    print STDERR "\nLinear Matches\n" if $verbose;
    foreach my $lin (@$linear) {
        print STDERR "before:  " . $lin->gffstring . "\n" if $verbose;
    }
    ### Is this merging features?  - May be cause of duplicate features if isn't
    my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::Finished_MiniEst2Genome(
        '-genomic'    => $self->unmasked,
        '-features'   => $linear,
        '-seqfetcher' => $self->seqfetcher,
    );
    $e2g->run;
    
    
    foreach my $fp ( $e2g->output ) {
        print STDERR "after:  " . $fp->gffstring . "\n" if $verbose;
        # source_tag and primary_tag have to be set to
        # something, or validate method in FeaturePair
        # (callled by RunnableDB) throws!
        #        $fp->feature2->source_tag('I_am_valid');
        #        $fp->feature2->primary_tag('I_am_valid');
        $self->add_output($fp);
    }
}


sub add_output {
    my ( $self, @feat ) = @_;
    
    my $ana = $self->analysis;
    foreach my $f (@feat) {
        $f->analysis($ana);
    }
    
    $self->{'_output'} ||= [];
    
    push ( @{ $self->{'_output'} }, @feat );
    
}

sub output {
    my ($self) = @_;
    
    if ( my $out = $self->{'_output'} ) {
        return @$out;
    }
    else {
        return;
    }
}
sub _make_seqfetcher {
    my ( $self, $indices ) = @_;
    
    my $seqfetcher;
    if ( ref($indices) eq "ARRAY"){
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(-db => $indices);
    }
    else{
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();
    }
    return $seqfetcher;
}
sub db_version_searched{
    my $self = shift;
    
    $self->{'_db_version_searched'} = shift if @_;
    print "I've decided to search db version: " . $self->{'_db_version_searched'} . "\n";
    return $self->{'_db_version_searched'};
}
1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

