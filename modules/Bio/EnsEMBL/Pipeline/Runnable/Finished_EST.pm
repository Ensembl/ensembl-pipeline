
### Bio::EnsEMBL::Pipeline::Runnable::Finished_EST

package Bio::EnsEMBL::Pipeline::Runnable::Finished_EST;

use strict;
use Data::Dumper;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;

@ISA = ('Bio::EnsEMBL::Pipeline::RunnableI');

sub new {
    my ($class,@args) = @_;

    #my $self = $class->SUPER::new(@args);
    my $self = bless {}, $class;
    my ($query,
        $unmasked,
        $analysis,
        $seqfetcher,
        ) = $self->_rearrange([qw{
            QUERY
            UNMASKED
            ANALYSIS
            SEQFETCHER
            }], @args);

        #$database,
        #$program,
        #$options,
        #$threshold,
        #$thres_type,
        #    DATABASE
        #    PROGRAM
        #    OPTIONS
        #    THRESHOLD
        #    THRESHOLD_TYPE
    
    die "No QUERY (masked genomic sequence) given" unless $query;
    die "No UNMASKED (genomic sequence) given"     unless $unmasked;

    $self->query            ($query);
    $self->unmasked         ($unmasked);
    $self->analysis         ($analysis);
    $self->seqfetcher       ($seqfetcher);
    #$self->database         ($database);
    #$self->program          ($program);
    #$self->options          ($options);
    #$self->threshold        ($threshold);
    #$self->threshold_type   ($thres_type);
    
    return $self;
}

sub query {
    my ($self, $seq) = @_;

    if ($seq) {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_query'} = $seq ;
    }
    return $self->{'_query'};
}

sub unmasked {
    my ($self, $seq) = @_;

    if ($seq) {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_unmasked'} = $seq ;
    }
    return $self->{'_unmasked'};
}

sub seqfetcher {
    my( $self, $seqfetcher ) = @_;
    
    if ($seqfetcher) {
        $self->{'_seqfetcher'} = $seqfetcher;
    }
    return $self->{'_seqfetcher'};
}

sub analysis {
    my( $self, $analysis ) = @_;
    
    if ($analysis) {
        $self->{'_analysis'} = $analysis;
    }
    return $self->{'_analysis'};
}

#sub database {
#    my( $self, $database ) = @_;
#    
#    if ($database) {
#        $self->{'_database'} = $database;
#    }
#    return $self->{'_database'} || 'dbEST';
#}
#
#sub program {
#    my( $self, $program ) = @_;
#    
#    if ($program) {
#        $self->{'_program'} = $program;
#    }
#    return $self->{'_program'} || 'wublastn';
#}
#
#sub options {
#    my( $self, $options ) = @_;
#    
#    if ($options) {
#        $self->{'_options'} = $options;
#    }
#    return $self->{'_options'} || 'Z=500000000';
#}
#
#sub threshold {
#    my( $self, $threshold ) = @_;
#    
#    if ($threshold) {
#        $self->{'_threshold'} = $threshold;
#    }
#    return $self->{'_threshold'} || 1e-4;
#}
#
#sub threshold_type {
#    my( $self, $threshold_type ) = @_;
#    
#    if ($threshold_type) {
#        $self->{'_threshold_type'} = $threshold_type;
#    }
#    return $self->{'_threshold_type'} || 'PVALUE';
#}

#sub _make_blast_paramters {
#    my( $self ) = @_;
#    
#    my @param;
#    foreach my $meth (qw{
#        query
#        database
#        program
#        options
#        threshold
#        threshold_type
#        })
#    {
#        push(@param, "-$meth", $self->$meth());
#    }
#    return @param;
#}

sub _make_blast_paramters {
    my( $self ) = @_;
    
    # Set parameters from analysis object, or use defaults
    my $ana = $self->analysis or $self->throw('analysis not set');
    
    my %param = (
        '-query'            => $self->query,
        '-database'         => $ana->db      || 'dbEST',
        '-program'          => $ana->program || 'wublastn',
        '-threshold_type'   => 'PVALUE',
        '-threshold'        => 1e-4,
        '-options'          => 'Z=500000000',
        );
    
    my( $arguments );
    foreach my $ele (split /\s*,\s*/, $ana->parameters) {
        if (my ($key, $value) = split (/\s*=>\s*/, $ele);
            if (defined $value) {
                if ($key eq '-threshold_type' or $key eq '-threshold') {
                    $param{$key} = $value;
                }
            } else {
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
    my( $self ) = @_;
    
    my $seq = $self->query;
    my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast
        ->new($self->_make_blast_paramters);
    $blast->run;
    my $features = [$blast->output];
    
    $self->run_est_genome_on_strand( 1, $features);
    $self->run_est_genome_on_strand(-1, $features);
}

sub run_est_genome_on_strand {
    my( $self, $strand, $feat ) = @_;
    
    my $hit_features = {};
    for (my $i = 0; $i < @$feat; $i++) {
        my $f   = $feat->[$i];
        my $hid = $f->hseqname or $self->throw("Missing hid");;
        next unless $f->strand == $strand;
        $hit_features->{$hid} ||= [];
        push(@{$hit_features->{$hid}}, $f);
    }
    
    my( $is_linear );
    if ($strand == 1) {
        $is_linear = sub {
            my( $x, $y ) = @_;
            
            return $x->hstart < $y->hstart;
        };
    } else {
        $is_linear = sub {
            my( $x, $y ) = @_;
            
            return $x->hstart > $y->hstart;
        };
    }
    
    my $count = 0;
    while (my ($hid, $flist) = each %$hit_features) {
        $count++;
        # Sort features by start/end in genomic
        @$flist = sort {$a->start <=> $b->start or $a->end <=> $b->end} @$flist;
        
        # Group into linear matches
        my @sets = ([$flist->[0]]);
        my $curr = $sets[0];
        for (my $i = 1; $i < @$flist; $i++) {
            my $prev = $flist->[$i - 1];
            my $this = $flist->[$i];
            if (&$is_linear($prev, $this)) {
                push(@$curr, $this);
            } else {
                $curr = [$this];
                push(@sets, $curr);
            }
        }
        
        foreach my $lin (@sets) {
            $self->do_mini_est_genome($lin);
        }
        last if $count >= 10;
    }
}

sub do_mini_est_genome {
    my( $self, $linear ) = @_;
    
    my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome(
        '-genomic'    => $self->unmasked,
        '-features'   => $linear,
        '-seqfetcher' => $self->seqfetcher,
        );
    $e2g->run;
    
    # Could possibly get matches out of MiniEst2Genome on
    # opposite strand to the features fed into it!
    # Should fix properly in MiniEst2Genome a la Finished_Est2Genome
    my $est_strand = $linear->[0]->strand;

    foreach my $e2g_gene ($e2g->output) {
        my @exons = $e2g_gene->sub_SeqFeature;
        foreach my $exon (@exons){
            my @supp_evidence = $exon->sub_SeqFeature;
            #display(@sub_Feats);                    
            foreach my $sp (@supp_evidence) {

                # Fix strand
                $sp->strand($est_strand);

                # source_tag and primary_tag have to be set to
                # something, or validate method in FeaturePair
                # (callled by RunnableDB) throws!
                $sp->feature2->source_tag ('I_am_valid');
                $sp->feature2->primary_tag('I_am_valid');
                
                $self->add_output($sp);
            }
        }
    }
}

sub add_output {
    my( $self, @feat ) = @_;
    
    my $ana = $self->analysis;
    foreach my $f (@feat) {
        $f->analysis($ana);
    }
    
    $self->{'_output'} ||= [];
    push(@{$self->{'_output'}}, @feat);
}

sub output {
    my( $self ) = @_;
    
    if (my $out = $self->{'_output'}) {
        return @$out;
    } else {
        return;
    }
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

