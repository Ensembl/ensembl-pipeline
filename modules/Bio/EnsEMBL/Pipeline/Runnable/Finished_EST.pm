
### Bio::EnsEMBL::Pipeline::Runnable::Finished_EST

package Bio::EnsEMBL::Pipeline::Runnable::Finished_EST;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = ('Bio::EnsEMBL::Pipeline::RunnableI');

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($query,
        $database,
        $program,
        $options,
        $threshold,
        $thres_type,
        ) = $self->_rearrange([qw{
            query
            database
            program
            options
            threshold
            threshold_type
            }], @args);
    
    die "No QUERY (genomic sequence) given" unless $query;

    $self->clone            ($query);
    $self->database         ($database);
    $self->program          ($program);
    $self->options          ($options);
    $self->threshold        ($threshold);
    $self->threshold_type   ($thres_type);
    
    return $self;
}

sub database {
    my( $self, $database ) = @_;
    
    if ($database) {
        $self->{'_database'} = $database;
    }
    return $self->{'_database'} || 'dbEST';
}

sub program {
    my( $self, $program ) = @_;
    
    if ($program) {
        $self->{'_program'} = $program;
    }
    return $self->{'_program'} || 'wublastn';
}

sub options {
    my( $self, $options ) = @_;
    
    if ($options) {
        $self->{'_options'} = $options;
    }
    return $self->{'_options'} || 'Z=500000000';
}

sub threshold {
    my( $self, $threshold ) = @_;
    
    if ($threshold) {
        $self->{'_threshold'} = $threshold;
    }
    return $self->{'_threshold'} || 1e-4;
}

sub threshold_type {
    my( $self, $threshold_type ) = @_;
    
    if ($threshold_type) {
        $self->{'_threshold_type'} = $threshold_type;
    }
    return $self->{'_threshold_type'} || 'PVALUE';
}

sub _make_blast_paramters {
    my( $self ) = @_;
    
    my @param;
    foreach my $meth (qw{
        database
        program
        options
        threshold
        threshold_type
        })
    {
        push(@param, "-$meth", $self->$meth());
    }
    return @param;
}

sub run {
    my( $self ) = @_;
    
    my $seq = $self->clone;
    my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast
        ->new($self->_make_blast_paramters);
    $blast->run;
    my @features = $blast->fetch_output;
}

sub output {
    
}


1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::Runnable::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

