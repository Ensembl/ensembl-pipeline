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

Bio::EnsEMBL::Pipeline::RunnableDB::SSAHA

=head1 SYNOPSIS

my $db    = Bio::EnsEMBL::DBLoader->new($locator);
my $ssaha = Bio::EnsEMBL::Pipeline::RunnableDB::SSAHA->new(
    -dbobj      => $db,
    -input_id   => $input_id
    -analysis   => $analysis
);

$ssaha->fetch_input();
$ssaha->run();
$ssaha->output();
$ssaha->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::SSAHA to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction of
appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for databse access.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::SSAHA;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::SSAHA;

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    return $self;
}


sub fetch_input {
    my( $self) = @_;

    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->db->get_Contig($contigid);
    my $genseq    = $contig->primary_seq() or
     $self->throw("Unable to fetch contig");
    $self->genseq($genseq);
}


sub runnable {
    my ($self, $runnable) = @_;

    #need to pass one peptide at a time
    $self->throw("Input must be fetched before run") unless ($self->genseq);

    if ($runnable) {

        #extract parameters into a hash
        my ($parameter_string) = $self->analysis->parameters();
        my %parameters;
        my ($thresh, $thresh_type, $arguments);

        if ($parameter_string) {

            $parameter_string =~ s/\s+//g;
            my @pairs = split (/,/, $parameter_string);
            foreach my $pair (@pairs)
            {
                my ($key, $value) = split (/=>/, $pair);
                if ($key eq '-threshold_type' && $value) {
                    $thresh_type = $value;
                }
                elsif ($key eq '-threshold' && $value) {
                    $thresh = $value;
                }
                else
	        # remaining arguments not of '=>' form
	        # are simple flags (like -p1)
                {
                    $arguments .= " $key ";
                }
            }
        }

        $parameters{'-query'} = $self->genseq;
        $parameters{'-database'} = $self->analysis->db_file;
        $parameters{'-ssaha'} = $self->analysis->program;

        my $runnable = Bio::EnsEMBL::Pipeline::Runnable::SSAHA->new(
	    %parameters
        );
    }

    return $self->runnable();

}


sub output {
    my ($self) = @_;

    my @runnable = $self->runnable;
    my @results;

    foreach my $runnable(@runnable){
      print STDERR "runnable = ".$runnable[0]."\n";
      push(@results, $runnable->output);
    }
    return @results;
}

1;
