#
# Written by Stephen Keenan
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FindPrimer

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $find_primer = Bio::EnsEMBL::Pipeline::RunnableDB::FindPrimer->new(
    -dbobj      => $db,
    -input_id   => $input_id
    -analysis   => $analysis
);
$find_primer->fetch_input();
$find_primer->run();
$find_primer->output();
$find_primer->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::FindPrimer to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction of
parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is required
for databse access.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::FindPrimer;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::FindPrimer;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::FindPrimer object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::FindPrimer object
    Args    :     -dbobj     A Bio::EnsEMBL::DB::Obj, 
                  -input_id  Contig input id , 
                  -analysis  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    $self->throw("Analysis object required") unless ($self->analysis);
    
    
    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::FindPrimer');
    return $self;
}

=head2 fetch_input

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);
    
    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);
    my $genseq    = $contig->primary_seq()
     or $self->throw("Unable to fetch contig");
    $self->genseq($genseq);
    
    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::FindPrimer->new(
        '-clone'   => $genseq,
        '-analysis'   => $self->analysis,
    );
    
    $self->runnable($runnable);

}

#get/set for runnable and args
sub runnable {
   my ( $self, $arg ) = @_;

    if ( defined($arg) ) {
        $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI")
          unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");

        $self->{_runnable} = $arg;
    }

    return $self->{_runnable};
}


__END__


=head1 AUTHOR

Stephen Keenan B<email> keenan@sanger.ac.uk

1;
