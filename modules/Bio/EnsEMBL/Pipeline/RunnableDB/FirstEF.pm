#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB:FirstEF

=head1 SYNOPSIS

my $fef = Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF->new(-input_id  => $input_id,
                                                           -analysis  => $analysis );
$fef->fetch_input;
$fef->run;
my @output = $fef->output;
$fef->write_output;

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::FirstEF to add
functionality for reading and writing to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::FirstEF;
use Bio::EnsEMBL::Pipeline::Config::General qw(SLICE_INPUT_ID_REGEX
					       PIPELINE_WORK_DIR);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(_things_)
    Function:   Creates a new RunnableDB/FirstEF object
    Returns :   Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF
    Args    :   lots

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    return $self; 
}


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;

    $self->fetch_sequence;

    my %parameters = $self->parameter_hash;

    my $repeatmasked = 1 if defined($parameters{'-repeatmasked'});

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::FirstEF->new(
		     -query        => $self->query,
		     -repeatmasked => $repeatmasked,
		     -analysis     => $self->analysis,
		     -firstef_bin  => $self->analysis->program_file,
		     -param_dir    => $parameters{'-parameters_dir'},
		     -parse_script => $parameters{'-parse_script'},
		     -work_dir     => $PIPELINE_WORK_DIR);

							
    $self->runnable($runnable);

    return 1;
}


sub write_output {
    my($self) = @_;

    my $db  = $self->db;
    my $sfa = $self->db->get_SimpleFeatureAdaptor;

    my @mapped_features;

    my $slice = $self->query;

    foreach my $f ($self->output) {

      $f->analysis($self->analysis);
      $f->slice($slice);
    }
    $sfa->store($self->output);

    return 1;
}



1;
