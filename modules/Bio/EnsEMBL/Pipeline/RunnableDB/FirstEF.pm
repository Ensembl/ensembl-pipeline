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
use Bio::EnsEMBL::Pipeline::Config::FirstEF;
use Bio::EnsEMBL::Pipeline::Config::FirstEF qw(  );

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
    
    # Check analysis exists
    $self->throw("Analysis object not passed to RunnableDB/FirstEF") 
      unless $self->analysis->isa("Bio::EnsEMBL::Analysis");

    # Check input id
    $self->throw("No input id passed to RunnableDB/FirstEF (eg. 1.500000-1000000)") 
      unless $self->input_id;
    
    # Fetch slice specified by input id
    unless ($self->input_id =~ /$FEF_INPUTID_REGEX/ ){
      $self->throw("Input id [".$self->input_id."] not compatible with ".
		   "FEF_INPUTID_REGEX [$FEF_INPUTID_REGEX]");
    }
    my $chr   = $1;
    my $start = $2;
    my $end   = $3;

    my $slice = $self->db->get_SliceAdaptor->fetch_by_chr_start_end($chr,$start,$end);
    $self->query($slice);

#print ">thing\n" . $self->query->get_repeatmasked_seq->seq . "\n";

    my $runnable = Bio::EnsEMBL::Pipeline::Runnable::FirstEF->new(
		     -query            => $self->query,
		     -repeatmasked_seq => $self->query->get_repeatmasked_seq,
		     -analysis         => $self->analysis,
		     -firstef_dir      => $FEF_APPLICATION_DIR,
		     -param_dir        => $FEF_PARAMETER_DIR);

							   
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
	$f->contig($slice);
	my @mapped = $f->transform;

        if (@mapped == 0) {
	    $self->warn("Couldn't map $f - skipping");
	    next;
        }
        if (@mapped == 1 && $mapped[0]->isa("Bio::EnsEMBL::Mapper::Gap")) {
	    $self->warn("$f seems to be on a gap - something bad has happened ...");
	    next;
        }

	push @mapped_features, $mapped[0];

    }
    $sfa->store(@mapped_features) if @mapped_features;

    return 1;
}



1;
