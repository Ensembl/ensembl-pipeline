#
#
# Author Thomas Down <td2@sanger.ac.uk>
#
# Based on CPG.pm by Val Curwen
# Modified by SCP to run on VC's
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::VC_EponineTSS

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);

my $eponine = Bio::EnsEMBL::Pipeline::RunnableDB::VC_EponineTSS->new ( 
                                   -db          => $db,
			           -input_id   => $input_id
                                   -analysis   => $analysis 
                                    );

$eponine->fetch_input();

$eponine->run();

$eponine->output();

$eponine->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::EponineTSS to add
functionality to read and write to databases. The appropriate
Bio::EnsEMBL::Analysis object must be passed for extraction
of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
is required for database access.

=head1 CONTACT

For general Ensembl comments mail to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::VC_EponineTSS;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::EponineTSS;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

 Title   : new
 Usage   : $self->new(-DBOBJ       => $db,
                      -INPUT_ID    => $id,
                      -ANALYSIS    => $analysis
	   );
                           
 Function: creates a Bio::EnsEMBL::Pipeline::RunnableDB::VC_EponineTSS object
 Returns : Bio::EnsEMBL::Pipeline::RunnableDB::VC_EponineTSS
 Args    :      -db      :     A Bio::EnsEMBL::DBSQL::DBAdaptor, 
                -input_id:  A virtual contig ('chr_name.start.end')
                -analysis:  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_fplist'}      = []; # ???   
    $self->{'_genseq'}      = undef;
    $self->{'_runnable'}    = undef;
    
    $self->throw("Analysis object required") unless ($self->analysis);
    
    $self->runnable('Bio::EnsEMBL::Pipeline::Runnable::EponineTSS');
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for Eponine from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my ($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $vc_str = $self->input_id;
    my ($chr, $start, $end, $sgp) =
     $vc_str =~ m!(\S+)\.(\d+)\.(\d+):?([^:]*)!;

    $self->db->assembly_type($sgp) if $sgp;

    my $slice = $self->db->get_SliceAdaptor->
     fetch_by_chr_start_end($chr, $start, $end);

    my $genseq = $slice->primary_seq() or
     $self->throw("Unable to fetch virtual contig");

    $self->vcontig($slice);
    $self->genseq($genseq);
}

#get/set for runnable and args
sub runnable {
    my ($self, $runnable) = @_;

    if ($runnable)
    {
        #extract parameters into a hash
        my ($parameter_string) = $self->parameters() ;
        my %parameters;
        if ($parameter_string)
        {
            $parameter_string =~ s/\s+//g;
            my @pairs = split (/,/, $parameter_string);
            
            foreach my $pair (@pairs)
            {
                my ($key, $value) = split (/=>/, $pair);
                $parameters{$key} = $value;
            }
        }
        $parameters{'-java'} = $self->analysis->program_file;
        #creates empty Bio::EnsEMBL::Runnable::EponineTSS object
        $self->{'_runnable'} = $runnable->new(
	      '-threshold' => $parameters{'-threshold'},
	      '-epojar'    => $parameters{'-epojar'},
	      '-java'      => $parameters{'-java'}
	);
    }
    return $self->{'_runnable'};
}


sub write_output {
    my($self) = @_;
    my $contig;

    my $db  = $self->db;
    my $sfa = $db->get_SimpleFeatureAdaptor;

    my $vc = $self->vcontig;
    my @mapped_features;

    foreach my $f ($self->output) {

	$f->analysis($self->analysis);
	$f->contig($vc);

	my (@mapped) = $f->transform;

	unless (@mapped == 1) {
	    print STDERR "Warning: can't map $f";
	    next;
	}

	push @mapped_features, $mapped[0];

	my $contig_id = $f->contig->dbID;
    }

    $sfa->store(@mapped_features);

    return 1;
}

1;
