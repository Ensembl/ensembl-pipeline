#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Blast

=head1 SYNOPSIS

my $db      = Bio::EnsEMBL::DBLoader->new($locator);
my $repmask = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Blast->new ( 
                                                    -dbobj      => $db,
			                                        -input_id   => $input_id
                                                    -analysis   => $analysis );
$repmask->fetch_input();
$repmask->run();
$repmask->output();
$repmask->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::Blast to add functionality 
for reading and writing to databases. This objects requires clone ids, 
Bio::EnsEMBL::Pipeline::RunnableDB::Blast acts on contigs.
The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Blast;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::Blast);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Args    :    -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_runnable'}    = [];
    
    my ( $threshold) = $self->_rearrange ([qw(THRESHOLD)], @args);    

    $self->throw("Analysis object required") unless ($self->analysis);
    
#    if ($threshold)
#    {
#        $self->threshold($threshold);
#    }
#    else
#    {
#        $self->threshold(1e-6);
#    }
    return $self;
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($obj);
    Function:   Gets or sets the value of dbobj
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=head2 threshold

    Title   :   threshold
    Usage   :   $obj->threshold($value);
    Function:   Get/set method for threshold score required for writing
                Feature/FeaturePair to database.
    Args    :   Optional value (depends on type of Analysis)

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for blast from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    print STDERR "Input id " . $self->input_id . "\n";
    print STDERR "Dbobj " . $self->dbobj . "\n";
    
    my $cloneid     = $self->input_id;
    my $clone       = $self->dbobj->get_Clone($cloneid);
    foreach my $contig  ($clone->get_all_Contigs())
    {       
        my $genseq = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");
	print STDERR "Passing in genseq $genseq\n";
        $self->runnable($genseq);
    }
}

=head2 runnable

    Title   :   runnable
    Usage   :   $self->runnable
    Function:   creates Blast runnable
    Returns :   none
    Args    :   none

=cut

sub runnable {
    my ($self, $genseq) = @_;
    if ($genseq)    {
	my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (   
					-query    => $genseq,
					-program  => $self->analysis->program,
					-database => $self->analysis->db,
					-threshold => 1, 
								     );
	
	push (@{$self->{'_runnable'}}, $blast);
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self,$dir) = @_;
    $self->throw("Runnable modules not set") unless ($self->runnable());
    foreach my $runnable ($self->runnable)
    {
#        $runnable->threshold($self->threshold);
        $runnable->run($dir);
    }
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::RepeatMasker->output()
    Returns :   An array of Bio::EnsEMBL::Repeat objects (FeaturePairs)
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    
    my @output;
    foreach my $runnable ($self->runnable)
    {
        push (@output, $runnable->output());
    }
    return @output;
}

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetchs output data from a frozen perl object
                stored in file $file_name
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self) = @_;
    
    $self->throw("fetch_input must be called before write_output\n") 
        unless ($self->runnable);

    my $db=$self->dbobj();
    foreach my $runnable ($self->runnable)
    {
        my $contig;
        my @features = $runnable->output();
        eval 
        {
	        $contig = $db->get_Contig($runnable->clone->display_id);
        };
        if ($@) 
        {
	        print STDERR "Contig not found, skipping writing output to db\n";
        }
        elsif (@features) 
        {
	        foreach my $feature (@features)
            {
                print STDERR ($feature->hseqname()."\t");
                $feature->analysis($self->analysis);
            }
            my $feat_adp=Bio::EnsEMBL::DBSQL::FeatureAdaptor->new($db);
	        $feat_adp->store($contig, @features);
        }
        return 1;
    } 
}

1;
