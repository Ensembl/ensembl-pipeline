#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
# RunnableDB::EPCR modified by Simon Potter to run on clones
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_EPCR

=head1 SYNOPSIS

  my $db     = Bio::EnsEMBL::DBLoader->new($locator);
  my $epcr   = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_EPCR->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );
  $epcr->fetch_input();
  $epcr->run();
  $epcr->output();
  $epcr->write_output(); #writes to DB

=head1 DESCRIPTION

This object wraps Bio::EnsEMBL::Pipeline::Runnable::EPCR to add
functionality for reading and writing to databases.
This object takes clone ids, while Bio::EnsEMBL::Pipeline::RunnableDB::EPCR
acts on contigs. This allows us to submit one Job per clone rather
than one per contig and should speed things up ...

The appropriate Bio::EnsEMBL::Pipeline::Analysis object must be passed for
extraction of appropriate parameters. A Bio::EnsEMBL::Pipeline::DBSQL::Obj is
required for databse access.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_EPCR;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::EPCR;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB::EPCR);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Clone_EPCR object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Clone_EPCR object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:  Clone input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);    
    $self->{'_fplist'}      = [];
    $self->{'_runnable'}    = [];
        
    $self->throw("Analysis object required") unless ($self->analysis);
    
    return $self;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    my $cloneid  = $self->input_id;
    my $clone    = $self->dbobj->get_Clone($cloneid);
    foreach my $contig ($clone->get_all_Contigs()) {
	my $genseq = $contig->get_repeatmasked_seq() or $self->throw("Unable to fetch contig");
	$self->runnable($genseq);
    }
}

#get/set for runnable and args
sub runnable {
    my ($self, $genseq) = @_;
    
    if ($genseq)
    {
        #extract parameters into a hash
        my ($parameter_string) = $self->analysis->parameters() ;
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
        $parameters {'-db'}      = $self->analysis->db_file();  
        $parameters {'-clone'}   = $genseq;
        push (@{$self->{'_runnable'}}, 
	      Bio::EnsEMBL::Pipeline::Runnable::EPCR->new(%parameters));
    }
    return @{$self->{'_runnable'}};
}

=head2 run
          
    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::CPG->run()
    Returns :   none                                             
    Args    :   none                                             
                    
=cut

sub run {
    my ($self) = @_;
    $self->throw("Runnable modules not set") unless ($self->runnable());
    foreach my $runnable ($self->runnable())
    {
        $runnable->run();                   
    }                    
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::CPG->output() 
                on each runnable  
    Returns :   An array of SeqFeatures representing cpg islands
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    my @output;
    foreach my $runnable ($self->runnable())
    {
	push (@output, $runnable->output());
    }   
    return @output;
}   

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetches output data from a frozen perl object
                stored in file $file_name
    Returns :   array of repeats (with start and end)
    Args    :   none

=cut

#cut & pasted from Vert_Est2genome
sub fetch_output {
    my($self,$output) = @_;
    
    $output || $self->throw("No frozen object passed for the output");
    
    my $object;
    open (IN,"<$output") || do {print STDERR ("Could not open output data file... skipping job\n"); next;};
    
    while (<IN>) 
    {
        $_ =~ s/\[//;
	    $_ =~ s/\]//;
	    $object .= $_;
    }
    close(IN);
    my @out;
   
    if (defined($object)) 
    {
        my (@obj) = FreezeThaw::thaw($object);
        foreach my $array (@obj) 
        {
	        foreach my $object (@$array) 
            {
	            push @out, $object;
	        }
        }
    }
    return @out;
}

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
	unless ($self->runnable());

    my $db=$self->dbobj();
    foreach my $runnable ($self->runnable())
    {
	my $contig;
	my @features = $runnable->output();
	eval 
	{
          $contig = $db->get_Contig($runnable->clone->display_id());
	};
	if ($@)
	{
	    print STDERR "Contig not found, skipping writing output to db: $@\n";
	}
	elsif (@features)
	{
            foreach my $f (@features) {
                $f->seqname($contig->internal_id);
            }
	    print STDERR "Writing features to database\n";
	    my $feat_adp = Bio::EnsEMBL::DBSQL::FeatureAdaptor->new($db);
	    $feat_adp->store($contig, @features);
	}
    }
    return 1;
}

1;
