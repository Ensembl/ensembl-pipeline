#Cared for by Tania & Kailan <tania/chris@fugu-sg.org>
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Synteny.pm

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Synteny->new(-input_id => $input_id
                                             -dbobj     => $db,
                                             -input_id  => $id
                                             );
    $obj->fetch_input
    $obj->run

    my @syntenyclusters = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Alan Chris. (Kailan)
Tania Oh
 
=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Syntenylite;

BEGIN {
  require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}


use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Syntenylite;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);



=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Synteny object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Synteny object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    $self->{'_syntenyhits'} = undef; #create key to store the aligned blocks 
    $self->throw("Analysis object required") unless (defined($self->analysis));

    return $self;
}


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   fetches protein hits from the compara db  
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_;


    $self->throw("No input id") unless defined($self->input_id);
    
    my ($c1,$c2);
    #if ($input_id =~ /(\S+)::(\S+)/ ) {
    if ($self->input_id =~ /(\S+)::(\S+)/ ) {
        $c1 = $1;
        $c2 = $2;
    }
    else {
        $self->throw("Input id not in correct format: got $self->input_id, should be parsable by (\w+)\:(\S+)\/(\w+)\:(\S+)");
    }

    $self->_c1_id($c1);
    $self->_c2_id($c2);

    #get the Synteny dbadaptor from the current ensembl dbadaptor
    my $newadp = $self->dbobj->get_ComparaDBAdaptor();
    my $hitadap = $newadp->get_SyntenyHitAdaptor();
    my $prot_hits_arrref = $hitadap->fetch_protein_match_fd_on_fugu_n_chr($c1,$c2); 
=head 
    if (!$prot_hits_arrref) {
    #if there are no hits
        exit 0;
    }
=cut
    my @prot_hits = @$prot_hits_arrref; 
    $self->protein_matches(@prot_hits);
    #might be a bad idea, but saves typing
    $self->adaptor ($hitadap);
    
}


=head2 protein_matches 


    Title   :	protein_matches
    Usage   :	$self->protein_matches(@protein_id);
    Function:
    Returns :	A simple  array of protein ids
    Args    :	none

=cut

sub protein_matches{
    my ($self, @proteins) = @_;

    if (@proteins)
    {
        
        #lazy match, just match the first ones 
        $proteins[0]->isa("Bio::EnsEMBL::Compara::Synteny_Hit") or
		$self->throw("there are no proteins available-in protein_matches, Synteny!!");
	push (@{$self->{'_protein_matches'}}, @proteins);
    }
    return @{$self->{'_protein_matches'}};



}





sub runnable {
    my ($self, @runnable) = @_;
    if (@runnable)
    {
        foreach my $runnable (@runnable)
        {
            $runnable->isa("Bio::EnsEMBL::Pipeline::RunnableI") or
                $self->throw("Input to runnable is not Bio::EnsEMBL::Pipeline::RunnableI");
        }
        push (@{$self->{'_runnable'}}, @runnable);
    }
    return @{$self->{'_runnable'}};
}

=head2 run

    Title   :   run
    Usage   :   $self->run();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Blast->run()
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self) = @_;
   #testing 
=head
     my $compute = scalar ($self->protein_matches); 
     my $computer = $compute/2;
=cut

    #coz it's 2X the no. of proteins
    if ((scalar ($self->protein_matches)) / 2  > 1 ){
   
        
        my @prot_objs = $self->protein_matches;
        my $prot_objs = \@prot_objs;

        #get the Synteny_hit objects
        my $contig = $self->_c1_id;
        my $chr = $self->_c2_id;
        
        my $runnable = Bio::EnsEMBL::Pipeline::Runnable::Syntenylite->new(-prot_hits => $prot_objs);
        $runnable->run();
        $self->runnable($runnable);
    } else {exit 0;} # should we return something??
}

=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::Synteny->output()
    Returns :   An array of Bio::EnsEMBL::Compara::Synteny_Hits objects (same as teh compara objects)
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    my @output;
    foreach my $run ($self->runnable) {
      push(@output,$run->output);
    }
    return @output;
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

    my $synclusadp= $self->dbobj->get_ComparaDBAdaptor->get_SyntenyClusterAdaptor;
    foreach my $op ($self->output){
        my $syntenycluster= $op;
        $synclusadp->store($syntenycluster);
   }
   

    return 1;
}

=head2 _c1_id

 Title   : _c1_id
 Usage   : $obj->_c1_id($newval)
 Function: Getset for _c1_id value
 Returns : value of _c1_id
 Args    : newvalue (optional)


=cut

sub _c1_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_c1_id'} = $value;
    }
    return $obj->{'_c1_id'};

}

=head2 _c2_id

 Title   : _c2_id
 Usage   : $obj->_c2_id($newval)
 Function: Getset for _c2_id value
 Returns : value of _c2_id
 Args    : newvalue (optional)


=cut

sub _c2_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_c2_id'} = $value;
    }
    return $obj->{'_c2_id'};

}

=head2 adaptor 

 Title   : adaptor 
 Usage   : $obj->adaptor($newval)
 Function: Getset for adaptor value
 Returns : value of adaptor 
 Args    : newvalue (optional)


=cut

sub  adaptor{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'adaptor'} = $value;
    }
    return $obj->{'adaptor'};

}



1;

