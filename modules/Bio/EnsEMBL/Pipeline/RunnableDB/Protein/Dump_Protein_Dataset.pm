
#
# BioPerl module for Prints.pm
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Prints.pm - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Dump_Protein_Dataset;


use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::DBSQL::Protein_Adaptor;
use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;





@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => @id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Protein::Prints object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Blast object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                -input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Pipeline::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->{'_fplist'}      = [];
    $self->{'_pepseq'}      = undef;
    $self->{'_runnable'}    = undef;            
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
     open (OUT,">>".$self->pep_file);

    $self->throw("No input id") unless defined($self->input_id);
    my $prot_adapt = $self->dbobj->get_Protein_Adaptor();

    my @ids = split(/,/,$self->input_id);

    foreach my $id (@ids) {
	print STDERR "ID: $id\n";
    	my $prot = $prot_adapt->fetch_Protein_by_dbid($id);
	my $pepseq    = $prot->seq;
	
	print STDERR "SEQ: $pepseq\n";

	$self->genseq($pepseq);
	print OUT ">".$id."\n$pepseq\n";
	
    }
}



sub run {
    my($self) = @_;
    return 1;
}


sub write_output {
    my($self) = @_;
    return 1;
}

=head2 pep_file

 Title   : pep_file
 Usage   : $obj->pep_file($newval)
 Function: 
 Returns : value of pep_file
 Args    : newvalue (optional)


=cut

sub pep_file{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'pep_file'} = $value;
    }
    return $obj->{'pep_file'};

}


1;










