#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer->new(
					     -input_id => $id,
					     -score => $score
								    );
    $obj->fetch_input();
    $obj->run();

    my @newfeatures = $obj->output();
    
    $obj->write_output();

=head1 DESCRIPTION

    This runnabledb creates the crossmatch runnable needed to 
    compare two contigs from databases from different organisms

    This runnabledb is more complex than usual in that it needs
    to be connected to several dbs:

    dbobj: this is where the output is finally written, 
    storing the hits between two Ensembl databases. 
    It is actually a Bio::EnsEMBL::Comparer::DBSQL::DBAdaptor 
    database, not a normal pipeline EnsEMBL database

    query_db and target_db: these are the two databases from which the 
    input dna is fetched

    The input id has the following format: 
    db_name:contig_id/db_name:contig_id

    min_score is an argument for the crossmatch runnable to set the 
    minimum score used in the crossmatch run

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
use Bio::PrimarySeq;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::FeaturePair;
use Bio::SeqIO;
use Bio::Root::RootI;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($score) = $self->_rearrange([qw(MIN_SCORE)],@args);
  bless $self, $class;
  $score ||= 100;
  $self->min_score($score);
  return $self; 
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input dna for crossmatch from the databases
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_; 

    $self->throw("No input id") unless defined($self->input_id);

    my $input_id  = $self->input_id;

    my ($db1,$db2,$c1,$c2);
    if ($input_id =~ /(\S+):(\S+)::(\S+):(\S+)/ ) {
	$db1 = $1;
	$c1 = $2;
	$db2 = $3;
	$c2 = $4;
    }
    else {
	$self->throw("Input id not in correct format: got $input_id, should be parsable by (\w+)\:(\S+)\/(\w+)\:(\S+)");
    }

    $self->_c1_id($c1);
    $self->_c2_id($c2);

    my $gadp = $self->dbobj->get_GenomeDBAdaptor();

    $db1 = $gadp->fetch_by_species_tag($db1)->ensembl_db();
    $db2 = $gadp->fetch_by_species_tag($db2)->ensembl_db();

    my $contig1 = $db1->get_Contig($c1);
    my $contig2 = $db2->get_Contig($c2);

    my $seq1 = Bio::PrimarySeq->new( -display_id => 'seq1', -seq => $contig1->seq);
    my $seq2 = Bio::PrimarySeq->new( -display_id => 'seq2', -seq => $contig2->seq);
    
    my $cross = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( 
								   -nocopy => 1,
								   -seq1 => $seq1,
								   -seq2 => $seq2,
								   -score => $self->min_score,
								   -minmatch => 14,
								   -masklevel => 80
								   );
    $self->runnable($cross);
}

=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the crossmatch analysis, producing feature pairs
    Returns :   Nothing, but $self->output is filled 
                contains the crossmatch feature pairs.
    Args    :   None

=cut

sub run {
    my ($self) = @_;

    $self->throw("Can't run - no runnable objects") unless $self->runnable;
    $self->runnable->run();
    $self->output($self->runnable->output);
}

=head2 output

    Title   :   output
    Usage   :   $self->output() or $self->output(@)
    Function:   Push an array of Bio::EnsEMBL::FeaturePairs if Args is given and
                Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   Array of Bio::EnsEMBL::FeaturePairs (optional)

=cut

sub output {
   my ($self,@args) = @_;
   if(@args) {
     push @{$self->{'output'}}, @args;
    }
    return @{$self->{'output'}};
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->{_output} into $self->dbobj
    Returns :   1
    Args    :   None

=cut

sub write_output {
    my ($self) = @_;

    my @features = $self->output;
    foreach my $f (@features) {
#        print "Got feature $f\n";

	print $f->seqname."\t".$f->start."\t".$f->end."\t".$f->strand."\tscore:".$f->score."\t".$f->hseqname."\t".$f->hstart."\t".$f->hend."\t".$f->hstrand."\n";
    }
    my $db = $self->dbobj();
    print STDERR "Going to write to ".$db->dbname."\n";

# Use write.t as an example for building the genomic align datastructure
# make sure you write back dnafrag objects first into the database.



    $self->throw("Have not delt with write back yet");

#    foreach my $f (@features) {
#		    my $contig = $contig_hash{$c};
#	    my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
#	    $feat_Obj->write($contig, @features);
#	}
#    }
#    return 1;
}

=head2 runnable

 Title   : runnable
 Usage   : $obj->runnable($runnable)
 Function: get/set for runnable
 Returns : path to runnable
 Args    : runnable (optional)


=cut

sub runnable {
   my ($self, $runnable) = @_;

   if( defined $runnable ) {
       $self->{'_runnable'} = $runnable;
   }
   return $self->{'_runnable'};
}


=head2 min_score

 Title   : min_score
 Usage   : $obj->min_score($newval)
 Function: Getset for min_score value
 Returns : value of min_score
 Args    : newvalue (optional)


=cut

sub min_score{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'min_score'} = $value;
    }
    return $obj->{'min_score'};

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


1;
