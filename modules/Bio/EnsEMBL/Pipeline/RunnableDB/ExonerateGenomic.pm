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

Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Exonerate->new(
					     -dbobj     => $db,
					     -input_id  => $id,
					     -analysis   => $analysis 		 
                                             );
    $obj->fetch_input();
    $obj->run();

    my @newfeatures = $obj->output();
    
    $obj->write_output();

=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateGenomic;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::FeaturePair;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;

    my ($db,$query,$target,$exonerate,$exargs,$max_query,$max_target) = $self->_rearrange([qw(DB QUERYDB TARGETDB EXONERATE EXARGS MAXQUERY MAXTARGET)],@args);
    
    $query || $self->throw("Did not specify query EnsEMBL db");
    $target || $self->throw("Did not specify target EnsEMBL db");

    $self->db($db);
    $self->exonerate($exonerate) if defined $exonerate;
    $self->querydb($query);
    $self->targetdb($target);
    #Apparently exonerate likes whole genome on one side, 
    #500 Kb on the other
    $max_query ||=  1000000000;
    $max_target ||= 500000;
    $self->max_query($max_query);
    $self->max_target($max_target);
    
    #Thanks God, Val hardcoded it ;)
    $exargs = " -w 14 -t 65 -H 100 -D 15 -m 500 " unless defined $exargs;
    $self->exonerate_args($exargs) if defined $exargs;
    
    return $self; 
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for exonerate from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_; 
  
    my @query_files = $self->_split_db_to_files('query');
    my @target_files = $self->_split_db_to_files('target');

    #Create Exonerate Runnable between all 
    foreach my $qf (@query_files) {
	foreach my $tf (@target_files) {
	    $self->throw("Can't run Exonerate without both genomic and EST sequences") 
		unless (scalar(@query_files) && scalar(@target_files));
	    #my $executable =  $self->analysis->program_file();
	    my $exonerate = new Bio::EnsEMBL::Pipeline::Runnable::Exonerate(
			    '-genomic'  => $qf,
			    '-est'      => $tf,
			    '-args'     => $self->exonerate_args,
			    );
	    push(@{$self->{'_runnables'}},$exonerate);
	}
    }
    return 1;
}


=head2 _split_db_to_files

 Title   : _split_db_to_files
 Usage   :
 Function: Split query and target db sequences in chunks up to max 
           query/target length and write them to file, keep array 
           of filenames 
 Example :
 Returns : 
 Args    : 'query' or 'target'


=cut

sub _split_db_to_files {
    my ($self,$switch) = @_;

    my @genseqfiles;
    my $length=0;
    my $c=1;
    my $filename = "/tmp/exonerate_".$switch."_$c.fa"; 
    my $seqout = Bio::SeqIO->new(-file => ">$filename" , '-format' => 'Fasta');
    my ($db,$max);
    #Does not need input id, because it has to compare all against all
    if ($switch eq 'query') {
	$db = $self->querydb;
	$max = $self->max_query;
    }
    else {
	$db = $self->targetdb;
	$max = $self->max_target
    }
    my @contigids = $db->get_all_Contig_id();
    while (my $cid = shift (@contigids)) {
	my $contig;
	eval {
	    $contig = $db->get_Contig($cid);
	};
	if ($@) {
	    print STDERR "Couldn't fetch contig $@\n";
	    next;
	}

	$length += $contig->length;
	if ($length > $max) {
	    push (@genseqfiles,$filename);
	    $c++;
	    $length = $contig->length;
	    $filename = "/tmp/exonerate_".$switch."_$c.fa";
	    $seqout = Bio::SeqIO->new(-file => ">$filename" , '-format' => 'Fasta');
	    $seqout->write_seq($contig->get_repeatmasked_seq());
	    last;
	}
	else {
	    $seqout->write_seq($contig->get_repeatmasked_seq());
	} 
    }
    push (@genseqfiles,$filename);

    return @genseqfiles;
}

=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the exonerate analysis, producing Bio::EnsEMBL::Gene predictions
    Returns :   Nothing, but $self{_output} contains the predicted genes.
    Args    :   None

=cut

sub run {
    my ($self) = @_;

    #Since this method, called once, runs all runnables, and generates 
    #all feature pairs, it is necessary to call write_output at the end 
    #of each run, rather than waiting to accumulate all output, and then 
    #writing it to the db

    #It is not orthodox, any better ideas?
    #(The right way of doing it is commented out)

    $self->throw("Can't run - no runnable objects") unless scalar(@{$self->{'_runnables'}});
    my @output;
    foreach my $runnable (@{$self->{'_runnables'}}) {
	$runnable->run('ungapped');
	my @fps = $runnable->output;
	#push (@output,@fps);
	$self->write_output(\@fps);
    }
    #$self->output(\@output);
}

=head2 output

    Title   :   output
    Usage   :   $self->output()
    Function:   Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   None

=cut

sub output{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'output'} = $value;
    }
    return $obj->{'output'};

}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->{_output} into $self->dbobj
    Returns :   1
    Args    :   None

=cut

sub write_output {
    my($self,$out) = @_;

    #foreach my $f (@$out) {
#	print $f->seqname."\texonerate\tsimilarity\t".$f->start."\t".$f->end."\t".$f->score."\t".$f->strand."\t".$f->hseqname."\t".$f->hstart."\t".$f->hend."\t".$f->hstrand."\n";
#    }
#    return 1;
    my @features;
    my %contig_hash;

    if ($out) {
	@features = @$out;
	print STDERR "Found out var, with ".scalar(@features)." features\n";
    }
    else {
	@features = @{$self->output};
    }
    
    foreach my $f (@features) {
      my $c = $f->seqname();
      my $targetid = $f->hseqname;
      $f->hseqname("mouse_$targetid");

      my $contig = $contig_hash{$c};
      unless($contig) {
	#cache contigs that have been obtained for speed purposes
	eval {
	  $contig = 
	    $self->db->get_RawContigAdaptor->fetch_by_name($c);
	  $contig_hash{$c} = $contig;
	};
	    
	if ($@) {
	  print STDERR "Contig not found, skipping writing output to db: $@\n";
	  return;
	}

	#attach contigs to features
	$f->attach_seq($contig);
      }
    }
    
    my $feat_adp = $self->db->get_DnaAlignFeatureAdaptor;
    $feat_adp->store(@features);
    
    return 1;
}

=head2 querydb

 Title   : querydb
 Usage   : $obj->querydb($newval)
 Function: Getset for querydb value
 Returns : value of querydb
 Args    : newvalue (optional)


=cut

sub querydb{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'querydb'} = $value;
    }
    return $obj->{'querydb'};

}

=head2 targetdb

 Title   : targetdb
 Usage   : $obj->targetdb($newval)
 Function: Getset for targetdb value
 Returns : value of targetdb
 Args    : newvalue (optional)


=cut

sub targetdb{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'targetdb'} = $value;
    }
    return $obj->{'targetdb'};

}

=head2 exonerate

 Title   : exonerate
 Usage   : $obj->exonerate($exonerate)
 Function: get/set for exonerate
 Returns : path to exonerate
 Args    : exonerate (optional)


=cut

sub exonerate {
   my ($self, $exonerate) = @_;

   if (!defined $self->{'_exonerate'}){
      $self->{'_exonerate'} = "";
   }

   if( defined $exonerate ) {
      $self->{'_exonerate'} = $exonerate;
    }
    return $self->{'_exonerate'};

}


=head2 exonerate_args

 Title   : exonerate_args
 Usage   : $obj->exonerate_args($exargs)
 Function: get/set for arguments to exonerate
 Returns : value of exonerate_args
 Args    : exargs (optional)


=cut

sub exonerate_args {
   my ($self, $exargs) = @_;

   if (!defined $self->{'_exonerate_args'}){
      $self->{'_exonerate_args'} = "";
   }

   if( defined $exargs ) {
      $self->{'_exonerate_args'} = $exargs;
    }
    return $self->{'_exonerate_args'};

}

=head2 max_query

 Title   : max_query
 Usage   : $obj->max_query($newval)
 Function: Getset for max_query value
 Returns : value of max_query
 Args    : newvalue (optional)


=cut

sub max_query{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'max_query'} = $value;
    }
    return $obj->{'max_query'};

}

=head2 max_target

 Title   : max_target
 Usage   : $obj->max_target($newval)
 Function: Getset for max_target value
 Returns : value of max_target
 Args    : newvalue (optional)


=cut

sub max_target{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'max_target'} = $value;
    }
    return $obj->{'max_target'};

}



1;
