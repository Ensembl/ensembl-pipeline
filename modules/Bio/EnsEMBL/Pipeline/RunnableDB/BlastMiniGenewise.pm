#!/usr/local/bin/perl

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

Bio::EnsEMBL::Pipeline::RunnableDB::AlignFeature

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::MiniGenewise->new(
					     -dbobj     => $db,
					     -input_id  => $id
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastMiniGenewise;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    my( $dbobj,$input_id ) = $self->_rearrange(['DBOBJ',
						'INPUT_ID'], @args);
       
    $self->throw("No database handle input")                 unless defined($dbobj);
    $self->dbobj($dbobj);

    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);
    
    return $self; # success - we hope!
}
sub input_id {
	my ($self,$arg) = @_;

   if (defined($arg)) {
      $self->{_input_id} = $arg;
   }

   return $self->{_input_id};
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($db)
    Function:   Get/set method for database handle
    Returns :   Bio::EnsEMBL::Pipeline::DB::ObjI
    Args    :   

=cut

sub dbobj {
    my( $self, $value ) = @_;    
    if ($value) {

        $value->isa("Bio::EnsEMBL::DB::ObjI") || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DB::ObjI");
        $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
}

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output($file_name);
    Function:   Fetchs output data from a frozen perl object
                stored in file $file_name
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub fetch_output {
    my($self,$output) = @_;
    
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self,@features) = @_;

    my $db=$self->dbobj();

    @features || $self->throw("No frozen object passed for the output");
    my $contig;
    eval {
	$contig = $db->get_Contig($self->input_id);
    };
    if ($@) {
	print STDERR "Contig [" . $self->input_id . "] not found, skipping writing output to db ($@\n";
    }
    else {
	my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	$feat_Obj->write($contig,@features);
    }
    return 1;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for est2genome from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    print STDERR "Fetching input \n";
    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);
    my $genseq    = $contig->primary_seq;
    my @features  = $contig->get_all_SimilarityFeatures_above_score('blastp',1);

    my %idhash;

    foreach my $f (@features) {
	if (defined($f->hseqname)) {
	    $idhash{$f->hseqname} = 1;
	}
    }

    my @ids = keys %idhash;
    
    print STDERR "Feature ids are @ids\n";

    $self->{_genseq} = $genseq;
    $self->{_ids}    = \@ids;
}
     
sub runnable {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	if ($arg->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
	    $self->{_runnable} = $arg;
	} else {
	    $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI");
	}
    }
    return $self->{_runnable};
}

sub run {
    my ($self) = @_;

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise('-genomic'  => $self->{_genseq},
									   '-ids'      => $self->{_ids});
								     
    
    $self->runnable($runnable);
    $self->runnable->run;
    $self->convert_output;

}

sub convert_output {
    my ($self) =@_;
    my @tmpf = $self->runnable->output;
    my @genes;
    my $count = 1;
    my $time  = time; chomp($time);

    foreach my $tmpf (@tmpf) {
	my $gene   = new Bio::EnsEMBL::Gene;
	my $tran   = new Bio::EnsEMBL::Transcript;
	my $transl = new Bio::EnsEMBL::Translation;

	$gene->id($self->input_id . ".genewise.$count");
	$gene->created($time);
	$gene->modified($time);
	$gene->version(1);

	$tran->id($self->input_id . ".genewise.$count");
	$tran->created($time);
	$tran->modified($time);
	$tran->version(1);

	$transl->id($self->input_id . ".genewise.$count");
	$transl->version(1);



	$gene->add_Transcript($tran);
	$tran->translation($transl);

	push(@genes,$gene);

	my $excount = 1;
	my @exons;

	foreach my $subf ($tmpf->sub_SeqFeature) {
	    my $exon = new Bio::EnsEMBL::Exon;
	    $exon->id($self->input_id . ".genewise.$count.$excount");
	    $exon->created($time);
	    $exon->modified($time);
	    $exon->version(1);

	    $exon->start($subf->start);
	    $exon->end  ($subf->end);
	    $exon->strand($subf->strand);
	    
	    # This is dummy phase
	    $exon->phase(0);
	    push(@exons,$exon);

	    $excount++;
	}
	$count++;
	if ($exons[0]->strand == -1) {
	    @exons = sort {$b->start <=> $a->start} @exons;
	} else {
	    @exons = sort {$a->start <=> $b->start} @exons;
	}

	foreach my $exon (@exons) {
	    $tran->add_Exon($exon);
	}

	$transl->start_exon_id($exons[0]->id);
	$transl->end_exon_id  ($exons[$#exons]->id);

	my @subf = $tmpf->sub_SeqFeature;

	if ($#subf > 0) {
	    my $i;
	    for ($i = 0; $i <= $#subf-1; $i++) {
		$self->check_splice($subf[$i],$subf[$i+1]);
	    }
	}
    }

    $self->{_output} = \@genes;

}

sub check_splice {
    my ($self,$f1,$f2) = @_;
    
    my $splice1 = substr($self->{_genseq}->seq,$f1->end,2);
    my $splice2 = substr($self->{_genseq}->seq,$f2->start-3,2);
    
    if (abs($f2->start - $f1->end) > 50) {
	print ("Splices are " . $f1->hseqname . " [" . 
	                        $splice1      . "][" . 
	                        $splice2      . "] " . 
	       ($f2->start - $f1->end)        . "\n");
    }
}


sub output {
    my ($self) = @_;

    if (defined($self->{_output})) {
	return @{$self->{_output}};
    } 
	
}


1;


