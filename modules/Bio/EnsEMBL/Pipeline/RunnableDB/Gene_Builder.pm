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

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome->new(
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

package Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableDBI;
use Bio::EnsEMBL::Pipeline::GeneBuilder;
use Bio::EnsEMBL::DB::ConvertibleVirtualContig;
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $dbobj,$input_id,$vcontig,$extend ) = $self->_rearrange([qw(DBOBJ INPUT_ID VCONTIG EXTEND)], @args);
       
    $self->throw("No database handle input")           unless defined($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::DB::ObjI");
    $self->dbobj($dbobj);

    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);

    $vcontig = 1 unless defined($vcontig);
    
    $self->vcontig($vcontig);
    $self->extend($extend);

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
    
    $output || $self->throw("No frozen object passed for the output");
    
    my $object;
    open (IN,"<$output") || do {print STDERR ("Could not open output data file... skipping job\n"); next;};
    
    while (<IN>) {
    $_ =~ s/\[//;
	$_ =~ s/\]//;
	$object .= $_;
    }
    close(IN);
    my @out;
   
    if (defined($object)) {
    my (@obj) = FreezeThaw::thaw($object);
    foreach my $array (@obj) {
	foreach my $object (@$array) {
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
    Returns :   array of exons (with start and end)
    Args    :   none

=cut


sub write_output {
    my($self,@genes) = @_;

    my $db = $self->dbobj();
   
    my %contighash;
    my $gene_obj = $db->gene_Obj;

        return unless ($#genes >= 0);


	      my $sth = $db->prepare("lock tables gene write, exon write, transcript write, exon_transcript write, translation write,dna read,contig read,clone read,feature read,analysis read");
	      $sth->execute;

	  foreach my $gene (@genes) {
	    (my $gcount = $gene_obj->get_new_GeneID($GENE_ID_SUBSCRIPT))
		=~ s/$GENE_ID_SUBSCRIPT//;
	    (my $tcount = $gene_obj->get_new_TranscriptID($TRANSCRIPT_ID_SUBSCRIPT))
		=~ s/$TRANSCRIPT_ID_SUBSCRIPT//;
	    (my $pcount = $gene_obj->get_new_TranslationID($PROTEIN_ID_SUBSCRIPT))
		=~ s/$PROTEIN_ID_SUBSCRIPT//;
	    (my $ecount = $gene_obj->get_new_ExonID($EXON_ID_SUBSCRIPT))
		=~ s/$EXON_ID_SUBSCRIPT//;
	    
	    $gene->id($GENE_ID_SUBSCRIPT . $gcount);

	    print (STDERR "Writing gene " . $gene->id . "\n");

            # Convert all exon ids and save in a hash
            my %namehash;

            foreach my $ex ($gene->each_unique_Exon) {
                  print STDERR "Exon id " . $ex . " " . $ex->id . " " . $EXON_ID_SUBSCRIPT . "\n";
                  $namehash{$ex->id} = $EXON_ID_SUBSCRIPT.$ecount;
                  $ex->id($EXON_ID_SUBSCRIPT.$ecount);
	   	  $ecount++;
            }

	    print (STDERR "Transcripts are\n");
	    foreach my $tran ($gene->each_Transcript) {
		$tran->id             ($TRANSCRIPT_ID_SUBSCRIPT . $tcount);
		$tran->translation->id($PROTEIN_ID_SUBSCRIPT . $pcount);

		my $translation = $tran->translation;

		$tcount++;
		$pcount++;

		print (STDERR "Transcript  " . $tran->id . "\n");
		print (STDERR "Translation " . $tran->translation->id . "\n");

		foreach my $ex ($tran->each_Exon) {
                    my @sf = $ex->each_Supporting_Feature;
                    print STDERR "Supporting features are " . scalar(@sf) . "\n";

                    if ($namehash{$translation->start_exon_id} ne "") {
		      $translation->start_exon_id($namehash{$translation->start_exon_id});
                    }
                    if ($namehash{$translation->end_exon_id} ne "") {
		      $translation->end_exon_id  ($namehash{$translation->end_exon_id});
                    }
		    print(STDERR "Exon         " . $ex->id . "\n");
		}
		
	    }

	    $gene_obj->write($gene);

	}
	      $sth = $db->prepare("unlock tables");
	      $sth->execute;

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

    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->dbobj->get_Contig($contigid);

    if ($self->vcontig) {
	my $focus = int(($contig->golden_start + $contig->golden_end)/2);

	$self->focuscontig($contig);

	$contig = new Bio::EnsEMBL::DB::ConvertibleVirtualContig(-focuscontig   => $contig,
						      -focusposition => $focus,
						      -ori           => 1,
						      -left          => $self->extend,
						      -right         => $self->extend);

	$contig = $contig->extend_maximally;

	print(STDERR "Contig length is " . $contig->length . "\n");

	if ($contig->length > 50000000) {
	    $self->throw("Aborting - virtual contig loo tong");
	}

    }

    $contig->_dump_map(\*STDOUT);
    $contig->primary_seq;

    print STDERR "Length of primary seq is",$contig->primary_seq->length,"\n";

    my $analysis = $self->dbobj->get_OldAnalysis(8);
           
    my $genebuilder = new Bio::EnsEMBL::Pipeline::GeneBuilder(-contig => $contig,
							      -analysis => $analysis,
							      );



    $self->addgenebuilder($genebuilder,$contig);
    
}
sub focuscontig {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_contig} = $arg;
    }

    return $self->{_contig};
}

sub vcontig {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_vcontig} = $arg;
    }

    return $self->{_vcontig};
}

sub extend {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_extend} = $arg;
    }

    return $self->{_extend} || 400000;
}

sub addgenebuilder {
    my ($self,$arg,$contig) = @_;

    if (defined($arg) && defined($contig)) {
	$self->{_genebuilder}{$contig->id} = $arg;
    } else {
	$self->throw("Wrong number of inputs [$arg,$contig]\n");
    }
}

sub get_genebuilders {
    my ($self) = @_;

    return $self->{_genebuilder};
}

sub check_gene {
   my ($self,$gene) = @_;

   foreach my $tran ($gene->each_Transcript) {
      my $seq = $tran->translate->seq;

      if ($seq =~ /\*/) {
        $self->throw("Stop codons in gene " . $gene->id . " transcript " . $tran->id . " - exiting");
      }
   }
}
	
sub run {
    my ($self) = @_;

    my $genebuilders = $self->get_genebuilders;
    my @gene;

    $self->{_output} = [];
    
    foreach my $contig (keys %$genebuilders) {
        my $vc = $genebuilders->{$contig}->contig;
	print(STDERR "Building for $contig\n");
	$genebuilders->{$contig}->build_Genes;
	my @vcgenes = @{$genebuilders->{$contig}{_genes}};
        print STDERR "Genes before conversion\n";
        $genebuilders->{$contig}->print_Genes(@vcgenes);
        print STDERR "Converting coordinates";
        foreach my $g (@vcgenes) {
           my $newgene = $vc->convert_Gene_to_raw_contig($g);
           #$self->check_gene($newgene);
	   push(@gene,$newgene);
        }
         $vc = undef;
         $genebuilders->{$contig} = undef;
    }
    
    foreach my $gene (@gene) {
	foreach my $exon ($gene->each_unique_Exon) {
	    $exon->{_gsf_seq} = undef;
	}
    }
	    
    push(@{$self->{_output}},@gene);
}

sub output {
    my ($self) = @_;

    if (!defined($self->{_output})) {
	$self->{_output} = [];
    }
    return @{$self->{_output}};
}


1;






