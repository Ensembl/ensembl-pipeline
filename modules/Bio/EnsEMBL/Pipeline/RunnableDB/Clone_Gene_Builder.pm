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

Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Gene_Builder

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Gene_Builder->new(
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

package Bio::EnsEMBL::Pipeline::RunnableDB::Clone_Gene_Builder;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::GeneBuilder;
#use Bio::EnsEMBL::DB::ConvertibleVirtualContig;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Utils::GTF_handler;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_GOLDEN_PATH
					 );
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
			   -SEQFETCHER  => $sf,
			   -GOLDEN_PATH => $path,
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Pipeline::Analysis (optional) 
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $vcontig,$extend, $path ) = $self->_rearrange([qw(VCONTIG EXTEND GOLDEN_PATH)], @args);
       
    $vcontig = 1 unless defined($vcontig);
    
    $self->vcontig($vcontig);
    $self->extend($extend);

    # golden path
    if(!defined $path){
      $path = $GB_GOLDEN_PATH;
    }

    $path = 'UCSC' unless (defined $path && $path ne '');
    $self->dbobj->static_golden_path_type($path);


    return $self;
}

sub input_id {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_input_id} = $arg;
    }
    
    return $self->{_input_id};
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

    # not ideal, but temporarily write to the same database
    my $db = $self->dbobj;

    # sort out analysis
    my $genetype = 'ensembl';
    my $anaAdaptor = $db->get_AnalysisAdaptor;
    my @analyses = $anaAdaptor->fetch_by_logic_name($genetype);
    
    my $analysis_obj;
    
    if(scalar(@analyses) > 1){
      $self->throw("panic! > 1 analysis for $genetype\n");
    }
    elsif(scalar(@analyses) == 1){
      $analysis_obj = $analyses[0];
    }
    else{
      # make a new analysis object
      $analysis_obj = new Bio::EnsEMBL::Analysis
	(-db              => 'NULL',
	 -db_version      => 1,
	 -program         => $genetype,
	 -program_version => 1,
	 -gff_source      => $genetype,
	 -gff_feature     => 'gene',
	 -logic_name      => $genetype,
	 -module          => 'GeneBuilder',
	);
    }
    
    my %contighash;
    my $gene_adaptor = $db->get_GeneAdaptor;

    # this now assummes that we are building on a single VC.
    my $genebuilders = $self->get_genebuilders;
    my ($contig)     = keys %$genebuilders;
    my $vc = $genebuilders->{$contig}->contig;

    @genes = $genebuilders->{$contig}->each_Gene;

    return unless ($#genes >= 0);
    my @newgenes;

    foreach my $gene (@genes) { 
      eval {
	print STDERR "Converting gene " . $gene->{'temporary_id'} . "\n";
	$gene->analysis($analysis_obj);
	my $newgene = $vc->convert_Gene_to_raw_contig($gene);

	push(@newgenes,$newgene);
	$newgene->type($genetype);
	$newgene->analysis($analysis_obj);
      };
      if ($@) {
	print STDERR "ERROR converting gene " . $gene->{'temporary_id'} . " " . $gene->type . " \n[$@]\n";
      }
    }

  GENE: foreach my $gene (@newgenes) {	
      # do a per gene eval...
      eval {
	$gene_adaptor->store($gene);
	print STDERR "wrote gene " . $gene->dbID . " \n";
      }; 
      if( $@ ) {
	print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
      }
    
    }
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for genebuilder from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;

    $self->throw("No input id") unless defined($self->input_id);

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();

    my $fpcid  = $self->input_id;
    
    print STDERR "FPC contig $fpcid\n";
    my $contig    = $stadaptor->fetch_VirtualContig_by_fpc_name($fpcid);

    $contig->primary_seq;

    print STDERR "Length of primary seq is ",$contig->primary_seq->length,"\n";

    my $genebuilder = new Bio::EnsEMBL::Pipeline::GeneBuilder(-contig => $contig,
							      -input_id => $self->input_id,
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
    #my @gene;

    $self->{'_output'} = [];
    
    my @vcgenes;
    foreach my $contig (keys %$genebuilders) {
        my $vc = $genebuilders->{$contig}->contig;
	print(STDERR "Building for $contig\n");
	$genebuilders->{$contig}->build_Genes;
	@vcgenes = @{$genebuilders->{$contig}{_genes}};
        print STDERR "Genes before conversion\n";
        $genebuilders->{$contig}->print_Genes(@vcgenes);
        print STDERR "Converting coordinates\n";
    }
    
	    
    push(@{$self->{'_output'}},@vcgenes);
}

1;
