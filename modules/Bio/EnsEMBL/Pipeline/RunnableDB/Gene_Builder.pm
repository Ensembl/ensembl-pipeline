#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# written by Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder->new(
								    -db        => $db,
								    -input_id  => $id,
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

# Object preamble
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::GeneBuilder;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases   qw (
							       GB_FINALDBNAME
							       GB_FINALDBHOST
							       GB_FINALDBUSER
							       GB_FINALDBPASS
							       GB_FINALDBPORT
							       GB_COMB_DBHOST
							       GB_COMB_DBNAME
							       GB_COMB_DBUSER
							       GB_COMB_DBPASS
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder qw (
							       GB_VCONTIG
							       GB_FINAL_GENETYPE
							      );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General     qw (
							       GB_INPUTID_REGEX
							      );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

############################################################

=head2 new

    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
			   -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis,
			   -VCONTIG     => 1,
			   );

                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Gene_Builder object
    Args    :   -dbobj:      A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
                -vcontig:    determines whether it is running on virtual contigs
                             or RawContigs
                -extend:     determines the extension of the virtual contig
                             note: not implemented yet!
                -golden_path: determines the name of the golden path to use
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
    my( $use_vcontig) = $self->_rearrange([qw(VCONTIG)], @args);
       
    if (! defined $use_vcontig) {
	$use_vcontig = $GB_VCONTIG;
    }  
    
    $self->use_vcontig($use_vcontig);

    return $self;
}

############################################################

sub input_id {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_input_id} = $arg;
    }
    
    return $self->{_input_id};
}

############################################################

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut
    
    
sub write_output {
    my($self,@genes) = @_;
    
    # write genes out to a different database from the one we read genewise genes from.
    my $dbname = $GB_FINALDBNAME;
    my $dbhost = $GB_FINALDBHOST;
    my $dbuser = $GB_FINALDBUSER;
    my $dbpass = $GB_FINALDBPASS;
    my $dbport = $GB_FINALDBPORT;

    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						'-host'   => $dbhost,
						'-user'   => $dbuser,
						'-dbname' => $dbname,
						'-pass'   => $dbpass,
						'-port'   => $dbport,
						'-dnadb'  => $self->db,
					       );
    # sort out analysis
    my $genetype = $GB_FINAL_GENETYPE;
    unless ( $genetype ){
	$self->throw("Please, define GB_FINAL_GENETYPE in Config::GeneBuild::GeneBuilder");
    }
    my $analysis = $self->analysis;
    unless ($analysis){
	$self->throw("an analysis logic name must be defined in the command line");
    }
    
    my %contighash;
    my $gene_adaptor = $db->get_GeneAdaptor;
    
    # this now assummes that we are building on a single VC.
    my $genebuilders = $self->get_genebuilders;
    
    foreach my $contig ( keys %$genebuilders ){
        my $vc = $genebuilders->{$contig}->query;
	
	@genes = $genebuilders->{$contig}->final_genes;
	
	return unless ($#genes >= 0);
	my @newgenes;
	
	foreach my $gene (@genes) { 
	    $gene->analysis($analysis);
	    $gene->type($genetype);
	    
	    # coordinate transformation
	    if ($self->use_vcontig) {
		eval {
		    $gene->transform;
		};
		if ($@) {
		    $self->warn("Cannot convert gene to raw contigs:\n$@");
		    foreach my $tran (@{$gene->get_all_Transcripts}){
		      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($tran);
		    }
		}
	    }
	    
	    # store
	    eval {
		$gene_adaptor->store($gene);
		print STDERR "wrote gene " . $gene->dbID . " \n";
	    }; 
	    if( $@ ) {
		$self->warn("NABLE TO WRITE GENE:\n$@");
		foreach my $tran (@{$gene->get_all_Trascripts}){
		  Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($tran);
		}
	    }
	}
	
    }
    
}

############################################################

=head2 fetch_input

    Function:   It fetches the slice or contig according to the input_id, 
                and it defines the database where the
                previous annotations are stored and create a Bio::EnsEMBL::Pipeline::GeneBuilder
                object for that genomic, input_id and db
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);
    
    my $contigid  = $self->input_id;
    my $slice;
    
    if ($self->use_vcontig) {
	my $slice_adaptor = $self->db->get_SliceAdaptor();
	
	$contigid =~/$GB_INPUTID_REGEX/;
	
	my $chr   = $1;
	my $start = $2;
	my $end   = $3;
	print STDERR "Chr $chr - $start : $end\n";
	$slice   = $slice_adaptor->fetch_by_chr_start_end($chr,$start,$end);
    }
    else {
	# not sure this is the correct call:
	$slice = $self->db->get_Contig($contigid);
    }
    # database where all the genewise and combined genes are:
    my $genes_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						      '-host'   => $GB_COMB_DBHOST,
						      '-user'   => $GB_COMB_DBUSER,
						      '-dbname' => $GB_COMB_DBNAME,
						      '-pass'   => $GB_COMB_DBPASS,
						      '-dnadb'  => $self->db,
						      );

    print STDERR "reading genewise and combined genes from $GB_COMB_DBNAME : $GB_COMB_DBHOST\n";
    
    my $genebuilder = new Bio::EnsEMBL::Pipeline::GeneBuilder(
							      '-slice'   => $slice,
							      '-input_id' => $self->input_id,
							      );
    $genebuilder->genes_db($genes_db);
 
    # store the object and the piece of genomic where it will run
    $self->addgenebuilder($genebuilder,$slice);
    
}

############################################################

sub use_vcontig {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_vcontig} = $arg;
    }

    return $self->{_vcontig};
}

############################################################

sub addgenebuilder {
    my ($self,$arg,$contig) = @_;
    
    if (defined($arg) && defined($contig)) {
	$self->{_genebuilder}{$contig->id} = $arg;
    } 
    else {
	$self->throw("Wrong number of inputs [$arg,$contig]\n");
    }
}

############################################################

sub get_genebuilders {
    my ($self) = @_;
    
    return $self->{_genebuilder};
}

############################################################
	
sub run {
    my ($self) = @_;
    
    # get a hash, with keys = contig/slice and value = genebuilder object
    my $genebuilders = $self->get_genebuilders;
    
    my @genes;
    foreach my $contig (keys %{ $genebuilders } ) {
	my $query = $genebuilders->{$contig}->query;
	
	print(STDERR "GeneBuilding for $contig\n");
	
	$genebuilders->{$contig}->build_Genes;
	
	@genes = $genebuilders->{$contig}->final_genes;
    }
    
    $self->output( @genes );
}

############################################################

# override the evil RunnableDB output method:

sub output{
    my ($self, @genes ) = @_;
    unless ( $self->{_output} ){
	$self->{_output} = [];
    }
    if (@genes){
	push( @{$self->{_output}}, @genes );
    }
    return @{$self->{_output}};
}

############################################################



1;
