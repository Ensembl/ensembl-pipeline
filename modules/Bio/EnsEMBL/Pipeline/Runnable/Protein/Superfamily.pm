
#
# Ensembl module for Superfamily
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Superfamily - DESCRIPTION of Object

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


package Bio::EnsEMBL::Pipeline::Runnable::Protein::Superfamily;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);



=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::CPG->new (-QUERY => $seq);
    Function:   Initialises CPG object
    Returns :   a CPG Object
    Args    :   A Bio::Seq object (-QUERY), any arguments (-LENGTH, -GC, -OE) 

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    
    $self->{'_flist'}     = [];    # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef; # location of Bio::Seq object
    $self->{'_program'}   = undef; # location of Superfamily executable
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_database'}  = undef; # Location of the database
    $self->{'_results'}   = undef; # file to store results of Superfamily
    $self->{'_threshold'} = undef; # Value of the threshod
    $self->{'_protected'} = [];    # a list of files protected from deletion ???
    



    my ($query, $analysis) = $self->_rearrange([qw(QUERY 
						   ANALYSIS)], 
					       @args);

    print STDERR "ANALYSIS: $analysis\n";
    
    $self->query ($query) if ($query);       
    $self->analysis ($analysis) if ($analysis);
	
    return $self; # success - we hope!
}

######
#Get set methods
######

=head2 query

 Title   : query
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub query{
      my ($self, $seq) = @_;
    if ($seq) {
	eval {
	    $seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")
	};

	if (!$@) {
	    $self->{'_sequence'} = $seq ;
	    $self->queryname ($self->query->id);
	    $self->filename ($self->query->id.".$$.seq");
	    $self->results ($self->filename.".out");
	}
	else {
	    print STDERR "WARNING: The input_id is not a Seq object but if its a peptide fasta file, it should go fine\n";
	    $self->{'_sequence'} = $seq ;
	    $self->filename ("$$.tmp.seq");
	    
	    $self->results ("Superfamily.$$.out");
	    
	}
    }
    return $self->{'_sequence'};
}


=head2 analysis

 Title   : analysis
 Usage   : $obj->analysis($newval)
 Function: 
 Returns : value of analysis
 Args    : newvalue (optional)


=cut

sub analysis{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'analysis'} = $value;
    }
    return $obj->{'analysis'};

}


###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
 my ($self, $dir) = @_;

    # check query
    my $seq = $self->query || $self->throw("Query required for Program");

    # set directory if provided
    $self->workdir ('/tmp') unless ($self->workdir($dir));
    $self->checkdir;

    # reset filename and results as necessary (adding the directory path)
    my $tmp = $self->workdir;
    my $input = $tmp."/".$self->filename;
    $self->filename ($input);
    $tmp .= "/".$self->results;
    $self->results ($tmp);


 eval {
	$seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")
	};
	

    if (!$@) {
	#The inputId is a sequence file...got the normal way...

	# write sequence to file
	$self->writefile;        

	# run program
	$self->run_analysis;

	# parse output
	$self->parse_results;
	$self->deletefiles;
    }
    else {
	#The query object is not a seq object but a file.
	#Perhaps should check here or before if this file is fasta format...if not die
	#Here the file does not need to be created or deleted. Its already written and may be used by other runnables.

	$self->filename($self->query);

	# run program
	$self->run_analysis;

	# parse output
	$self->parse_results;
    }
  }

=head2 run_analysis

    Title   :   run_analysis
    Usage   :   $obj->run_analysis
    Function:   Runs the blast query
    Returns :   nothing
    Args    :   none

=cut

sub run_analysis {
    my ($self) = @_;

    print STDERR "RUNNING: ".$self->analysis->program." ".$self->filename." ".$self->results." /tmp ".$self->analysis->db_file." /acari/work1/mongin/superfamily/model.tab /usr/local/ensembl/sam/bin/hmmscore 0.02 y /acari/work1/mongin/superfamily"."\n";
	$self->throw("Failed during Superfamily run $!\n")
	    
	    unless (system ($self->analysis->program . ' ' .
			    $self->filename.' '.
			    $self->results.
			    ' /tmp '.
                            $self->analysis->db_file.
			    ' /acari/work1/mongin/superfamily/model.tab /usr/local/ensembl/sam/bin/hmmscore 0.02 y /acari/work1/mongin/superfamily') == 0);
}



=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses cpg output to give a set of features
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut
sub parse_results {
    my ($self) = @_;
   
    my $filehandle;
    my $resfile = $self->results();
    
    if (-e $resfile) {
	
	if (-z $self->results) {  
	    print STDERR "Superfamilyscan didn't find any hits\n";
	    return; }       
	else {
	    open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
	    }
    }
    my %superfamilyac;
    my $line;
    
    my @features;
    my %feat;
    while (<CPGOUT>) {
	chomp;
	my ($seqid,$a,$familyid,$start,$end,$evalue);
	my $tmp = "$a:$familyid:$start:$end:$evalue";

	print STDERR "$tmp\n";
	push(@{$feat{$seqid}},$tmp);
    }
	

    foreach my $id ( keys %feat ) {
	my @array = @{$feat{$id}};
	foreach my $tmp (@array) {
	    $self->create_feature($tmp,$id);
	}
    }
}
	


##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of features
    Returns :   Returns an array of features
    Args    :   none

=cut

sub output {
    my ($self) = @_;
    return @{$self->{'_flist'}};
}

=head2 create_feature

    Title   :   create_feature
    Usage   :   obj->create_feature($feature)
    Function:   Returns an array of features
    Returns :   Returns an array of features
    Args    :   none

=cut
sub create_feature {
    my ($self, $feat, $sequenceId) = @_;

    my @f = split (/:/,$feat);
    
    
    my $feat1 = new Bio::EnsEMBL::SeqFeature ( -start => $f[2],                   
					       -end => $f[3],        
					       -analysis => $self->analysis,
					       -seqname => $sequenceId,
					       -p_value => $f[4]);
    
    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start =>0,
					      -end => 0,
					      -analysis => $self->analysis,
					      -seqname => $f[1]);
    
    
    my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
						-feature2 => $feat2);
    
    if ($feature)
    {
	#$feat1->validate();
	
	# add to _flist
	push(@{$self->{'_flist'}}, $feature);
    }
}

