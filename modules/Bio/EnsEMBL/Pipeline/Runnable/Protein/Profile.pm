
#
# Ensembl module for Profile
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Profile - DESCRIPTION of Object

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


package Bio::EnsEMBL::Pipeline::Runnable::Protein::Profile;
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
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::CPG->new (-CLONE => $seq);
    Function:   Initialises CPG object
    Returns :   a CPG Object
    Args    :   A Bio::Seq object (-CLONE), any arguments (-LENGTH, -GC, -OE) 

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    
    $self->{'_flist'}     = [];    # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef; # location of Bio::Seq object
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of Profile
    $self->{'_threshold'} = undef; # Value of the threshod
    $self->{'_protected'} = [];    # a list of files protected from deletion ???
    
  
     my( $clone, $analysis) = $self->_rearrange([qw(QUERY
						    ANALYSIS
						    )],
						@args);
    
    
    $self->clone ($clone) if ($clone);       
    $self->analysis ($analysis) if ($analysis);
        	
    return $self; # success - we hope!
}

######
#Get set methods
######

=head2 clone

 Title   : clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub clone{
    my ($self, $seq) = @_;
    if ($seq) {
	eval {
	    $seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")
	    };
	
	if (!$@) {
	    $self->{'_sequence'} = $seq ;
	    $self->queryname ($self->clone->id);
	    $self->filename ($self->clone->id.".$$.seq");
	    $self->results ($self->filename.".out");
	}
	else {
	    print STDERR "WARNING: The input_id is not a Seq object but if its a peptide fasta file, it should go fine\n";
	    $self->{'_sequence'} = $seq ;
	    $self->filename ("$$.tmp.seq");
	    
	    $self->results ("profile.$$.out");
	    
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

    # check clone
    my $seq = $self->clone || $self->throw("Clone required for Program");

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
	#The clone object is not a seq object but a file.
	#Perhaps should check here or before if this file is fasta format...if not die
	#Here the file does not need to be created or deleted. Its already written and may be used by other runnables.

	#pfscan runs over one protein at a time
	
	my $in  = Bio::SeqIO->new(-file => $seq, '-format' =>'Fasta');

	while ( my $tmpseq = $in->next_seq() ) {
	    print STDERR "SEQ: ".$self->workdir."/profile_tmp_seq\n";    
	    open (OUT,">".$self->workdir."/profile_tmp_seq");
	    print STDERR "AC: ".$tmpseq->display_id."\n";
	    #print STDERR "AC: ".$tmpseq->seq."\n";
	    
	    print OUT ">".$tmpseq->display_id."\n".$tmpseq->seq."\n";

	    close(OUT);
	$self->filename($self->workdir."/profile_tmp_seq");

	# run program
	$self->run_analysis;

	# parse output
	$self->parse_results($tmpseq->display_id);
	}
	$self->throw("Can't remove tmp sequence $!\n")
	    unless (system ("rm ".$self->workdir."/profile_tmp_seq") == 0);

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

    print STDERR "RUNNING: ".$self->analysis->program . ' -fz ' .$self->filename. ' ' .$self->analysis->db . ' > ' .$self->results."\n";

    $self->throw("Failed during Profile run $!\n")
	    
	unless (system ($self->analysis->program . ' -fz ' . 
			$self->filename. ' ' .
			$self->analysis->db . ' > ' .
			$self->results) == 0) ;
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
    my ($self,$sequenceId) = @_;
    
    my $filehandle;
    my $resfile = $self->results();
    
    if (-e $resfile) {
	
	if (-z $self->results) {  
	    print STDERR "pfscan didn't find any hits\n";
	    return; }       
	else {
	    open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
	    }
    }
    my %printsac;
    my $line;
    
    my @features;
    while (<CPGOUT>) {
	$line = $_;
	chomp $line;
	#print STDERR "$line\n";
	my ($nscore,$rawscore,$from,$to,$hfrom,$hto,$ac) = $line =~ /(\S+)\s+(\d+)\s*pos.\s+(\d*)\s*-\s+(\d*)\s*\[\s+(\d*),\s+(\S*)\]\s*(\w+)/;

	my $feat = "$ac,$from,$to,$hfrom,$hto,$nscore";
		
		push (@features,$feat);
	    }
		
    foreach my $feats (@features) {
	$self->create_feature($feats,$sequenceId);
	print STDERR "$feats\n";
    }
    @features = 0;
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
    
#my $feat = "$print,$start,$end,$percentageIdentity,$profileScore,$pvalue";
    my @f = split (/,/,$feat);
    
#$f[4] represents here the position of the match on the profile but if its values its -1, it means that the match is on the entire profile, thus we should calculate it.

    my $hto = $f[4];

    if ($f[4] =~ /-1/) {
#Calculate the lenght of the match using the values given for the sequence.
	$hto = $f[2] - $f[1] + 1;
    }


    my $feat1 = new Bio::EnsEMBL::SeqFeature ( -start => $f[1],                   
					       -end => $f[2],        
					       -score => $f[5],
					       -analysis => $self->analysis,
					       -seqname => $sequenceId,
					       -percent_id => 'NULL',
					       -p_value => 'NULL');
    
    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => 0,
					      -end => 0,
					      -analysis => $self->analysis,
					      -seqname => $f[0]);
    
    
    my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
						-feature2 => $feat2);
    
    if ($feature)
    {
	#$feat1->validate();
	
	# add to _flist
	push(@{$self->{'_flist'}}, $feature);
    }
}

