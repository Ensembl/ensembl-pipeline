
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
    $self->{'_program'}   = undef; # location of Profile executable
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_database'}  = undef; # Location of the database
    $self->{'_results'}   = undef; # file to store results of Profile
    $self->{'_threshold'} = undef; # Value of the threshod
    $self->{'_protected'} = [];    # a list of files protected from deletion ???
    
  
    print STDERR "args: ", @args, "\n";
    
    my( $query, $program, $database, $threshold, $workdir) = $self->_rearrange([qw(QUERY
										   PROGRAM
										   DATABASE
										   THRESHOLD
										   WORKDIR
										   )],
									       @args);
    
    #$self->clone($sequence) if ($sequence);       
  
    if ($query) {
	$self->clone($query);
    } else {
	$self->throw("No query sequence given");
    }
    
    if ($program) {   
	$self->program($program); }
    else {   
	$self->program($self->locate_executable('pfscan')); }
    
    if ($threshold) {
	$self->threshold($threshold);
    }
    
    if ($database) {
	$self->database($database);
    } else {
	$self->throw("No database given");
    }
    if ($workdir) {
	$self->workdir($workdir);
    }
	
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
    my ($self,$seq) = @_;
    
    if ($seq)
    {
	unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
	{
	    $self->throw("Input isn't a Bio::SeqI or Bio::PrimarySeqI");
	}
	$self->{'_sequence'} = $seq ;
	
	$self->filename($self->clone->id.".$$.seq");
	$self->results($self->filename.".out");
    }
    return $self->{'_sequence'};
}

=head2 program

 Title   : program
 Usage   : $obj->program($newval)
 Function: 
 Returns : value of program
 Args    : newvalue (optional)


=cut

sub program{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_program'} = $value;
    }
    return $obj->{'_program'};
}

=head2 database

 Title   : database
 Usage   : $obj->database($newval)
 Function: 
 Returns : value of database
 Args    : newvalue (optional)


=cut

sub database{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_database'} = $value;
    }
    return $obj->{'_database'};

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

    my $seq = $self->clone || $self->throw("Query seq required for Blast\n");

    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();

    #write sequence to file
    $self->writefile(); 
    $self->run_analysis();

    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
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

    print STDERR "RUNNING: ".$self->program . ' -fz ' .$self->filename. ' ' .$self->database . ' > ' .$self->results."\n";

    $self->throw("Failed during Profile run $!\n")
	    
	unless (system ($self->program . ' -fz ' . 
			$self->filename. ' ' .
			$self->database . ' > ' .
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
    my ($self) = @_;
    
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
    
    my $sequenceId;
    my @features;
    while (<CPGOUT>) {
	$line = $_;
	chomp $line;
	#print STDERR "$line\n";
	my ($nscore,$rawscore,$from,$to,$hfrom,$hto,$ac) = $line =~ /(\S+)\s+(\d+)\s*pos.\s+(\d*)\s*-\s+(\d*)\s*\[\s+(\d*),\s+(\S*)\]\s*(\w+)/;
	#print STDERR "$nscore,$rawscore,$from,$to,$hfrom,$hto,$ac\n";


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
    
    #create analysis object
    my $analysis_obj = Bio::EnsEMBL::Analysis->new
	(   -db              => "PROSITE",
	    -db_version      => 1,
	    -program         => "pfscan",
	    -program_version => 1,
	    -gff_source      => "profile",
	    -gff_feature     => "domain");
    
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
					       -analysis => $analysis_obj,
					       -seqname => $self->clone->id,
					       -percent_id => 'NULL',
					       -p_value => 'NULL');
    
    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $f[3],
					      -end => $hto,
					      -analysis => $analysis_obj,
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

