# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

 Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil

=head1 SYNOPSIS

 my $seqstream = Bio::SeqIO->new ( -file => $clonefile,
                                   -fmt => 'Fasta',
                                 );
 $seq = $seqstream->next_seq;

 my $ncoils = Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil->new ( -CLONE => $seq);
 $ncoils->workdir ($workdir);
 $ncoils->run;
 my @results = $ncoils->output;

=head1 DESCRIPTION

 Coil takes a Bio::Seq (or Bio::PrimarySeq) object
 and runs ncoils on it (detecting coiled coils). 
 The resulting output file is parsed to produce a set of features.

=head1 CONTACT

 Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods. 
 Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil;

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


=head2 new

 Title    : new
 Usage    : my $ncoils =  Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil->new
                          ( -program    => '/usr/local/pubseq/bin/ncoils',
                            -clone      => $clone,
                            -analysi    => $analysis,
                          );
 Function : initialises Coil object
 Returns  : a Coil object
 Args     : a Bio::Seq object, the path to the program binaries and an analysisId
 Throws   :

=cut

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new (@_);    
  
    $self->{'_flist'}     = [];           # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef;        # location of Bio::Seq object
    $self->{'_program'}   = undef;        # location of executable
    $self->{'_workdir'}   = undef;        # location of tmp directory
    $self->{'_filename'}  = undef;        # file to store Bio::Seq object
    $self->{'_results'}   = undef;        # file to store results of program run
    $self->{'_protected'} = [];           # a list of files protected from deletion
  
    my ($clone, $analysis, $program) = $self->_rearrange([qw(CLONE 
						             ANALYSIS
                                                             PROGRAM)], 
					                 @args);
  
    
    $self->clone ($clone) if ($clone);

    if ($analysis) {
        $self->analysis ($analysis);
    }
    else {
        $self->throw("Coil needs an analysis");
    }

    $self->program ($self->find_executable ($self->analysis->program_file));
  
    return $self;
}


###################
# get/set methods 
###################

=head2 clone

 Title    : clone
 Usage    : $self->clone ($clone);
 Function : get/set method for the Sequence object; assigns clone, filename and results
 Example  :
 Returns  : a Bio::Seq or Bio::PrimarySeq object
 Args     : a Bio::Seq or Bio::PrimarySeq object (optional)
 Throws   :

=cut

sub clone {
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
	    
	    $self->results ("coil.$$.out");
	    
	}
    }
    return $self->{'_sequence'};
}


=head2 analysis

 Title    : analysis
 Usage    : $self->analysis ($analysis);
 Function : get/set method for the analysis
 Example  :
 Returns  : analysis
 Args     : analysis (optional)
 Throws   :

=cut

sub analysis {
    my $self = shift;
    if (@_) {
        $self->{'_analysis'} = shift;
    }
    return $self->{'_analysis'};
} 


=head2 program

 Title    : program
 Usage    : $self->program ('/usr/local/pubseq/bin/ncoils');
 Function : get/set method for the path to the executable
 Example  :
 Returns  : File path
 Args     : File path (optional)
 Throws   :

=cut

sub program {
    my $self = shift;
    if (@_) {
        $self->{'_program'} = shift;
    }
    return $self->{'_program'};
}

###################
# analysis methods
###################

=head2 run

 Title    : run
 Usage    : $self->run ($workdir, $args)
 Function : runs program and populates @{$self->{'_flist'}} (array of features)
 Example  :
 Returns  :   
 Args     : workdir (optional)
 Throws   :

=cut

sub run {
    my ($self, $dir) = @_;

   
    # check clone
    my $seq = $self->clone || $self->throw("Clone required for Coil\n");

   
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
	$self->run_program;

	# parse output
	$self->parse_results;
	$self->deletefiles;
    }
    else {
	#The clone object is not a seq object but a file.
	#Perhaps should check here or before if this file is fasta format...if not die
	#Here the file does not need to be created or deleted. Its already written and may be used by other runnables.
	
	$self->filename($self->clone());
	
	# run program
	$self->run_program;

	# parse output
	$self->parse_results;
    }

}


=head2 run_program

 Title    : run_program
 Usage    : $self->program
 Function : makes the system call to program
 Example  :
 Returns  : 
 Args     :
 Throws   :

=cut

sub run_program {
    my ($self) = @_;

    my $coilsdir='/usr/local/ensembl/data/coils';
    $ENV{'COILSDIR'}=$coilsdir;
    
    # run program
    print STDERR "running ".$self->program."\n";
    $self->throw ("Error running ".$self->program." on ".$self->filename) 
        unless ((system ($self->program." -f < ".$self->filename." > ".$self->results)) == 0); 
}


=head2 parse_results

 Title    :  parse_results
 Usage    :  $self->parse_results ($filename)
 Function :  parses program output to give a set of features
 Example  :
 Returns  : 
 Args     : filename (optional, can be filename, filehandle or pipe, not implemented)
 Throws   :

=cut

sub parse_results {
    my ($self) = @_;
    my $filehandle;
    my $resfile = $self->results;
    
    if (-e $resfile) {
        # it's a filename
        if (-z $self->results) {  
	    print STDERR $self->program." didn't find anything\n";
	    return;
        }       
        else {
            open (OUT, "<$resfile") or $self->throw ("Error opening $resfile");
            $filehandle = \*OUT;
      }
    }
    else {
        # it'a a filehandle
        $filehandle = $resfile;
    }
    
    # parse
    my %result_hash = _read_fasta ($filehandle);
    close $filehandle;   

    foreach my $id (keys %result_hash) {
		
        my $pep = reverse ($result_hash{$id});
        my $count = my $switch = 0;
        my ($start, $end);
        while (my $aa = chop $pep) {
            $count++;
            if (!$switch && $aa eq "x") {
                $start = $count;
                $switch = 1;
   	    }
            elsif ($switch && $aa ne "x") {
                $end = $count-1;
                my (%feature);
	        $feature{name} = $id;
       	        $feature{start} = $start;
	        $feature{end} = $end;
    	        ($feature{source}) = $self->program =~ /([^\/]+)$/;
	        $feature{primary}= 'coiled_coil';
	        ($feature{program}) = $self->program =~ /([^\/]+)$/;
                $feature{logic_name} = 'coiled_coil';
  	        $self->create_feature (\%feature);
                $switch = 0;
	    }
        }
    }
}


=head2 create_feature

 Title    : create_feature
 Usage    : $self->create_feature ($feature)
 Function : creates a Bio::EnsEMBL::SeqFeature object from %feature,
            and pushes it onto @{$self->{'_flist'}}
 Example  :
 Returns  :
 Args     :
 Throws   :

=cut

sub create_feature {
    my ($self, $feat) = @_;

    my $analysis = $self->analysis;

    # create feature object
    my $feat1 = Bio::EnsEMBL::SeqFeature->new ( -seqname     => $feat->{name},
                                                  -start       => $feat->{start},
                                                  -end         => $feat->{end},
                                                  -source_tag  => $feat->{source},
                                                  -primary_tag => $feat->{primary},
                                                  -analysis    => $analysis,
						  -percent_id => 'NULL',
						  -p_value => 'NULL',
                                                ); 
  
    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => 0,
					      -end => 0,
					      -analysis => $analysis,
					      -seqname => 'ncoils');
    
    
    my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
						-feature2 => $feat2);
    
    if ($feature) {
	#$feature->validate_prot_feature;
	# add to _flist
	push (@{$self->{'_flist'}}, $feature);
    }
}


=head2 output

 Title    : output
 Usage    : $self->output
 Function : returns an array of feature objects
 Example  :
 Returns  : an array of Bio::EnsEMBL::SeqFeature objects
 Args     :
 Throws   :

=cut

sub output {
    my ($self) = @_;
    my @list = @{$self->{'_flist'}};
    return @{$self->{'_flist'}};
}

1;

#############################
# subroutines
#############################
sub _read_fasta {
    local (*FILE) = @_;
    my ($id , $seq , %name2seq);
    while (<FILE>) {
        chomp;
        if (/^>(\S+)/) {
            my $new_id = $1;
            if ($id) {
                $name2seq{$id} = $seq;
            }
            $id = $new_id ; $seq = "" ;
        } 
        elsif (eof) {
            if ($id) {
                $seq .= $_ ;
                $name2seq{$id} = $seq;
            }
        }
        else {
            $seq .= $_ ;
        }
    }
    return %name2seq;
}
