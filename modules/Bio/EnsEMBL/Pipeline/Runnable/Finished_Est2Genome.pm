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

Bio::EnsEMBL::Pipeline::Runnable::Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new(
                                             -genomic => $genseq,
                                             -est     => $estseq 
                                             );
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Est2Genome->new()

=head1 DESCRIPTION

Object to store the details of an est2genome run.
Stores the est2genome matches as an array of Bio::EnsEMBL::FeaturePair

=head2 Methods:

 new,
 genomic_sequence,
 est_sequence,
 run,
 output.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::Runnable::Est2Genome;

#use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Est2Genome);


###########
# Analysis methods
##########

=head2 run

  Title   : run
  Usage   : $self->run()
            or
            $self->run("genomic.seq", "est.seq")
  Function: Runs est2genome and stores results as FeaturePairs
  Returns : TRUE on success, FALSE on failure.
  Args    : Temporary filenames for genomic and est sequences

=cut

# this is a MESS
sub run {
    my ($self, @args) = @_;
    
    # some constant strings
    my $source_tag  = "est2genome";
    my $dirname     = "/tmp/";
    
    #flag for est strand orientation
    my $estOrientation; 
    
    #check inputs
    my $genomicseq = $self->genomic_sequence ||
        $self->throw("Genomic sequence not provided");
    my $estseq = $self->est_sequence ||
        $self->throw("EST sequence not provided");
    #print "running est2genome\n";
    #extract filenames from args and check/create files and directory
    my ($genname, $estname) = $self->_rearrange(['genomic', 'est'], @args);
    my ($genfile, $estfile) = $self->_createfiles($genname, $estname, $dirname);
    #use appropriate Bio::Seq method to write fasta format files
    {
        my $genOutput = Bio::SeqIO->new(-file => ">$genfile" , '-format' => 'Fasta')
                    or $self->throw("Can't create new Bio::SeqIO from $genfile '$' : $!");
        my $estOutput = Bio::SeqIO->new(-file => ">$estfile" , '-format' => 'Fasta')
                    or $self->throw("Can't create new Bio::SeqIO from $estfile '$' : $!");

        #fill inputs
        $genOutput->write_seq($self->{'_genomic_sequence'}); 
        $estOutput->write_seq($self->{'_est_sequence'});

    }
    
             
    my $est_genome_command = "est2genome -reverse -genome $genfile -est $estfile -space 500000 -out stdout |";
        
    #print STDERR "running for " . $estseq->display_id . "\n";
    eval {
      #print (STDERR "Running command $est_genome_command\n");
      open (ESTGENOME, $est_genome_command) 
	or $self->throw("Can't open pipe from '$est_genome_command' : $!");
      
      #Use the first line to get gene orientation
      my $firstline = <ESTGENOME>;
      #print STDERR "firstline: \t$firstline\n";
      # put the gene on the minus strand iff splice sites imply reversed gene
      if ($firstline =~ /REVERSE/) { 
	#print STDERR "***reversed gene***\n"; 
	$estOrientation = -1; 
      }
      else {$estOrientation = 1}
     
      #read output
      while (<ESTGENOME>) {
	if ($_ =~ /Segmentation fault/) {
	  $self->warn("Segmentation fault from est_genome\n");
	  close (ESTGENOME) or $self->warn("problem running est_genome: $!\n");
	  return(0);
	}
	elsif ($_ =~ /^Segment/) {  # We only care about Segments
	  
            # "gen" = genomic sequence
            my ($primary, $score, $percent_id,
                $gen_start, $gen_end, $gen_id,
                $est_start, $est_end, $est_id) = (split)[0,1,2,
                                                         3,4,5,
                                                         6,7,8];
          
            ### Skip puny little bits ###
            next unless $score > 6;

            # Switch the starts and ends if we have a reverse strand gene
            if ($estOrientation == -1) {
                  ($gen_start, $gen_end,   $est_start, $est_end)
                = ($gen_end,   $gen_start, $est_end,   $est_start);
            }

            # Which  strand is the genomic match on?
            my $gen_strand = 1;
            if ($gen_start > $gen_end) {
                ($gen_start, $gen_end) = ($gen_end, $gen_start);
                $gen_strand = -1;
            }
            
            # Which strand is the est match on?
            my $est_strand = 1;
            if ($est_start > $est_end) {
                ($est_start, $est_end) = ($est_end, $est_start);
                $est_strand = -1;
            }
            
            # For Segements, make sure that the est_match
            # is on the forward strand, because est_strand is not
            # stored in the database and we can't view the
            # alignment for est_strand == -1 matches.
            if ($est_strand == -1) {
                $est_strand = 1;
                $gen_strand *= -1;
            }

            #print STDERR "genomic_strand = $gen_strand; est_strand = $est_strand\n";

 	    my $fp = $self->_create_FeaturePair(
                $score, $percent_id,
		$gen_start, $gen_end, $gen_id, 
		$est_start, $est_end, $est_id,
		$source_tag, 
		$gen_strand, $est_strand,   # est_strand is NOT stored in the db!
		$primary);
            $self->add_output($fp);
       }    

      }
      if(!close(ESTGENOME)){
	$self->warn("problem running est_genome: exit $?\n");
	return(0);
      }
    };

    $self->_deletefiles($genfile, $estfile);

    if ($@) {
        $self->throw("Error running est_genome:\n$@");
    } else {
        return 1;
    }
}





sub _create_FeaturePair {
    my ($self, $f1score, $f1percent_id, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id, $f1source, $f1strand, $f2strand, $f1primary) = @_;
    
    #print "creating feature pair ".$f1primary." ".$f1source." \n";
    my $analysis_obj    = new Bio::EnsEMBL::Analysis
                                (-db              => "none",
                                 -db_version      => "none",
                                 -program         => "est_genome",
                                 -program_version => "none",
                                 -gff_source      => $f1source,
                                 -gff_feature     => $f1primary,);
    #create features
    my $feat1 = new Bio::EnsEMBL::SeqFeature  (-start      =>   $f1start,
                                              -end         =>   $f1end,
                                              -seqname     =>   $f1id,
                                              -strand      =>   $f1strand,
                                              -score       =>   $f1score,
					      -percent_id  =>   $f1percent_id, 
					      -source_tag  =>   $f1source,
                                              -primary_tag =>   $f1primary,
                                              -analysis    =>   $analysis_obj );
     
   my $feat2 = new Bio::EnsEMBL::SeqFeature  (-start       =>   $f2start,
                                              -end         =>   $f2end,
                                              -seqname     =>   $f2id,
                                              -strand      =>   $f2strand,
                                              -score       =>   $f1score,
					      -percent_id  =>   $f1percent_id, 
                                              -source_tag  =>   $f1source,
                                              -primary_tag =>   $f1primary,
                                              -analysis    =>   $analysis_obj );

    #create featurepair
    my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $feat1,
                                             -feature2 => $feat2) ;
    return $fp;
}

sub output {
    my ($self) = @_;
    if (!defined($self->{'_output'})) {
	$self->{'_output'} = [];
    }
    
    return @{$self->{'_output'}};
}

sub add_output {
    my( $self, @feat_pairs ) = @_;
    
    push(@{$self->{'_output'}}, @feat_pairs);
}


1;
