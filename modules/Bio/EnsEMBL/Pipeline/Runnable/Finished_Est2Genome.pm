#
#
# 
#
# Copyright GRL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome->new(
									 -genomic => $genseq,
									 -est     => $estseq 
									);
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome->new()

=head1 DESCRIPTION

needs to parse and output the data differenty to the standard est2genome so implements only methods where this has changed

=head2 Methods:

run
_createfeature
convert_output

=head1 CONTACT

ensembl-dev <ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::Finished_Est2Genome;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

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
    
   
    my $source_tag  = "est2genome";
    my $dirname     = "/tmp/";
    
    
    my $estOrientation; 
    
    #check inputs
    my $genomicseq = $self->genomic_sequence ||
        $self->throw("Genomic sequence not provided");
    my $estseq = $self->est_sequence ||
        $self->throw("EST sequence not provided");
    
    my ($genname, $estname) = $self->_rearrange(['genomic', 'est'], @args);
    my ($genfile, $estfile) = $self->_createfiles($genname, $estname, $dirname);
    my $output_file = $estfile.".output";    
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
        
    #The -reverse switch ensures correct numbering on EST seq in either orientation
    my $est_genome_command = "est_genome  -reverse -genome $genfile -est $estfile |";
    #my $est_genome_command = "est_genome  -reverse -genome $genfile -est $estfile | tee -a $output_file | "; 
    #debug line so ouput of est2genome can be looked at
    eval {
      #print (STDERR "Running command $est_genome_command\n");
      open (ESTGENOME, $est_genome_command) 
	or $self->throw("Can't open pipe from '$est_genome_command' : $!");
      
      #Use the first line to get gene orientation
      my $firstline = <ESTGENOME>;
  
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
	elsif ($_ =~ /^(Segment|Exon|Span)/) {
	  
            # "gen" = genomic sequence
            my ($primary, $score, $percent_id,
                $gen_start, $gen_end, $gen_id,
                $est_start, $est_end, $est_id) = (split)[0,1,2,
                                                         3,4,5,
                                                         6,7,8];
          
            ### Skip puny little bits? ###
            #next unless $score > 6;

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
            
            if ($primary eq 'Segment') {
                # For Segements, make sure that the est_match
                # is on the forward strand, because est_strand is not
                # stored in the database and we can't view the
                # alignment for est_strand == -1 matches.
                if ($est_strand == -1) {
                    $est_strand = 1;
                    $gen_strand *= -1;
                }
            } else {
                # For Spans and Exons, put the genomic feature
                # feature on the correct strand for the gene
                if ($estOrientation == -1 and $gen_strand == 1) {
                    $gen_strand = -1;
                    $est_strand *= -1;
                }
            }
            
	    $source_tag, 
	    $gen_strand, $est_strand,
 	    $self->_createfeatures ($score, $percent_id,
				    $gen_start, $gen_end, $gen_id, 
				    $est_start, $est_end, $est_id,
				    $source_tag, 
				    $gen_strand, $est_strand,   # est_strand is NOT stored in the db!
				    $primary);
       }    

      }
      if(!close(ESTGENOME)){
	$self->warn("problem running est_genome: $!\n");
	return(0);
      }

      $self->convert_output;

    };
    $self->_deletefiles($genfile, $estfile);
    #$self->_deletefiles($genfile, $estfile, $output_file);
    if ($@) {
        $self->throw("Error running est_genome [$@]\n");
    } else {
        return 1;
    }
}



=head2 _createfeatures

  Args   : $self, 
  Function  : create a feature pair using the data provided abouyt the two features
  Returntype: feature pair
  Exceptions: none
  Caller    : 
  Example   : 

=cut




sub _createfeatures {
    my ($self, $f1score, $f1percent_id, $f1start, $f1end, $f1id, $f2start, $f2end, $f2id, $f1source, $f1strand, $f2strand, $f1primary) = @_;
    
   
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

    push(@{$self->{'_fplist'}}, $fp);
    #print "the new fp has ".$fp->source_tag." ".$fp->primary_tag."\n";
}



=head2 convert_output

  Arg [1]   : self 
  Function  : splits the feature pair list in to spans, exons and segments and returns an array of segments
  Returntype: array of featurepairs
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub convert_output {
  my ($self) = @_;
  my @genes;
  my @exons;
  my @supp_feat;

  # split the different features up
  foreach my $f(@{$self->{'_fplist'}}){
    if ($f->primary_tag eq 'Span'){
      push(@genes, $f);
    }
    elsif($f->primary_tag eq 'Exon'){
      push(@exons, $f);
    }
    elsif($f->primary_tag eq 'Segment'){
      push(@supp_feat, $f);
    }
  }
  
 
  
  push(@{$self->{'_output'}},@supp_feat);

}


1;
