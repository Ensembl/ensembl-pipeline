
#
#
# Cared for by Laura Clarke  <lec@sanger.ac.uk>
#
# Copyright Laura Clarke
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM

=head1 SYNOPSIS

  #create and fill Bio::Seq object
  my $seqfile = '/nfs/disk65/mq2/temp/bA151E14.seq'; 
  my $seq = Bio::Seq->new();
  my $seqstream = Bio::SeqIO->new(-file => $seqfile, -fmt => 'Fasta');
  $seq = $seqstream->next_seq();
  #create Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM object
  my $halfwisehmm = Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM->new (-QUERY => $seq);
  $halfwisehmm->workdir($workdir);
  $halfwisehmm->run();
  my @genes = $halfwisehmm->output();
  my @exons = $halfwisehmm->output_exons();
  my $seqfeature = $halfwisehmm->output_singlefeature();

=head1 DESCRIPTION

This package is based on the genscan runnable.
HalfwiseHMM takes a Bio::Seq (or Bio::PrimarySeq) object and runs HalfwiseHMM on it. The
resulting output is parsed to produce a set of Bio::SeqFeatures. 

=head2 Methods:

=over 4

=item new($seq_obj)

=item    HalfwiseHMM($path_to_HalfwiseHMM)

=item    workdir($directory_name)

=item    run()

=item    output()

=back

=head1 CONTACT

ensembl-dev.ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 


=cut

package Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM;


use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter;
use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use Bio::Tools::BPlite;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;



BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM->new (-QUERY => $seq);
    Function:   Initialises HalfwiseHMMobject
    Returns :   a HalfwiseHMMObject
    Args    :   A Bio::Seq object 
                (blastdb, hmmdb location of exes etc  optional)

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_query'} = undef; #location of Bio::Seq object
    $self->{'_blastprog'} = undef; #location of blast program
    $self->{'_blastdb'} = undef; #location fo blast database
    $self->{'_blastobj'} = undef; #location of blast object
    $self->{'_threshold'} = undef; #threshold value for blast analysis
    $self->{'_blastfplist'} = [];        # an array of feature pairs (the output)
    $self->{'_pfamlist'}    = []; #list of pfam hits
    $self->{'_uniquepfamids'} = [];
    $self->{'_featurelist'} = []; #array of features
    $self->{'_genewise'} = undef; #location of genewise prog
    $self->{'_options'} = undef; #options for genewise program
    $self->{'_workdir'}   = undef;     # location of temp directory
    $self->{'_filename'}  = undef;     # file to store Bio::Seq object
    $self->{'_hmmfilename'} = undef; #file to store protein profile hmms 
    $self->{'_errorfile'} = undef; #file to store any errors hmmfetch throws 
    $self->{'_results'}   = undef;     # file to store results of analysis
    $self->{'_hmmfetch'} = undef; #location of hmmfetch program
    $self->{'_hmmdb'} = undef; #location of hmmdb
    
    print STDERR "args: ", @args, "\n";
    my ($query, $blastprog, $blastdb, $threshold, $genewise, $options, $hmmfetch, $hmmdb) 
	= $self->_rearrange([qw (QUERY BLASTPROG BLASTDB THRESHOLD GENEWISE OPTIONS HMMFETCH HMMDB)], @args);
    
    $self->query($query) if ($query);
    if ($blastprog){
       $self->blast_loc($self->find_executable($blastprog));
    } else {
      $self->blast_loc($self->find_executable('wublastx'));
    }
    if ($blastdb) {
	$self->blastdb($blastdb);
    } else {#note database name must always contain halfwise
	$self->blastdb("/usr/local/ensembl/halfwise.db");
	#$self->blastdb("/nfs/disk100/pubseq/halfwise/halfwise.db");
    }
    if (defined($genewise )){
	$self->genewise_loc($self->find_executable($genewise));
    } else {
	$self->genewise_loc($self->find_executable('genewisedb'));
    }
    if (defined($threshold)) {
	$self->threshold($threshold);
    }
    else
    {
	$self->threshold(1e-6);
    }
    if ($options) {
      $self->options($options);
    } else {
      $self->options(' -init wing -cace -cut 25 -aln 200 -gff -quiet ');  
    }

    if ($hmmfetch){
	$self->hmmfetch($hmmfetch);
    } else {
	$self->hmmfetch('hmmfetch');
    }
    
    if($hmmdb){
	$self->hmmdb($hmmdb);
    } else {
	$self->hmmdb('/usr/local/ensembl/data/Pfam');
	#$self->hmmdb('/nfs/disk100/pubseq/Pfam/DB/Pfam');
    }
    return $self; # success - we hope!
}




#################
#GET/SET METHODS#
#################

=head2 query

    Title   :   query
    Usage   :    $HalfwiseHMM->query($seq);
    Function:   sets the sequence the halfwisehmm object will run on
  and checks it is a Bio::Seq
    Returns :   a seq
    Args    :   A Bio::Seq object 
                

=cut

sub query {
    my ($self, $seq) = @_;
    if ($seq)
    {
        unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")) 
        {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_query'} = $seq ;
        $self->filename($self->query->id.".$$.seq");
        $self->results($self->filename.".hlf");
	$self->hmmfilename($self->filename."dbhmm.$$");
	$self->errorfile($self->filename."err.$$");
    }
    return $self->{'_query'};
}

=head2 blast_loc

 Title    : blast_loc
 Usage    : $self->blast_loc($self->find_executable('wublastx'))
 Function : checks exceutable is found in correct location
 Returns  : executable location
 Args     : executable location

=cut

sub blast_loc {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("Blast not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_blastprog'} = $location ;
    }
    #print STDERR "blast location is ".$self->{'_blastprog'}."\n";
    return $self->{'_blastprog'};
}


=head2 blastdb

 Title    : blastdb
 Usage    : $self->blastdb($db);
 Function : gets and sets balst db location
 Returns  : blastdb location
 Args     : blastdb location

=cut

sub blastdb {
    my ($self, $db) = @_;

    if (defined($db)) {
	$self->{'_blastdb'} = $db ;
    }
    return $self->{'_blastdb'};
}


=head2 genewise_loc

 Title    : genewise_loc
 Usage    : $self->genewise_loc($self->find_executable('wublastx'))
 Function : checks exceutable is found in correct location
 Returns  : executable location
 Args     : executable location

=cut


sub genewise_loc {
    my ($self, $location) = @_;
    if ($location)
    {
        $self->throw("genewise not found at $location: $!\n") 
                                                    unless (-e $location);
        $self->{'_genewise'} = $location ;
    }
    return $self->{'_genewise'};
}


=head2 options

 Title    : options
 Usage    : $self->options($options);
 Function : gets and sets the options which may be passed to genewise
 Returns  : the options
 Args     : the options

=cut


=head2 hmmfetch 

 Title    : hmmfetch
 Usage    : $self->hmmfetch($location);
 Function : gets and sets the location of the hmmfetch command
 Returns  : the hmmfetch command
 Args     : the hmmfetch comand

=cut

sub hmmfetch {
    my ($self, $args) =@_;

    if (defined($args)){
	$self->{'_hmmfetch'} = $args;
    }
    return $self->{'_hmmfetch'};
}

=head2 hmmdb

 Title    : hmmdb
 Usage    : $self->hmmdb($db);
 Function : gets and sets the location of the hmmdb
 Returns  : the location of the hmmdb
 Args     : the location of the hmmdb

=cut

sub hmmdb {
    my ($self, $args) =@_;

    if (defined($args)){
	$self->{'_hmmdb'} = $args;
    }
    return $self->{'_hmmdb'};
}


=head2  hmmfilename

 Title    : hmmfilename
 Usage    : $self->hmmfilename($filename);
 Function : gets and sets the name of the hmmdb file
 Returns  : the name of the hmmdb file
 Args     : the name of the hmmdb file

=cut

sub hmmfilename {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_hmmfilename'} = $args;
    }

    return $self->{'_hmmfilename'};

}


=head2  errorfile

 Title    : errorfile
 Usage    : $self->errorfile($filename);
 Function : gets and sets the name of the hmmdb file
 Returns  : the name of the hmmdb file
 Args     : the name of the hmmdb file

=cut

sub errorfile {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_errorfile'} = $args;
    }

    return $self->{'_errorfile'};

}



##########
#ANALYSIS#
##########


=head2  run
 
 Title    : Run
 Usage    : $self->run
 Function : runs blast and halfwise
 Returns  : nothing
 Args     : workdir

=cut

sub run {

    my ($self) = @_;
    

    $self->workdir('/tmp') unless $self->workdir();
    $self->checkdir();

    #print STDERR "write sequence to file\n";
    #system ('pwd');
    $self->writefile();
    #print STDERR "written sequence to file\n";
    my $filesize = (-s $self->filename);

    if (!defined $filesize || $filesize == 0)
      {
	print STDERR "results file ".$self->filename." has not been created or is empty\n";
      } #else {
	#print $self->filename." is ".$filesize." big\n";
      #}
    #system ('pwd');
    my $blast = $self->create_blast;
    $blast->run;
    #print "run blast\n";
    $self->writefile();
    #print "recreated sequence file";
    $filesize = (-s $self->filename);
    
    if (!defined $filesize || $filesize == 0)
      {
	print STDERR "sequence file ".$self->filename." has not been created or is empty\n";
      }#else {
	#print $self->filename." is ".$filesize." big\n";
      #}
    my $fp;
    my @featurepairs = $blast->output;
    #print STDERR "blast output = ".@featurepairs."\n";
    foreach $fp (@featurepairs)
    {
      #print "fp = ".$fp."\n";
      #print "gff = ".$fp->gffstring."\n";
      #print "seqname ".$fp->hseqname."\n";
      my $pfamid = $fp->hseqname;
      $self->add_pfam_id($pfamid);
    }
    $self->find_unique_ids;
    #print "found unique ids \n";
    $self->get_hmm();
    #print "got hmms\n";
    $self->run_genewise();
    $self->parse_results();
    $self->deletefiles();
    unlink ($self->hmmfilename);
}


=head2 create_blast

 Title    : create_blast
 Usage    : $self->create_blast();
 Function : creates the blast object so the blast analysis can be run
 Returns  : nothing
 Args     : nothing

=cut

sub create_blast {
    
    my ($self) = @_;
    eval{
      $self->query;
    };
    if($@){
      $self->throw("cannot create a blast object without a seq object : ".$@."\n");
    }
    my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (-QUERY => $self->query,
							      -PROGRAM => $self->blast_loc,
							      -DATABASE => $self->blastdb,
							      -THRESHOLD => 1e-6,
							      );
    $self->{'_blastobj'} = $blast;

    return $self->{'_blastobj'};
}


=head2 blast_output

 Title    : blast_output
 Usage    : $self->blast_output;
 Function : gets and sets the blast output
 Returns  : the array of features produced by blast
 Args     : the blast runnable

=cut

sub blast_output{

    my ($self, $blast) = @_;

    if (defined($blast))
    {
	@{$self->{'_blastfplist'}} = $blast->output;
    }
    return @{$self->{'_blastfplist'}};
}


=head2  add_pfam_id

 Title    : add_pfam_id
 Usage    : $self->add_pfam_id($id);
 Function : adds a pfam id to the array of ids
 Returns  : the array of ids
 Args     : a pfam id;

=cut

sub add_pfam_id{
    my ($self, $id) = @_;
    if (defined($id))
    {
	push(@{$self->{'_pfamlist'}}, $id);
    }
    return @{$self->{'_pfamlist'}};
}


=head2 

 Title    : each_pfam_id
 Usage    : @pfamids = $self->each_pfam_id;
 Function : returns the array of pfam ids
 Returns  : an array of pfam ids
 Args     : none

=cut

sub each_pfam_id{
    my ($self) = @_;
    return @{$self->{'_pfamlist'}};
}


=head2 

 Title    :find_unique_ids
 Usage    :$self->find_unique_ids;
 Function :ensures all the pfam ids wich will be passed to hmmfetch will be unique so the same analysis isnt run multiple times
 Returns  : nothing
 Args     : nothing

=cut

sub find_unique_ids{
  my ($self) = @_;
  my @pfam = $self->each_pfam_id;
  my %hash;
  foreach my $id(@pfam){
    if(!defined($hash{$id})){
      $hash{$id} = 1;
    }
  }
  foreach my $id(keys %hash){
      $self->add_unique_pfam_id($id);
    }
}     
     

=head2 add_unique_pfam_id

 Title    : add_unique_pfam_id
 Usage    : $self->add_unique_pfam_id($id);
 Function : sets up an array of unique pfam ids 
 Returns  : the array of unique ids
 Args     : a pfam id preferable unique

=cut

sub add_unique_pfam_id{
    my ($self, $id) = @_;
    if (defined($id))
    {
     # print "id = ".$id."\n";
	push(@{$self->{'_uniquepfamids'}}, $id);
    }
    return @{$self->{'_uniquepfamids'}};
}


=head2 each_unique_pfam_id

 Title    : each_unique_pfam_id
 Usage    : @ids = $self->each_unique_pfam_id
 Function : returns all the unique pfma ids
 Returns  : an array of pfam ids
 Args     : none

=cut

sub each_unique_pfam_id{
    my ($self) = @_;
    return @{$self->{'_uniquepfamids'}};
}


=head2 

 Title    : get_hmm
 Usage    : $self->get_hmm;
 Function : gets all the hmms for the pfam matchs found by blast
 Returns  : none
 Args     : none

=cut

sub get_hmm{
  my ($self) = @_;
  #open(TEMPDB, ">".$self->hmmfilename) or die "couldn't open ".$self->hmmfilename." : $!\n";
  my $count = 0;
  foreach  my $id  ( $self->each_unique_pfam_id) {
    $id =~s/^(.*?\;).*/$1/;
    $count++;
    #print  "Loading $id\n";
    #print "command = ".$self->hmmfetch." ".$self->hmmdb." ".$id." 2>> ".$self->errorfile."1>>".$self->hmmfilename." \n ";
    system($self->hmmfetch." ".$self->hmmdb." ".$id." 2>> ".$self->errorfile."1>>".$self->hmmfilename)
    #open(GETZ,  $self->hmmfetch." ".$self->hmmdb." ".$id." | ") or die "couldn't open hmmfetch : $! \n";
    #while(<GETZ>) {
     # print TEMPDB $_;
    #}
   # close(GETZ) or die "couldn't close hmmfetch : $!\n";
  }
  #close(TEMPDB) or die "couldn't open hmmfile : $!\n";
  #print "starting error analysis\n";
  if (-e $self->errorfile){
    open(ERROR, $self->errorfile) or die "couldn't open error file : $!";
    while(<ERROR>){
      if(/no such hmm/i){
	print STDERR "error message : ".$_."\n";
      }
    }
    close(ERROR) or die "couldn't close error file : $!";
  }
}


=head2 run_genewise

 Title    : run_genewise
 Usage    : $self->run_genewise
 Function : runs genewisedb wotht the pfam hmms and the clone
 Returns  : none
 Args     : none

=cut

sub run_genewise {

 my ($self) = @_;
 #print STDERR "trying to run genewise\n";
 #print STDERR "systems command is ".$self->genewise_loc." -pfam ".$self->hmmfilename." -dnas ".$self->filename ." ".$self->options." > ".$self->results."\n";
 my $filesize = (-s $self->hmmfilename); #check that hmmfile has been creeated

 if (!defined $filesize || $filesize == 0)
   {
     print STDERR "results file ".$self->hmmfilename." has not been created or is empty\n";
   }
 #system ('pwd'); 
 system ($self->genewise_loc.' -pfam /tmp/'.$self->hmmfilename.' -dnas /tmp/'.$self->filename .' '.$self->options.' > '.$self->results);
    
 $filesize = (-s $self->results);#check that results file exists

 if (!defined $filesize || $filesize == 0)
   {
     print STDERR "results file ".$self->results." has not been created or is empty\n";
   }

}


=head2 parse_results

 Title    :parse_results
 Usage    :$self->parse_results
 Function :parse the results produced by genewisedb
 Returns  :none
 Args     :none

=cut

sub parse_results{
  
  my ($self) = @_;

  my $filehandle;
  if (ref ($self->results) !~ /GLOB/){
    open (FH, "<".$self->results) or $self->throw ("couldn't open file".$self->results.":$!\n");
    $filehandle = \*FH;
  } else {
    $filehandle = $self->results;
  }
  #print "opened results file\n";
  while(<$filehandle>){

    my %feature;
    my $count=0;
    my $match=0;
    if(/^>/){
      
      while(<$filehandle>){
	
	#print "1 line to be parsed = ".$_."\n"; 
	if(/>Results/i){
	    $match =0;
	    $count++;
	    #print "another pfam domain matched setting count to ".$count." \n";
	  }elsif(/genewise-prediction/i){ 
	    #print "2 line to be parsed = ".$_."\n"; 
	    my @elements = split;
	    # my $count =0;
	    # foreach my $element(@elements){
	    #   print "element ".$count." = ".$element." \n";
	    #   $count++;
	    # }
	    
	    $self->throw("unable to properly parse genewisedb output there are ".scalar(@elements)." elements : $!\n")
	      unless(scalar(@elements) == 9);
	    #print "3 line to be parsed = ".$_."\n";   
	    
	    if(/match/i){
	      
	      $feature{score} = $elements[5];
	      
	    } elsif (/cds/i) {
	      if($elements[6] eq '+'){
		$feature{strand} = 1;
		$feature{start} = $elements[3];
		$feature{end} = $elements[4];
	      }elsif($elements[6] eq '-'){
		$feature{strand} = -1;
		$feature{start} = $elements[4];
		$feature{end} = $elements[3];
	      }
	      $feature{frame} = $elements[7];
	      $feature{seqname} = $elements[0];
	      $elements[8] =~ s/-genewise-prediction-/\./i;
	      $feature{hseqname} = $elements[8].".".$count.".".$match;
	      #print "$count = ".$count."\n";
	      #print "hseqname = $feature{hseqname}\n";
	      $match++;
	      $feature{type} = $elements[2];
	      #     print "checking hash contains data\n";
	      #     foreach my $key(keys %feature){
	      #	my $value = $feature{$key};
	      #	print $key."  = ".$value."\n";
	      #}
	      $self->create_feature(\%feature);
	      #print "created feature\n";
	    }
	  }
      }
    }
  }
}



=head2 create_feature

 Title    :create_feature
 Usage    :$self->create_feature(\%feats);
 Function :creates features from an hash of results
 Returns  : none
 Args     : reference to hash of features

=cut

sub create_feature{

  my ($self, $feat) = @_;
  #print "creating feature\n";

  
    #create analysis object
    my $analysis_obj = Bio::EnsEMBL::Analysis->new
                        (   -db              => undef,
                            -db_version      => undef,
                            -program         => 'halfwise',
                            -program_version => '1.0',
                            -gff_source      => 'halfwise',
                            -gff_feature     => 'similarity');

  #create and fill Bio::EnsEMBL::Seqfeature objects
  #print "seqname = ". $feat->{'seqname'}."\n start = ".$feat->{'start'}."\n end = ".$feat->{'end'}."\n";
  #print "strand >". $feat->{'strand'}."\n score >". $feat->{'score'}."\n frame = ".$feat->{'frame'}."\n";
  #print "hseqname > ".$feat->{'hseqname'}."\n source_tag > 'genewisedb'\n  primary_tag => 'similarity'\n, analysis = ". $analysis_obj."\n"; 
  my $feature = Bio::EnsEMBL::SeqFeature->new
    (   -seqname => $feat->{'hseqname'},
	-start   => $feat->{'start'},
	-end     => $feat->{'end'},
	-strand  => $feat->{'strand'},
	-score   => $feat->{'score'},
	-source_tag  => 'genewisedb',
	-primary_tag => 'similarity',
	-analysis => $analysis_obj);  
  
  #print "feature = ".$feature."\n";
  push(@{$self->{'_featurelist'}}, $feature);

}





################
#OUTPUT METHODS#
################


=head2 

 Title    : ouptut
 Usage    : $self->output
 Function : turns the features into feature pairs
 Returns  : an array of feature pairs
 Args     : none

=cut

sub output {

  my ($self) = @_;

  my @feat;

  my $analysis = Bio::EnsEMBL::Analysis->new(   -db              => undef,
						-db_version      => undef,
						-program         => 'halfwise',
						-program_version => '1.0',
						-gff_source      => 'halfwise',
						-gff_feature     => 'similarity',
						-logic_name      => 'halfwise',
						);

  
  

  my @feature = @{$self->{'_featurelist'}};
  
  foreach my $feat (@feature) {
    my $f = new Bio::EnsEMBL::SeqFeature(-seqname => $feat->seqname,
					 -start   => $feat->start,
					 -end     => $feat->end,
					 -strand  => $feat->strand,
					 -frame => $feat->frame,
					 -score   => $feat->score,
					 -source_tag => 'halfwise',
					 -primary_tag => 'similarity',
					 -analysis     => $analysis);
    my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $feat->seqname,
					  -start   => $feat->start,
					  -end     => $feat->end,
					  -strand  => $feat->strand,
					  -frame => $feat->frame,
					  -score   => $feat->score,
					  -source_tag => 'halfwise',
					  -primary_tag => 'similarity',
					  -analysis     => $analysis);
    

    $f->analysis($analysis);
    $f2->analysis($analysis);

    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f,
					   -feature2 => $f2);
    
    push(@feat, $fp);
  }
  
  
  return @feat;

}
1;
