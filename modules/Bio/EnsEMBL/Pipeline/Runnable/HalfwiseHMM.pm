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

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM->new(
					     -genomic => $seq,
					     -features => \@features,
								);
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

Finds which pfam domains a particular swissprot hit matchs, finds the related hmms and runs genewiseHmm

=head1 CONTACT

lec@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM;


use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter;
use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Bio::Tools::BPlite;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::Finished_GenewiseHmm;
use DB_File;
use Fcntl;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2  new

    Arg      : Bio:Seq obj, reference to array of swall features, location of hmmfetch program, location of hmmdb, location of dbm index
    Function : Make a new HalfwiseHMM object defining the above variables
    Exception: Will throw an exception if no genomic sequence is defined or no features are passed
    Caller   : 
    Example  : $halfwise = Bio::EnsEMBL::Pipeline::Runnable::HalfwiseHMM->new(genomic => $seq
                                                                              features => \@features);

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_query'} = undef;
    $self->{'_input_features'} = [];
    $self->{'_output_features'} = [];
    $self->{'_hmmfetch'} = undef;
    $self->{'_hmmdb'} = undef;
    $self->{'_hmmfilename'} = undef; #file to store protein profile hmms 
    $self->{'_errorfile'} = undef; #file to store any errors hmmfetch throw
    $self->{'_dbm_file'} = undef;
    #print STDERR "args = @args\n";
    my ($genomic, $features, $hmmfetch, $hmmdb, $dbmfile, $memory) = $self->_rearrange([qw(QUERY
											   FEATURES
											   HMMFETCH
											   HMMDB 
											   DBMFILE
											   MEMORY)], @args);
    $self->throw("No genomic sequence input") unless defined($genomic);
    $self->query($genomic);
    $self->throw("No features input") unless defined($features);
    $self->add_input_features($features);
    
    if ($hmmfetch){
	$self->hmmfetch($hmmfetch);
    } else {
	$self->hmmfetch('hmmfetch');
    }
    
    if($hmmdb){
	$self->hmmdb($hmmdb);
    } else {
	$self->hmmdb('/usr/local/ensembl/data/blastdb/Ensembl/Pfam_ls');
	
    }
    if($dbmfile){
      $self->dbm_file($dbmfile);
    } else {
      $self->dbm_file('/usr/local/ensembl/data/swiss2prot_7.0.index');	
    }
    $self->memory ($memory)    if (defined($memory));
    return $self; # success - we hope!
}




#################
#GET/SET METHODS#
#################

=head2 query

  Arg      : Bio:Seq object
  Function : get/set method for query
  Exception: throws an exception if not passed a Bio:PrimarySeq
  Caller   : 

=cut

sub query {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_query'} = $value;
    }
    return $self->{'_query'};
}

=head2 add_input_features

    Arg      : reference to array of features
    Function : get/set method for input features
    Exception: throws an exception if not passed an array ref
    Caller   : 

=cut

sub add_input_features{
    
    my ($self, $features) = @_;
    
    if (ref($features) eq "ARRAY") {
      push(@{$self->{'_input_features'}},@$features);
    } else {
      $self->throw("[$features] is not an array ref.");
    }
      
}

=head2 all_input_features

    Arg      : none
    Function : returns all input features
    Exception: none
    Caller   : 

=cut


sub all_input_features{

  my ($self) = @_;

  return @{$self->{'_input_features'}};

}


#methods that probably won't be needed


=head2 hmmfetch

    Arg      : location of hmmfetch program
    Function : gets/sets location of hmmfetch program
    Exception: none
    Caller   : 

=cut

sub hmmfetch {
    my ($self, $args) =@_;

    if (defined($args)){
	$self->{'_hmmfetch'} = $args;
    }
    return $self->{'_hmmfetch'};
}

=head2 hmmdb

    Arg      : location of hmmdb 
    Function : gets/sets location of hmmdb
    Exception: none
    Caller   : 

=cut

sub hmmdb {
    my ($self, $args) =@_;

    if (defined($args)){
	$self->{'_hmmdb'} = $args;
    }
    return $self->{'_hmmdb'};
}


=head2 hmmfilename

    Arg      : hmm filename
    Function : gets/sets hmmfilename
    Exception: none
    Caller   : 

=cut

sub hmmfilename {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_hmmfilename'} = $args;
    }

    return $self->{'_hmmfilename'};

}
=head2 dbm_file

    Arg      : dbm filename
    Function : gets/sets dbm filename
    Exception: none
    Caller   : 

=cut


sub dbm_file {
    my ($self, $args) = @_;
    if (defined ($args)){
	$self->{'_dbm_file'} = $args;
    }

    return $self->{'_dbm_file'};

}

=head2 memory

    Arg      : value memory to be set to, this is an option of how much memory genewise can use when it is run
    Function : gets/sets memory
    Exception: none
    Caller   : 

=cut

sub memory {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_memory'} = $arg;
    }

    return $self->{'_memory'} || 400000;
}

##########
#ANALYSIS#
##########

=head2 run

    Arg      : directory if to be set to anything other than tmp 
    Function : runs functions which run genewisehmm
    Exception: none
    Caller   : 

=cut

sub run {

    my ($self, $dir) = @_;
    #print "running halfwise\n";
    my %swissprot_ids = $self->get_swissprot_ids();
    #print "got swiss prot ids \n";
    my @keys = keys(%swissprot_ids);
    
    #foreach my $key (@keys){
    #  print "hid = ".$key." strand = ".$swissprot_ids{$key}."\n";
    #}

    my %pfam_ids = $self->get_pfam_ids(\%swissprot_ids);
    #print "got pfam ids \n";
    $self->create_genewisehmm(\%pfam_ids, $dir);
    
}

=head2 get_swissprot_ids

    Arg      : none
    Function : gets swissprot ids and strands from input features and returns a hash of these
    Exception: none
    Caller   : 

=cut

### 1, -1 and 0 are used to repesent strand, 1 and -1 have the same meaning as in teh ensembl databases 1 is the forward strand and -1 is the reverse 
### 0 is used if a particular id has hit on both strands
sub get_swissprot_ids{

  my ($self) = @_;

  my @features = $self->all_input_features();
  my %swissprot_ids;
  foreach my $f(@features){
    #print "swissprot id = ".$f->hseqname."\n";
    if(!$swissprot_ids{$f->hseqname}){
     
      $swissprot_ids{$f->hseqname} = $f->strand;
    
    }elsif(defined$swissprot_ids{$f->hseqname}){
      
      if($swissprot_ids{$f->hseqname} != $f->strand){
	$swissprot_ids{$f->hseqname} = 0;
      }else{
	$swissprot_ids{$f->hseqname} = $f->strand;
      }

    }
    
  }

  return %swissprot_ids;
    
}

=head2 get_pfam_ids

    Arg      : reference to a hash of swissprotids and stands
    Function : gets all pfam ids for a particular swissprot id and puts them in a hash along with the strand the swissprot id is found on and returns its
    Exception: warns if swissprot id has no pfam domains
    Caller   : 

=cut

sub get_pfam_ids{

  my ($self, $swissprot_ids) = @_;

  my %swiss_pfam;
  my %pfams;
  my $permission = O_RDONLY;
  tie(%swiss_pfam, "DB_File", $self->dbm_file, $permission, 0400) or die "couldn't open ".$self->dbm_file." : $!\n";;
  
  foreach my $swiss_id(keys(%$swissprot_ids)){
    #print "searching for pfam ids for ".$swiss_id."\n";  
    my $strand = $swissprot_ids->{$swiss_id};
    
    my $pfam = $swiss_pfam{$swiss_id};
    #print "found these pfam domains ".$pfam." \n";
    if(!$pfam){
      #$self->warn("the swissprot ".$swiss_id ."entry has no pfam domains :!")
    }
    my @pfams = split /,/, $pfam;
    #print "pfams = @pfams\n";
    foreach my $pfam_id(@pfams){
      #print "trying to create a hash entry for pfam".$pfam_id."\n";
      if(!defined $pfams{$pfam_id}){
	#print "defining pfam hash key\n";
	$pfams{$pfam_id} = $strand;
    
     }elsif(defined$pfams{$pfam_id}){
	if($pfams{$pfam_id} != $strand){
	 $pfams{$pfam_id} = 0;
	}else{
	 $pfams{$pfam_id} = $strand;
       }
      }
    }
  }
 # my @ids = keys(%pfams);
 # foreach my $id(@ids){
 #   print $id."\n";
 # }
  return %pfams;
}
 
=head2 create_genewisehmm

    Arg      : reference to pfam_id hash and directory if it is to be set to anything other than temp
    Function : runs through each pfam id runs, get_hmm and creates and runs the GenewiseHmm runnable
    Exception: throws an exception if anything other than 1, -1 or 0 is found for strand
    Caller   : 

=cut

sub create_genewisehmm{

 my ($self, $pfam_ids, $dir) = @_;
 #print "getting hmms\n";
 $self->workdir('/tmp') unless ($self->workdir($dir)); 
 $self->checkdir();
 #print "running HalfwiseHMM\n";
 #system("pwd");
 my @ids = keys(%$pfam_ids);
 my $genewisehmm;
 my @features;
 my $memory = $self->memory;
 my %pfam_run;
 #print "there are ".scalar(@ids)." pfam ids to run\n";
 foreach my $id (@ids){
   
   my $strand = $pfam_ids->{$id};

   $self->get_hmm($id);
    if (-z $self->hmmfilename){
    $self->warn("hmm file not created :$!");
    next;
  }
   #print $id."\n";
   if($pfam_run{$id}){
     #print "running twice on ".$id."\n";
   }else{
     $pfam_run{$id} = 1;
   }
   if($strand == 1){
     @features = [];
     #print "running on ".$strand."\n";
     $genewisehmm = $self->run_genewisehmm($strand, $memory);
     $genewisehmm->run();
     @features = $genewisehmm->output();
     #$self->display_genes(\@features);
     #print "adding ".scalar(@features)." to output\n";
     unlink $self->hmmfilename();
     $self->add_output_features(\@features);
   }elsif($strand == -1){
     @features = [];
     #print "running on ".$strand."\n";
     $genewisehmm = $self->run_genewisehmm($strand, $memory);
     $genewisehmm->run();
     @features = $genewisehmm->output();
     #$self->display_genes(\@features);
     #print "adding ".scalar(@features)." to output\n";
     unlink $self->hmmfilename();
     $self->add_output_features(\@features);
  }elsif($strand == 0){
    @features = [];
    #print "running on ".$strand."\n";
    $genewisehmm = $self->run_genewisehmm(1, $memory);
    $genewisehmm->run();
    @features = $genewisehmm->output();
    #$self->display_genes(\@features); 
    #print "adding ".scalar(@features)." to output\n";
    $self->add_output_features(\@features);
    $genewisehmm = $self->run_genewisehmm(-1, $memory);
    $genewisehmm->run();
    @features = [];
    @features = $genewisehmm->output();
    #$self->display_genes(\@features);
    #print "adding ".scalar(@features)." to output\n";
    unlink $self->hmmfilename();
    $self->add_output_features(\@features);
  }else{
    $self->throw("strand = ".$strand." something funnies going on : $!\n");
  }

 }


}

=head2 run_genewisehmm

    Arg      : the strand of the hit and the memory to be used by genewisehmm
    Function : runs the new method of GenewiseHmm and returns the runnable
    Exception: none
    Caller   : 

=cut

sub run_genewisehmm{

  my ($self, $strand, $memory) = @_;
  my $reverse;
  if($strand == -1){
    $reverse = 1;
  }else{
    $reverse = undef;
  }
  #print "creating genewisehmm\n";
  my $genewisehmm = Bio::EnsEMBL::Pipeline::Runnable::Finished_GenewiseHmm->new('-query' => $self->query(),
                                                                    '-memory' => $memory,
                                                                    '-hmmfile' => $self->hmmfilename(),
                                                                    '-reverse' => $reverse);

  return $genewisehmm

}


=head2 get_hmm

    Arg      : id of hmm to be fetched normally a pfam id 
    Function : runs hmmfetch on given id
    Exception: thows if no hmmfile is created
    Caller   : 

=cut

sub get_hmm{
  
  my ($self, $id) = @_;
  #print "getting hmms\n";
  $self->hmmfilename($id.".".$$.".hmm");
  #print "getting hmm for id ".$id."\n";
  my $command =  $self->hmmfetch." ".$self->hmmdb." ".$id." > ".$self->hmmfilename;

  #print "command = ".$command."\n";
  #system ('pwd');
  eval{
    system($command);
  };
  if($@){
    $self->warn("hmmfetch threw error : $@\n");
  }
  
  if (-z $self->hmmfilename){
    $self->warn("hmm file not created :$!");
  }elsif(-e $self->hmmfilename){
    open(ERROR, $self->hmmfilename) or die "couldn't open error file : $!";
    while(<ERROR>){
      if(/no such hmm/i){
	print STDERR "error message : ".$_."\n";
      }
    }
    close(ERROR) or die "couldn't close error file : $!";
  }
  
}

=head2 add_output_features

    Arg      : array ref to array of output features
    Function : adds the array of features to the output feature array
    Exception: throws an exception if not passed an array ref
    Caller   : 

=cut

sub add_output_features{
    
    my ($self, $features) = @_;
    #print "adding ".scalar(@$features)." to output\n";
    if (ref($features) eq "ARRAY") {
      push(@{$self->{'_output_features'}},@$features);
    } else {
      $self->throw("[$features] is not an array ref.");
    }
      
    #print "total feature no = ".scalar(@{$self->{'_output_features'}})."\n";

}




################
#OUTPUT METHODS#
################


=head2 output

    Arg      : none
    Function : returns the array of output features
    Exception: none
    Caller   : 

=cut


sub output{
  my ($self) = @_;
  #print "outputing data\n";
  #print "returning ".scalar(@{$self->{'_output_features'}})." to output\n";
  return @{$self->{'_output_features'}};

}

sub display_genes {
  my ($self, $result) = @_;
  #Display output
  my @results = @$result;
  foreach my $obj (@results)
    {
      print STDERR ("gene:  ".$obj->gffstring. "\n");
      if ($obj->sub_SeqFeature)
	{
	  foreach my $exon ($obj->sub_SeqFeature)
	    {
	      print STDERR "Exon:  ".$exon->gffstring."\n";
	      if ($exon->sub_SeqFeature)
		{
		  foreach my $sub ($exon->sub_SeqFeature){
		    print STDERR "supporting features:  ".$sub->gffstring."\n";
		  }
		}
	    }
	}
    }
}


1;
