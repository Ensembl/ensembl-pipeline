#####note this is not a functional module yet drastic changes are being made to the way it is running donnot try and run this file#####



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

    $self->{'_genomic_seq'} = undef;
    $self->{'_input_features'} = [];
    $self->{'_forward_features'} = [];
    $self->{'_reverse_features'} = [];
    $self->{'_output_features'} = [];
    $self->{'_pfam_ids') = [];
    $self->{'_hmmfetch'} = undef;
    $self->{'_uniquepfamids'} = [];
    $self->{'_hmmdb'} = undef;
    print "args = @args\n";
    my ($genomic, $features, $hmmfetch, $hmmdb) = $self->_rearrange([qw(GENOMIC
									FEATURES
									HMMFETCH
									HMMDB)], @args);
    $self->throw("No genomic sequence input") unless defined($genomic);
    $self->genomic_sequence($genomic);
    $self->add_input_features($features);
    
    if ($hmmfetch){
	$self->hmmfetch($hmmfetch);
    } else {
	$self->hmmfetch('hmmfetch');
    }
    
    if($hmmdb){
	$self->hmmdb($hmmdb);
    } else {
	$self->hmmdb('/usr/local/ensembl/data/blastdb/Ensemb/Pfam');
	
    }
    
    $self->sort_input_features();

    return $self; # success - we hope!
}




#################
#GET/SET METHODS#
#################
=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
    my( $self, $value ) = @_;    
    if ($value) {
        #need to check if passed sequence is Bio::Seq object
        $value->isa("Bio::PrimarySeqI") || $self->throw("Input isn't a Bio::PrimarySeqI");
        $self->{'_genomic_sequence'} = $value;
    }
    return $self->{'_genomic_sequence'};
}


sub add_input_features{
    
    my ($self, $features) = @_;
    
    if (ref($features) eq "ARRAY") {
      push(@{$self->{'_input_features'}},@$ids);
    } else {
      $self->throw("[$ids] is not an array ref.");
    }
      


}

sub all_input_features{

  my ($self) = @_;

  return @{$self->{'_input_features'}};

}


sub sort_input_features{

  my ($self) = @_;

  @features = $self->all_input_features();

  

}

sub forward_features{
}

sub reverse_features{

  my ($self, $f) = @_;
 
  if (defined($f))
    {
	push(@{$self->{'_reverse_features'}}, $f);
    }

  return @{$self->{'_reverse_features'}};
}

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

  
  

 
  }
  
  
  return @feat;

}
1;
