=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Finished_HalfwiseHMM

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Finished_HalfwiseHMM->new(
                                                                          -genomic => $seq,
					                                 -features => \@features,
								         );
    $obj->run

    my @newfeatures = $obj->output;

=head1 DESCRIPTION

Finds which pfam domains a particular swissprot hit matchs, finds the related hmms and runs genewiseHmm

=head1 CONTACT

rds@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Finished_HalfwiseHMM;


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

use Fcntl;
use Data::Dumper;

our @ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head1  new

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
    my ($genomic, $features, $hmmfetch, $hmmdb, $pfam_ids, $memory, $options, $program) = $self->_rearrange([qw(QUERY
											   FEATURES
											   HMMFETCH
											   HMMDB 
											   PFAMIDS
											   MEMORY
											   OPTIONS
											   PROGRAM)], @args);
    $self->throw("No genomic sequence input") unless defined($genomic);
    $self->query($genomic);
    $self->throw("No features input") unless defined($features);
    $self->add_input_features($features);
    $self->options($options);
    $self->program($program);
    $self->pfam_ids($pfam_ids);
    
    if ($hmmfetch){
	$self->hmmfetch($hmmfetch);
    } else {
	$self->hmmfetch('/usr/local/ensembl/bin/hmmfetch');
    }
    
    if($hmmdb){
	$self->hmmdb($hmmdb);
    }else{
        # Should not get here as the default is hard coded into RunnableDB.
	$self->hmmdb('/data/blastdb/Ensembl/Pfam_ls');
    }
    $self->memory($memory) if (defined($memory));
    return $self; # success - we hope!
}

#################
#GET/SET METHODS#
#################

sub pfamDB{
    my ($self, $dbobj) = @_;
    $self->{'_pfamDB'} = $dbobj if $dbobj;
    $self->throw("Not a Bio::EnsEMBL::DBSQL::DBAdaptor")
        unless $self->{'_pfamDB'}->isa("Bio::EnsEMBL::DBSQL::DBConnection");
    return $self->{'_pfamDB'};
}

sub program{
    my ($self, $program) = @_;
    $self->{'_program'} = $program if $program;
    return $self->{'_program'};
}

=head1 query

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

=head1 add_input_features

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


=head1 hmmfetch

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

=head1 hmmdb

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


=head1 hmmfilename

    Arg      : hmm filename
    Function : gets/sets hmmfilename
      This is the hmm file containing just the hmms for the pfam_ids
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

=head1 memory

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

=head1 run

    Arg      : directory if to be set to anything other than tmp 
    Function : runs functions which run genewisehmm
    Exception: none
    Caller   : 

=cut

sub run {

    my ($self, $dir) = @_;
    #print "running halfwise\n";
    my $pfam_ids = $self->pfam_ids();
    #$self->create_genewisehmm($pfam_ids, $dir);
    #$self->create_genewisehmm_individually($pfam_ids, $dir);
    $self->create_genewisehmm_complete($pfam_ids, $dir);
    #$self->throw(Dumper($pfam_ids));
}

sub pfam_ids{
    my ($self, $ids) = @_;
    if($ids){
        $self->{'_pfam_ids_strand'} = $ids;
    }
    return $self->{'_pfam_ids_strand'};
}

=head1 create_genewisehmm_individually

    Arg      : reference to pfam_id hash and directory if it is to be set to anything other than temp
    Function : runs through each pfam id runs, get_hmm and creates and runs the GenewiseHmm runnable
    Exception: throws an exception if anything other than 1, -1 or 0 is found for strand
    Caller   : 

=cut

sub create_genewisehmm_individually{
    my ($self, $pfam_ids, $dir) = @_;
    print STDERR "getting hmms\n"; ##########
    $self->workdir('/tmp') unless ($self->workdir($dir)); 
    $self->checkdir();
    my $memory = $self->memory;

    print STDERR "there are ".scalar(keys(%$pfam_ids))." pfam ids to run\n"; ##########

    while(my ($pfam_id, $strand) = each(%$pfam_ids)){
        print STDERR "doing the hmm for id: $pfam_id\n"; ##########
        $self->get_hmm($pfam_id);
        if (-z $self->hmmfilename){
            $self->warn("hmm file not created :$!");
            next;
        }
        my @strands = ();
        if($strand == 1){
            push(@strands, $strand);
        }elsif($strand == -1){
            push(@strands, $strand);
        }elsif($strand == 0){
            push(@strands, -1);
            push(@strands, 1);
        }else{
            $self->throw("strand = ".$strand." something funnies going on : $!\n");
        }
        foreach my $s (@strands){
            $self->run_genewisehmm($s);
            #print STDERR "running on strand: $s\n"; ##########
            #my $genewisehmm = $self->get_GenewiseHMM($s, $memory);
            #$genewisehmm->run();
            #my @features = $genewisehmm->output();
            #$self->display_genes(\@features); ##########
            #print STDERR "adding ".scalar(@features)." to output\n"; ##########
            #$self->add_output_features(\@features);
        }
        print STDERR "removing hmmfile: ".$self->hmmfilename()."\n"; ##########
        unlink $self->hmmfilename();
    }
    $self->throw("finished this bit");
    return undef;
}

=head1 create_genewisehmm_individually

    Arg      : reference to pfam_id hash and directory if it is to be set to anything other than temp
    Function : creates a hmm databases of pfam ids, then creates and runs the GenewiseHmm runnable
    Exception: 
    Caller   : 

=cut

sub create_genewisehmm_complete{
    my ($self, $pfam_ids, $dir) = @_;
    print STDERR "getting hmms\n"; ##########
    $self->workdir('/tmp') unless ($self->workdir($dir)); 
    $self->checkdir();

    print STDERR "there are ".scalar(keys(%$pfam_ids))." pfam ids in database\n"; ##########
    return unless scalar(keys(%$pfam_ids));
    print STDERR "doing the hmm for ids: ". join(" ", keys(%$pfam_ids)) . "\n"; ##########
    $self->make_hmmdb($pfam_ids);
    if (-z $self->hmmfilename){
        $self->warn("hmm file not created :$!");
        $self->throw("hmm file not created :$!");
    }
    $self->run_genewisehmm(0);
    print STDERR "removing hmmfile: ".$self->hmmfilename()."\n"; ##########
    unlink $self->hmmfilename();
    #$self->throw("finished this bit");
    return undef;
}

=head1 run_genewisehmm

    Arg       :
    Function  :
    Exception :
    Caller    : $self->create_genewisehmm_complete, $self->create_genewisehmm_individually
    
=cut

sub run_genewisehmm{
    my ($self, $strand) = @_;
    my $memory = $self->memory;
    print STDERR "running on strand: $strand\n"; ##########
    my $genewisehmm = $self->get_GenewiseHMM($strand, $memory);
    $genewisehmm->run();
    my @features = $genewisehmm->output();
    $self->display_genes(\@features); ##########
    print STDERR "adding ".scalar(@features)." to output\n"; ##########
    $self->add_output_features(\@features);
}

=head1 get_GenewiseHMM

    Arg      : the strand of the hit and the memory to be used by genewisehmm
    Function : runs the new method of GenewiseHmm and returns the runnable
    Exception: none
    Caller   : 

=cut

sub get_GenewiseHMM{
  my ($self, $strand, $memory) = @_;
  my $reverse = ($strand == -1 ? 1 : undef);
  print STDERR "creating genewisehmm strand $strand reverse $reverse\n"; ##########
  print STDERR "OPTIONS To Genewise: ".$self->options()."\n"; ##########
#  $genewisehmm->set_environment("/usr/local/ensembl/data/wisecfg/");
  $ENV{WISECONFIGDIR} = "/usr/local/ensembl/data/wisecfg/";

  my $genewisehmm =
    Bio::EnsEMBL::Pipeline::Runnable::Finished_GenewiseHmm->new('-query'    => $self->query(),
                                                                '-memory'   => $memory,
                                                                '-hmmfile'  => $self->hmmfilename(),
                                                                '-reverse'  => $reverse,
                                                                '-genewise' => $self->program(),
                                                                '-options'  => $self->options()
    );
  return $genewisehmm;

}


=head1 get_hmm

    Arg      : id of hmm to be fetched normally a pfam id 
    Function : runs hmmfetch on given id
    Exception: thows if no hmmfile is created
    Caller   : 

=cut

sub _get_hmm{
  
  my ($self, $id) = @_;
  #print "getting hmms\n";
  $self->hmmfilename($id.".".$$.".hmm");
  #print "getting hmm for id ".$id."\n";
  my $command =  $self->hmmfetch." ".$self->hmmdb." ".$id." > ".$self->hmmfilename;

  print STDERR "command = ".$command."\n";
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

sub make_hmmdb{
    my ($self, $pfam_ids) = @_;
    my @pfamIds = keys(%$pfam_ids);

    print STDERR "getting the hmms for ids: ". join(" ", @pfamIds) . "\n"; ##########
    my $filename = substr(join("_", @pfamIds),0,20);
    my $hmm_file = $self->hmmfilename("$filename.$$.hmmdb");

    print STDERR "I'm trying to make the hmmdb in file '" . $self->hmmfilename() . "' \n";

    foreach my $id(@pfamIds){
        my $acc_sv = $id;#"${id}.$pfam_acc_sv->{$id}";
        print STDERR "Getting hmm for pfam id '".$acc_sv."' \n";

        my $command =  $self->hmmfetch." ".$self->hmmdb." ".$acc_sv ." 2> ${hmm_file}.err >> $hmm_file";
        eval{
            print STDERR $command . "\n";
            system($command);
        };
        if($@){
            $self->warn("hmmfetch threw error : $@\n");
        }
        if (-z $hmm_file){
            $self->warn("hmm file not appended :$!");
        }
        if(-e "${hmm_file}.err"){
            open(my $error, "${hmm_file}.err") or die "couldn't open error file : $!";
            while(<$error>){
                if(/no such hmm/i){
                    print STDERR "error message : ".$_."\n";
                }
            }
            close($error) or die "couldn't close error file : $!";
        }
    }
  
}
=head1 add_output_features

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


=head1 output

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
