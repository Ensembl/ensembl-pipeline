
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

or something similar????

=head1 DESCRIPTION

Finds which pfam domains a particular swissprot hit matchs, finds the related hmms and runs genewiseHmm

=head1 CONTACT

rds@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::Finished_HalfwiseHMM;


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
use BlastableVersion;
#use DB_File;
use Fcntl;
use Data::Dumper;

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
    my ($genomic, $features, $hmmfetch, $hmmdb, $pfamdb, $memory, $options, $program) = $self->_rearrange([qw(QUERY
											   FEATURES
											   HMMFETCH
											   HMMDB 
											   PFAMDB
											   MEMORY
											   OPTIONS
											   PROGRAM)], @args);
    $self->throw("No genomic sequence input") unless defined($genomic);
    $self->query($genomic);
    $self->throw("No features input") unless defined($features);
    $self->add_input_features($features);
    $self->throw("No pfam database") unless $pfamdb;
    $self->pfamDB($pfamdb);
    $self->options($options);
    $self->program($program);
    
    if ($hmmfetch){
	$self->hmmfetch($hmmfetch);
    } else {
	$self->hmmfetch('/usr/local/ensembl/bin/hmmfetch');
    }
    
    if($hmmdb){
	$self->hmmdb($hmmdb);
    } else {
	$self->hmmdb('/data/blastdb/Ensembl/Pfam_ls');
	
    }
    $self->memory ($memory)    if (defined($memory));
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
sub pfam_db_version{
    my ($self) = @_;
    unless($self->{'_pfam_db_version'}){
        my $db = $self->pfamDB();
        $self->{'_pfam_db_version'} = $db->get_meta_value_by_key('version');
    }
    return $self->{'_pfam_db_version'};
}
sub pfam_ls_version{
    my ($self) = @_;
    unless($self->{'_pfam_ls_version'}){
	my $ver = BlastableVersion->new($self->hmmdb);
	$self->{'_pfam_ls_version'} = $ver->version();
    }
    return $self->{'_pfam_ls_version'};
}
sub program{
    my ($self, $program) = @_;
    $self->{'_program'} = $program if $program;
    return $self->{'_program'};
}
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
    my $swissprot_ids = $self->get_swissprot_ids();
    #print "got swiss prot ids \n";
    #my @keys = keys(%$swissprot_ids);
    
    #foreach my $key (@keys){
    #  print "hid = ".$key." strand = ".$swissprot_ids{$key}."\n";
    #}

    my $pfam_ids = $self->get_pfam_ids($swissprot_ids);
    #$self->throw(Dumper($pfam_ids));
    #print "got pfam ids \n";
    #$self->create_genewisehmm($pfam_ids, $dir);
    #$self->create_genewisehmm_individually($pfam_ids, $dir);
    $self->create_genewisehmm_complete($pfam_ids, $dir);
    #$self->throw(Dumper($pfam_ids));
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
    my $swissprot_ids; #{}
    
    foreach my $f(@features){
        my $hname = $f->hseqname();
        my $strand = $f->strand();
        #print "swissprot id = ".$hname."\n";
        if(!defined($swissprot_ids->{$hname})){
            $swissprot_ids->{$hname} = $strand;
        }else{
            $swissprot_ids->{$hname} = ($strand eq $swissprot_ids->{$hname} ? $strand : 0);
        }
      
    }
    return $swissprot_ids;
}

=head2 get_pfam_ids

    Arg      : reference to a hash of swissprotids and stands
    Function : gets all pfam ids for a particular swissprot id and puts them in a hash along with the strand the swissprot id is found on and returns its
    Exception: warns if swissprot id has no pfam domains
    Caller   : 

=cut

sub get_pfam_ids{
    my ($self, $swissprot_ids) = @_;
    my $pfam_accs;
    my $pfam_lookup;
    # CREATE TEMPORARY TABLE
    my $tbl_name = "pipeline_tmp_$$";
    my $create_table = qq{CREATE TEMPORARY TABLE $tbl_name(
                            pfamseq_id varchar(12) NOT NULL PRIMARY KEY,
                            strand enum('1','0','-1') DEFAULT '0'
                            )TYPE = HEAP
			}; # should this be HEAP?
                            # There's never gonna be that many matches
                            # to exceed tmp_table_size = 10048576 ???
    my $db = $self->pfamDB();
    my $sth = $db->prepare($create_table);
    $sth->execute();
    $sth->finish();
    # INSERT
    my (@binds, @values);
    my $sql = qq{INSERT IGNORE INTO $tbl_name (pfamseq_id, strand) VALUES };
    while (my ($swiss_id, $strand) = each(%$swissprot_ids)){
        #print STDERR "$swiss_id : $strand\n";
	push(@binds, $swiss_id, $strand);
	push(@values, qq{ (?, ?)});
    }
    if(scalar(@values)){
	$sql .= join(", ", @values);
	warn $sql;
	$sth = $db->prepare($sql);
	$sth->execute(@binds);
	$sth->finish();
    }
    # SELECT
#     my $select = qq{SELECT a.pfamA_acc, t.strand, t.pfamseq_id, a.pfamA_id, a.description
#                         FROM pfamA a, pfamA_reg_full f, pfamseq p, $tbl_name t
#                         WHERE p.auto_pfamseq = f.auto_pfamseq
#                         && a.auto_pfamA = f.auto_pfamA
#                         && p.pfamseq_id = t.pfamseq_id
#                         && f.in_full = 1};
    my $select = qq{SELECT a.pfamA_acc, t.strand, t.pfamseq_id, a.pfamA_id, p.description
			FROM pfamA_reg_full f, pfamseq p, $tbl_name t,  pfamA a
			WHERE f.auto_pfamseq = p.auto_pfamseq
			&& p.pfamseq_acc     = t.pfamseq_id
			&& f.in_full         = 1
			&& a.auto_pfamA      = f.auto_pfamA;};
    $sth = $db->prepare($select);
    $sth->execute();
    my ($pfam_acc, $strand, $swall, $pfam_id, $description);
    $sth->bind_columns(\($pfam_acc, $strand, $swall, $pfam_id, $description));
    while(my $row = $sth->fetchrow_arrayref()){
        print STDERR "$swall : $strand : $pfam_acc : $pfam_id\n";
        # making pfam_lookup for later ...
        $pfam_lookup->{$pfam_id} = [$pfam_acc, $description];
        if(defined($pfam_accs->{$pfam_acc})){
            $pfam_accs->{$pfam_acc} = ($strand eq $pfam_accs->{$pfam_acc} ? $strand : 0);
        }else{
            $pfam_accs->{$pfam_acc} = $strand;
        }
    }
    $self->pfam_lookup($pfam_lookup); # store the lookup somewhere useful.
    return $pfam_accs;
}

# get/set for lookup multi-valued hash { pfam_id => [pfam_acc, pfam_desc], ... }
# can append multiple to the lookup (see { %{$self->{'_pfam_lookup'}}, %{$hash_ref} })
sub pfam_lookup{
    my ($self, $hash_ref) = @_;
    if(ref($hash_ref) eq 'HASH'){
        $self->{'_pfam_lookup'} ||= {};
        $self->{'_pfam_lookup'}   = { %{$self->{'_pfam_lookup'}}, %{$hash_ref} };
    }
    return $self->{'_pfam_lookup'};
}

=head2 create_genewisehmm_individually

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

=head2 create_genewisehmm_individually

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
    $self->get_hmmdb($pfam_ids);
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

=head2 run_genewisehmm

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

=head2 get_GenewiseHMM

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

sub get_hmmdb{
    my ($self, $pfam_ids) = @_;
    my @pfamIds = keys(%$pfam_ids);
    print STDERR "getting the hmms for ids: ". join(" ", @pfamIds) . "\n"; ##########
    my $filename = substr(join("_", @pfamIds),0,20);
    $self->hmmfilename("$filename.$$.hmmdb");

    foreach my $id(@pfamIds){
        print "getting hmm for id ".$id."\n";
        my $command =  $self->hmmfetch." ".$self->hmmdb." ".$id." >> ".$self->hmmfilename;
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
