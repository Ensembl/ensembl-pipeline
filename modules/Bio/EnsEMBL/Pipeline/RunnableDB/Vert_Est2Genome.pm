#!/usr/local/bin/perl

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

Bio::EnsEMBL::Pipeline::RunnableDB::AlignFeature

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome->new(
					     -dbobj     => $db,
					     -input_id  => $id
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Vert_Est2Genome;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::Runnable::AlignFeature;
use Bio::EnsEMBL::Analysis::MSPcrunch;
use Bio::SeqIO;
use Bio::Tools::Blast;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::Object );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $dbobj, $db2,$input_id ) = $self->_rearrange(['DBOBJ',
						      'DB2',
						      'INPUT_ID'], @args);
       
    $self->throw("No database handle input")           unless defined($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::Pipeline::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::Pipeline::DB::ObjI");
    if (defined($db2)) {
    $self->warn("[$db2] is not a Bio::EnsEMBL::DB::ObjI") unless $db2->isa("Bio::EnsEMBL::DBSQL::Obj");
    }
    $self->dbobj($dbobj);
    $self->db2($db2);

    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);
    
    return $self; # success - we hope!
}
sub input_id {
	my ($self,$arg) = @_;

   if (defined($arg)) {
      $self->{_input_id} = $arg;
   }

   return $self->{_input_id};
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($db)
    Function:   Get/set method for database handle
    Returns :   Bio::EnsEMBL::Pipeline::DB::ObjI
    Args    :   

=cut 

sub dbobj {
    my( $self, $value ) = @_;    
    if ($value) {

        $value->isa("Bio::EnsEMBL::Pipeline::DB::ObjI") || $self->throw("Input [$value] isn't a Bio::EnsEMBL::Pipeline::DB::ObjI");
        $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
}
sub db2 {
    my ($self,$value) = @_;

    if (defined($value)) {
	$self->{_db2} = $value;
    }
    return $self->{_db2};
}
=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output
    Function:   Fetchs output data from a frozen perl object
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub fetch_output {
    my($self,$output) = @_;
    
    $output || $self->throw("No frozen object passed for the output");
    
    my $object;
    open (IN,"<$output") || do {print STDERR ("Could not open output data file... skipping job\n"); next;};
    
    while (<IN>) {
	$_ =~ s/\[//;
	$_ =~ s/\]//;
	$object .= $_;
    }
    close(IN);
    my @out;
   
    if (defined($object)) {
    my (@obj) = FreezeThaw::thaw($object);
    foreach my $array (@obj) {
	foreach my $object (@$array) {
	    push @out, $object;
	}
    }
    }
    return @out;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
    my($self,@features) = @_;

    my $db=$self->dbobj();

    @features || $self->throw("No frozen object passed for the output");
    my $contig;
    eval {
	$contig = $db->Contig($self->input_id);
    };
    if ($@) {
	print STDERR "Contig not found, skipping writing output to db\n";
    }
    else {
	my $feat_Obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
	$feat_Obj->write($contig,@features);
    }
    return 1;
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetchs input data for est2genome from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    

    $self->throw("No input id") unless defined($self->input_id);

    my $contigid  = $self->input_id;
    my $contig    = $self->db2->get_Contig($contigid);
    my $genseq   = $contig->primary_seq;
    my @features = $contig->get_all_SimilarityFeatures;
    $self->{_genseq} = $genseq;
    $self->{_features} = [];
    push(@{$self->{_features}},@features);
}

sub get_organism {
    my ($self,$hid) = @_;

    $self->throw("No hid input") unless defined($hid);
    my $org;
    eval {
	open (ORG,"efetch $hid |");

	while (<ORG>) {
	    chomp;
	    if (/^OS\s+(.*)/) {
		$org = $1;
		close(ORG);
		return $org;
	    }
	}
	close(ORG);
    }; 
    if ($@) {
	$self->warn("Efetch failed for $hid. Skipping. Error was [$@]\n");
    }
    return $org;
}

sub runnable {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI") unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");
	
	$self->{_runnable} = $arg;
    }

    return $self->{_runnable};
}

sub run {
    my ($self) = @_;

    my @mrnafeatures;
    my @features =@{$self->{_features}};
    my $genseq   = $self->{_genseq};
    my %idhash;

    foreach my $f (@features) {
	if (defined($f->analysis) && $f->analysis->db eq "vert"  && $f->score > 1000) {
	  my $organism = $self->get_organism($f->hseqname);
	  if (!defined($idhash{$f->hseqname})) { #
#	  if ($organism eq "Homo sapiens (human)" && !defined($idhash{$f->hseqname})) { #
		push(@mrnafeatures,$f);
		$idhash{$f->hseqname} =1;
	    } else {
	      print("Ignoring feature " . $f->seqname . "\n");
	    }
	}
    }
    print STDERR "Number of features pre blast is " . scalar(@mrnafeatures) . "\n";

    my @seq         = $self->get_Sequences(@mrnafeatures);
    my $blastdb     = $self->make_blast_db(@seq);
    my @newfeatures = $self->run_blast($genseq,$blastdb);

    print STDERR "Number of features post blast is " . scalar(@newfeatures) . "\n";

    my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::AlignFeature(-genomic  => $genseq,
								      -features => \@newfeatures);

    # Get rid of large data objects
    
    $self->{_features} = undef;

    $self->runnable($runnable);
    $self->throw("Can't run - no runnable object") unless defined($self->runnable);
    $self->runnable->minirun;

    my @tmpf = $self->runnable->output;
   
    foreach my $tmpf (@tmpf) {

      $tmpf->source_tag("tblastn");
      $tmpf->primary_tag("similarity");
      $tmpf->strand($tmpf->hstrand);

      print $tmpf->gff_string . "\n";
    }
    if ($#tmpf > 0) {
      my $i;
      for ($i = 0; $i <= $#tmpf-1; $i++) {
	$self->check_splice($tmpf[$i],$tmpf[$i+1]);
      }
    }
  }

sub check_splice {
  my ($self,$f1,$f2) = @_;


  my $splice1 = substr($self->{_genseq}->seq,$f1->end,2);
  my $splice2 = substr($self->{_genseq}->seq,$f2->start-3,2);

#  print ("Start/end " . $f1->start . "\t" . $f1->end . "\t" . $f2->start . "\t" . $f2->end . "\n");

    print ("Splices are " . $f1->hseqname . " [" . $splice1 . "][" . $splice2 . "] " . ($f2->start - $f1->end) . "\n");
}

sub output {
    my ($self) = @_;

    $self->throw("Can't return output  - no runnable object") unless defined($self->runnable);

    return $self->runnable->output;
}

sub get_Sequences {
    my ($self,@pairs) = @_;

    my @seq;

    foreach my $pair (@pairs) {
	my $id = $pair->hseqname;
	if ($pair->analysis->db eq "vert") {
	    my $seq = $self->get_Sequence($id);
	    push(@seq,$seq);
	}
    }
    return @seq;
}

sub make_blast_db {
    my ($self,@seq) = @_;

    my $blastfile = $self->get_tmp_file('/tmp/','blast','fa');
    my $seqio = Bio::SeqIO->new(-format => 'Fasta',
			       -file   => ">$blastfile");
    print STDERR "seq io is " . $seqio . "\n";
    print STDERR "Blast db file is $blastfile\n";
    foreach my $seq (@seq) {
	print STDERR "Writing seq " . $seq->id ."\n";
	$seqio->write_seq($seq);
    }

    close($seqio->_filehandle);

    my $status = system("pressdb $blastfile");
    print (STDERR "Status from pressdb $status\n");

    return $blastfile;
}

sub get_tmp_file {
    my ($self,$dir,$stub,$ext) = @_;

    
    if ($dir !~ /\/$/) {
	$dir = $dir . "/";
    }

    $self->check_disk_space($dir);

    my $num = int(rand(10000));
    my $file = $dir . $stub . "." . $num . "." . $ext;

    while (-e $file) {
	$num = int(rand(10000));
	$file = $stub . "." . $num . "." . $ext;
    }			
    
    return $file;
}

sub check_disk_space {
    my ($self,$dir,$minimumkb) = @_;

    $self->throw("No directory entered") unless defined($dir);

    open(DF,"df -k $dir |");

    my @lines = <DF>;
    $self->throw("Wrong number of lines output from df") unless scalar(@lines) == 2;
    my @f = split(' ',$lines[1]);

    my $kbytes = $f[3];

    if ($kbytes > $minimumkb) {
	return 1;
    } else {
	return 0;
    }
}

sub run_blast {
    my ($self,$seq,$db) = @_;

    my $blastout = $self->get_tmp_file("/tmp/","blast","tblastn_vert.msptmp");
    my $seqfile  = $self->get_tmp_file("/tmp/","seq","fa");

    my $seqio = Bio::SeqIO->new(-format => 'Fasta',
			       -file   => ">$seqfile");

    $seqio->write_seq($seq);
    close($seqio->_filehandle);

    my $command  = "blastn $db $seqfile B=500 -hspmax 1000 2> /dev/null |MSPcrunch -d - >  $blastout";

    print ("Running command $command\n");
    my $status = system($command );

    print("Exit status of blast is $status\n");


    my $msp = new Bio::EnsEMBL::Analysis::MSPcrunch(-file => $blastout,
						    -type => 'DNA-DNA',
						    -source_tag => 'vert_eg',
						    -contig_id => $self->input_id,
						    );

    unlink $blastout;
    unlink $seqfile;
    unlink $db;

    my @pairs = $msp->each_Homol;

    foreach my $pair (@pairs) {
	$self->print_FeaturePair($pair);
    }
    return @pairs;
}

sub print_FeaturePair {
    my ($self,$pair) = @_;

    print STDERR $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . $pair->score . "\t" .
	$pair->strand . "\t" . $pair->hseqname . "\t" . $pair->hstart . "\t" . $pair->hend . "\t" . $pair->hstrand . "\n";
}
=head2 parse_Header

  Title   : parse_Header
  Usage   : my $newid = $self->parse_Header($id);
  Function: Parses different sequence headers
  Returns : string
  Args    : none

=cut

sub parse_Header {
    my ($self,$id) = @_;

    if (!defined($id)) {
	$self->throw("No id input to parse_Header");
    }

    my $newid = $id;

    if ($id =~ /^(.*)\|(.*)\|(.*)/) {
	$newid = $2;
	$newid =~ s/(.*)\..*/$1/;
	
    } elsif ($id =~ /^..\:(.*)/) {
	$newid = $1;
    }
    $newid =~ s/ //g;
    return $newid;
}
    
=head2 get_Sequence

  Title   : get_Sequence
  Usage   : my $seq = get_Sequence($id);
  Function: Fetches all sequences with ids in the array
  Returns : ref to hash of Bio::PrimarySeq keyed by id
  Args    : none

=cut

sub get_Sequence {
    my ($self,$id) = @_;

    if (defined($self->{_seq_cache}{$id})) {
	return $self->{_seq_cache}{$id};
    } 

    my $newid = $self->parse_Header($id);

    next ID unless defined($newid);
    print(STDERR "New id :  is $newid [$id]\n");

    open(IN,"pfetch -q $newid |") || $self->throw("Error fetching sequence for id [$newid]");

    my $seq;
	
    while (<IN>) {
	chomp;
	$seq .= $_;
    }
	
    if (!defined($seq) || $seq eq "no match") {
	open(IN,"efetch -q $newid |") || $self->throw("Error fetching sequence for id [$newid]");
           
	$seq = "";
	    
	while (<IN>) {
	    chomp;
	    $seq .= $_;
	}
    }

    if (!defined($seq)) {
	$self->throw("Couldn't find sequence for $newid [$id]");
    }
    
    my $seq = new Bio::Seq(-id  => $newid,
			   -seq => $seq);
	
    $self->{_seq_cache}{$id} = $seq;

    return $seq;

}


1;


