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

Importer

=head1 SYNOPSIS

  my $dbobj = new Bio::EnsEMBL::Pipeline::DBSQL::Obj(-host => $host,
                                                     -dbname => $dbname,
                                                     -user   => 'ensadmin');

  my $mirror_dir = '/nfs/disk100/humpub/th/unfinished_ana/';

  my $importer   = new Importer(-dbobj      => $dbobj,
                                -mirror_dir => $mirror_dir);

  $importer->importClones;

=head1 DESCRIPTION

Reads sequences from files in a directory, makes them into clone objects
and writes them into the database.  Before writing it detects whether the database 
already has the sequence and deletes it if it is a previous version.  If the
sequence being written is the same or a previous version to that in the database 
a warning is printed and the sequence isn't written.


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Importer;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object;

use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::Object;

use Bio::EnsEMBL::PerlDB::Contig;
use Bio::EnsEMBL::PerlDB::Clone;
use Bio::EnsEMBL::Pipeline::DBSQL::Obj;

@ISA = qw(Bio::Root::Object);

=head2 new

    Title   :   new
    Usage   :   my $importer   = new Importer(-dbobj      => $dbobj,
                                              -mirror_dir => $mirror_dir);
    Function:   Initializes module
    Returns :   
    Args    :   Bio::EnsEMBL::Pipeline::DBSQL::Obj,directory name

=cut

sub _initialize {
  my ($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@_);    

  $self->{_dbobj}        = undef;     
  $self->{_mirror_dir}   = undef;
  $self->{_clones}       = [];

  my( $dbobj,$mirror_dir, ) = 
    $self->_rearrange(['DBOBJ','MIRROR_DIR'], @args);

  if ($dbobj) {
    $self->dbobj($dbobj);
  } else {
    $self->throw("No database object input");
  }
  if ($mirror_dir) {
    $self->mirror_dir($mirror_dir);
  } else {
    $self->throw("No mirror directory input");
  }

    return $self; # success - we hope!
}


sub importClones {
  my ($self) = @_;

  my @files   = $self->getNewFiles;
  my $logfile = $self->mirror_dir . "import.log";

  open(LOG,">>$logfile") || $self->throw("Can't write to logfile $logfile");

 FILE: foreach my $file (@files) {
      next FILE if ($file =~ /_cc/); 

      $self->{_clones} = [];
      $self->{_clonehash} = {};

      eval {
	my $clones = $self->readFile($file);
	
	my @keys = keys(%$clones);

	print ("\nNumber of clones is " . scalar(@keys) . "\n");
	
	$self->makeClones     ($clones);
	$self->checkClones    ($clones);
	$self->writeClones    ($clones);

      };
      if ($@) {
	$self->warn("ERROR: Error reading file  $file [$@]\n");
      } else {
	print LOG $file . "\n";
      }
  }
  close (LOG);

}

sub checkClones {
  my ($self) = @_;

  my @analysis = $self->dbobj->get_AnalysisAdaptor->fetch_by_logic_name('SubmitContig');

  if ($#analysis != 0) {
    $self->throw("More than one or none SubmitContig logic name. Eeek!");
  }

  my @newclones;

  foreach my $clone ($self->clones) {
    my $ok = 1;

      if (!defined($clone->id)) {
	$self->warn("ERROR: Clone $clone inhas no id");
	$ok = 0;
      }

      my $file = $self->{_clonehash}{$clone->id}{file};
      
      if (!defined($clone->htg_phase)) {
	$self->warn("ERROR: Clone " . $clone->id . " in file $file has no htg_phase");
	$ok = 0;
      }
      if (!defined($clone->embl_version)) {
	$self->warn("ERROR: Clone " . $clone->id . " in file $file has no embl_version");
	$ok = 0;
      }
      if (!defined($clone->embl_id)) {
	$self->warn("ERROR: Clone " . $clone->id . " in file $file has no embl_id");
	$ok = 0;
      }
      
      my @contigs = $clone->get_all_Contigs;
      
      if (scalar(@contigs) == 0) {
	$self->warn("ERROR: No contigs for clone " . $clone->id . " in file $file \n");
	$ok = 0;
      }

      foreach my $contig (@contigs) {
	if (!defined($contig->id)) {
	  $self->warn("ERROR: Contig $contig in file $file has no id");
	  $ok = 0;
	}
	if (!defined($contig->seq)) {
	  $self->warn("ERROR: Contig " . $contig->id ."  in file $file has no seq");
	  $ok = 0;
	}
	if (!defined($contig->embl_offset)) {
	  $self->warn("ERROR: Contig " . $contig->id ."  in file $file has no embl_offset");
	  $ok = 0;
	}
	if (!defined($contig->length)) {
	  $self->warn("ERROR: Contig " . $contig->id ."  in file $file has no length");
	  $ok = 0;
	}
	if (!defined($contig->version)) {
	  $self->warn("ERROR: Contig " . $contig->id ."  in file $file has no version");
	  $ok = 0;
	}
	if (!defined($contig->embl_version)) {
	  $self->warn("ERROR: Contig " . $contig->id ."  in file $file has no embl_version");
	  $ok = 0;
	}
	if (!defined($contig->embl_order)) {
	  $self->warn("ERROR: Contig " . $contig->id ."  in file $file has no embl_order");
	  $ok = 0;
	}
      }

    if ($ok == 1) {
      my $oldclone;

      eval {
	$oldclone = $self->dbobj->get_Clone($clone->id);
      };

      if ($@) {

	if ($clone->htg_phase == 0) {
	  $self->warn("PHASE: Ignoring phase 0 clone " . $clone->id . "\n");
	} else {
	  push(@newclones,$clone);
	}
      } else {
	if (defined($oldclone)) {
	  my $oldversion = $oldclone->embl_version;
	  
	  if ($oldversion > $clone->embl_version) {
	    $self->warn("ERROR : Inconsistent clone versions for " . $clone->id . " : old - $oldversion new - " . $clone->embl_version);
	  }elsif ($oldversion == $clone->embl_version) {
	    $self->warn("ERROR : Identical clone versions for " . $clone->id . " : old - $oldversion new - " . $clone->embl_version. " in file " . $self->{_clonehash}{$clone->id}{file} . "\n");
	  } else {
	    
	    print STDERR "Deleting clone : Found new version for " . $clone->id . " old - $oldversion new - " . $clone->embl_version . "\n";
	    
	    my $std = $self->dbobj->get_AnalysisAdaptor;
	    
	    foreach my $contig ($oldclone->get_all_Contigs) {
	      $std->removeInputId($contig->id,'contig',$analysis[0]);
	    }
	    
	    $self->dbobj->delete_Clone($clone->id);
	    push(@newclones,$clone);
	  }

	}
      }
    }
  }

  $self->{_clones} = [];

  foreach my $clone (@newclones) {
    $self->clones($clone);
  }
}

sub writeClones {
  my  ($self) = @_;

  print STDERR "Writing clones\n";

  my $std = $self->dbobj->get_AnalysisAdaptor;
  my @analysis = $std->fetch_by_logic_name('SubmitContig');

  if ($#analysis != 0) {
    $self->throw("More than one or none SubmitContig logic name. Eeek!");
  }
  print STDERR "Clones are " . $self->clones . "\n";
  foreach my $clone ($self->clones) {
    print STDERR "\nWriting clone ". $clone->id . "\n";

    eval {
      $self->dbobj->write_Clone($clone);
    };
    if (!$@) {
      foreach my $contig ($clone->get_all_Contigs) {
	$std->submitInputId($contig->id,'contig',$analysis[0]);
      }
    } else {
      $self->warn("Couldn't write clone " . $clone->id . " [$@]");
    }
  }
}

=head2 getNewFiles

    Title   :   getNewFiles
    Usage   :   @files = $obj->getNewFiles
    Function:   Gets files in the mirror directory that are new.
                i.e. their names aren\'t in the import.log file
    Returns :   @string
    Args    :   none

=cut

sub getNewFiles {
  my ($self) = @_;
  
  my $dir     = $self->mirror_dir;
  my $logfile = $self->logfile;

  open(IN,"<$logfile") || $self->throw("Can't read $logfile file");

  my %oldfiles;

  while (<IN>) {
    chomp;
    $oldfiles{$_} = 1;
  }

  opendir(DIR,"$dir/embl_files") || $self->throw("Can't read directory $dir");

  my @allfiles = readdir(DIR);

  closedir(DIR);

  my @newfiles;

  foreach my $file (@allfiles) {
#    print STDERR "Looking at file $file\n";
    if ($file !~ /^\./  && 
	! -d $file      &&
	$file =~ /\.mm/ &&
	$oldfiles{$file} != 1) {

      push(@newfiles,$file);
    }
  }

  print STDERR "There are " . scalar(@newfiles) . " new files\n";

  return @newfiles;
      
  
}

=head2 readFile

    Title   :   readFile
    Usage   :   @seqs = $obj->readFile($file)
    Function:   Reads sequences from files
    Returns :   @Bio::SeqI
    Args    :   @string

=cut

sub readFile {
  my ($self,$file) = @_;

  my $clonehash = $self->{_clonehash};

  my $dir = $self->mirror_dir;
  my $count = 1;
      
  print STDERR "\nReading fasta file " . $file . "\n\n";
  
  my $ccfile = $file;
  $ccfile =~ s/mm_u/mm_u_cc/;
  
  if (! -e "$dir/embl_files/$ccfile") {
    $self->throw("No contig coordinate file [$dir/embl_files/$ccfile]");
  }
  
  open(IN,"gunzip -c $dir/embl_files/$file   |") || $self>throw("Can't unzip file [$file]\n");
  open(CC,"gunzip -c $dir/embl_files/$ccfile |") || $self->throw("Can't unzip file [$ccfile]\n");
  
  my $seqio = new Bio::SeqIO(-fh     => \*IN,
			     -format => 'fasta');
  
  while (my $fasta = $seqio->next_seq) {
    printf STDERR  "Sequence %20s %s\n", $fasta->id,$fasta->desc;
    my $id = $fasta->id;
    $id =~ s/\.(.*)//;
    
    $clonehash->{$id}{seq}     = $fasta;
    $clonehash->{$id}{version} = $1;
    $clonehash->{$id}{file}    = $file;
    
    my $desc = $fasta->desc;
    my ($id,$type,$phase) = split(' ',$desc);
    $clonehash->{$id}{id} = $id;
    if ($type eq "ROD") {
      $clonehash->{$id}{phase} = 4;
    } else {
      $phase =~ s/HTGS_PHASE//;
      $clonehash->{$id}{phase} = $phase;
    }
  }
  
  close(IN);
  
  my $accession;
  
  while (<CC>) {
    chomp;
    if ($_ =~ /^SV\s+(.*)\.(.*)/) {
      $accession = $1;
      
      $clonehash->{$1}{version} = $2;
      $clonehash->{$1}{contigs} = [];
    } elsif ($_ =~ /(\d+)\s+(\d+): contig of/) {
      my %contighash;
      
      $contighash{start} = $1;
      $contighash{end}   = $2;
      
      push(@{$clonehash->{$accession}{contigs}},\%contighash);
      
    }
  }
  close(CC);

  return ($clonehash);
  
}

=head2 readChromosomes

    Title   :   readChromosomes
    Usage   :   $obj->readChromosomes
    Function:   Import the chromosomes from a mirrored file
    Returns :   nothing
    Args    :   none

=cut

sub readChromosomes {
  my ($self,$seqs) = @_;

  my $dir = $self->mirror_dir;

  open(CHR,"<$dir/chrlist/mouse_chr.lis") || $self->throw("Couldn't open chromosome file [$dir/chrlist/mouse_chr.lis");

  while (<CHR>) {
    chomp;
    
    my ($acc,$ver,$id,$chr,@dum) = split(' ',$_);
    
  }
  close(CHR);
}

=head2 makeClones

    Title   :   makeClones
    Usage   :   $obj->makeClones($clonehashref)
    Function:   Turns the data in $clonehashref into clone objects
    Returns :   Nothing
    Args    :   hash reference

=cut

sub makeClones {
  my ($self,$clones) = @_;
  
  my @acc = keys (%$clones);
  
  foreach my $acc (@acc) {
    
    if (defined($clones->{$acc}{seq})) {
      print "\nProcessing $acc\n\n";
      
      my $clone     = new Bio::EnsEMBL::PerlDB::Clone;      
      
      $self->clones($clone);

      my $seq   = $clones->{$acc}{seq};
      my $ver   = $clones->{$acc}{version};
      my $id    = $clones->{$acc}{id};
      my $chr   = $clones->{$acc}{chr};
      my $phase = $clones->{$acc}{phase};

      $clone->htg_phase   ($phase);
      $clone->id          ($acc);
      $clone->embl_version($ver);
      $clone->embl_id     ($id);
      $clone->version     (1);

      printf (STDERR "\tFound %10s seq\n",    $clones->{$acc}{seq}->length);
      printf (STDERR "\tFound %10s version\n",$clones->{$acc}{version});
      printf (STDERR "\tFound %10s id\n",     $clones->{$acc}{id});
      printf (STDERR "\tFound %10s chr\n",    $clones->{$acc}{chr});
      printf (STDERR "\tFound %10s phase\n",  $clones->{$acc}{phase});
      
      if ($phase != 4) {
	my @contigs = @{$clones->{$acc}{contigs}};
	
	my $count = 1;
	print STDERR "\nContigs : " . scalar(@contigs) . "\n";
	foreach my $contig (@contigs) {

	  my $seqstr    = $seq->subseq($contig->{start},$contig->{end});
	  my $offset    = $contig->{start};
	  
	  my $newcontig    = new Bio::EnsEMBL::PerlDB::Contig;
	  
	  # my $padlen    = 5 - length($count);
	  # my $pad       =  "0" x $padlen;
	  # my $contigid  = $acc . ".$pad$count";
	  my $contigid  = $acc . ".$ver.$offset." . $contig->{end};
	  
	  $newcontig->id          ($contigid);
	  $newcontig->seq         (new Bio::Seq(-id => $id, -seq =>$seqstr));    
	  $newcontig->embl_offset ($offset);
	  $newcontig->version     (1);
	  $newcontig->embl_version($ver);
	  $newcontig->embl_order  ($count);
	  $newcontig->chromosome  ($chr);

#	  print (STDERR "\tContig " . $newcontig->id . "\t : " . $newcontig->embl_offset . "\t" . ($newcontig->offset+$newcontig->length-1) . "\n");
	  
	  $clone->add_Contig   ($newcontig);
	  $count++;
	}
      } else {
	  my $newcontig    = new Bio::EnsEMBL::PerlDB::Contig;
	  
	  # my $contigid  = $acc . ".00001";
	  my $contigid  = $acc . ".1.1." . length $seq;


	  $newcontig->id          ($contigid);
	  $newcontig->seq         ($seq);    
	  $newcontig->embl_offset (1);
	  $newcontig->version     (1);
	  $newcontig->embl_version($ver);
	  $newcontig->chromosome  ($chr);
	  $newcontig->embl_order  (1);

#	  print (STDERR "\tContig " .$newcontig->id . "\t : " . $newcontig->embl_offset . "\t" . ($newcontig->offset+$newcontig->length-1) . "\n");

	  $clone->add_Contig   ($newcontig);

	}
    }
  }
}
  

=head2 dbobj

    Title   :   dbobj
    Usage   :   $db = $obj->dbobj
    Function:   Get/set method for database handle
    Returns :   Bio::EnsEMBL::Pipeline::DB::ObjI
    Args    :   Bio::EnsEMBL::Pipeline::DB::ObjI

=cut

sub dbobj {
  my ($self,$arg) = @_;
  

  if (defined($arg)) {
    if (!($arg->isa("Bio::EnsEMBL::Pipeline::DBSQL::Obj"))) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::DBSQL::Obj");
    }
    $self->{_dbobj} = $arg;
  }

  return $self->{_dbobj};
}


=head2 mirror_dir

    Title   :   mirror_dir
    Usage   :   $dir = $obj->mirror_dir
    Function:   Get/set method for mirror directory name
    Returns :   string
    Args    :   string

=cut

sub mirror_dir {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    if (! -d $arg) {
      $self->throw("Directory [$arg] doesn't exist");
    }
    $arg .= "/" if ($arg !~ /\/$/);
    $self->{_mirror_dir} = $arg;
  }

  return $self->{_mirror_dir};
}

=head2 clones

    Title   :   clones
    Usage   :   @clones = $obj->clones
    Function:   Get/set method for adding a clone and returning all clones
    Returns :   Bio::EnsEMBL::DB::CloneI
    Args    :   Bio::EnsEMBL::DB::CloneI

=cut

sub clones {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    if (!($arg->isa("Bio::EnsEMBL::DB::CloneI"))) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::DB::CloneI");
    }
    push(@{$self->{_clones}},$arg);
  }

  return @{$self->{_clones}};
}


sub logfile {
  my ($self) = @_;

  $self->throw("No mirror directory defined") unless defined($self->mirror_dir);

  return $self->mirror_dir . "import.log";
}

1;
