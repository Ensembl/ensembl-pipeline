#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
# Last modified by SCP 15/06/2001
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Importer

=head1 SYNOPSIS

  my $dbobj = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(
      '-host'   => $host,
      '-dbname' => $dbname,
      '-user'   => 'ensadmin'
  );

  my $mirror_dir = '/nfs/disk100/humpub/th/unfinished_ana/';

  my $importer   = new Importer(
      '-dbobj'      => $dbobj,
      '-mirror_dir' => $mirror_dir
  );

  $importer->retain_old_versions(1);
  $importer->importClones;

=head1 DESCRIPTION

Reads sequences from files in a directory, makes them into clone objects
and writes them into the database. Before writing it looks to see what
previous versions of that clone exist. If it is a previous version and
the option retain_old_versions is set, it will keep the existing version,
otherwise it will delete it. If the sequence being written is the same or
a previous version to that in the database a warning is printed and the
sequence isn't written.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut
#' # make emacs happy


package Bio::EnsEMBL::Pipeline::Importer;

use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::EnsEMBL::Root;

use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;

use Bio::EnsEMBL::RawContig;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

use FileHandle;

@ISA = qw(Bio::EnsEMBL::Root);

=head2 new

    Title   :   new
    Usage   :   my $importer   = new Importer(
		    '-dbobj'      => $dbobj,
                    '-mirror_dir' => $mirror_dir,
                    '-species'    => $species
		);
    Function:   Initializes module
    Returns :   
    Args    :   Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor,
                directory name,
                species

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    

    $self->{'_dbobj'}        = undef;     
    $self->{'_mirror_dir'}   = undef;
    $self->{'_clones'}       = [];
    $self->{'_chromosome'}   = undef;
    $self->{'_retain_old_versions'} = undef;

    my( $dbobj,$mirror_dir,$species) = 
	$self->_rearrange([qw(DBOBJ MIRROR_DIR SPECIES)], @args);

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
    if ($species) {
	$self->species($species);
    } else {
	$self->throw("No species input");
    }

    return $self;
}


sub importClones {
  my ($self) = @_;

  my @files;
  # import from a specific list of files
  if (defined $self->files) {
    my $tmp = $self->files;
    foreach my $file (@{$tmp}) {
      if (-e ($self->mirror_dir . '/embl_files/' . $file)) {
        push @files, $file;
      }
      else {
        print STDERR "Can't find file $file\n";
      }
    }
    $self->throw("Couldn't find any of the specified files") 
     unless scalar @files > 0;
  }
  else {
    @files = $self->getNewFiles;
  }
  my $logfile = $self->logfile;

  open(LOG,">>$logfile") || $self->throw("Can't write to logfile $logfile");
  autoflush LOG;

  # $self->readChromosomes;

  FILE: foreach my $file (@files) {
    next FILE if ($file =~ /_cc/); 

    $self->{'_clones'} = [];
    $self->{'_clonehash'} = {};

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

  my $sic = $self->dbobj->get_StateInfoContainer;

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

      my $file = $self->{'_clonehash'}{$clone->id}{'file'};
      
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
	    $self->warn("ERROR : Identical clone versions for " . $clone->id . " : old - $oldversion new - " . $clone->embl_version. " in file " . $self->{'_clonehash'}{$clone->id}{'file'} . "\n");
	  } else {
            if ($self->retain_old_versions) {
	      print STDERR "Found new version for " . $clone->id . " old - $oldversion new - " . $clone->embl_version . "; keeping old version\n";
            }
            else
            {
	      print STDERR "Deleting clone : Found new version for " . $clone->id . " old - $oldversion new - " . $clone->embl_version . "\n";
	    
	      foreach my $contig ($oldclone->get_all_Contigs) {
	        $sic->delete_input_id($contig->id,'contig',$analysis[0]);
	      }
	    
	      $oldclone->delete;
            }
	    push(@newclones,$clone);
	  }

	}
      }
    }
  }

  $self->{'_clones'} = [];

  foreach my $clone (@newclones) {
    $self->clones($clone);
  }
}

sub writeClones {
  my  ($self) = @_;

  print STDERR "Writing clones\n";

  my @analysis = $self->dbobj->get_AnalysisAdaptor->fetch_by_logic_name('SubmitContig');

  if ($#analysis != 0) {
    $self->throw("More than one or none SubmitContig logic name. Eeek!");
  }
  print STDERR "Clones are " . $self->clones . "\n";
  foreach my $clone ($self->clones) {
    print STDERR "\nWriting clone ". $clone->id . "\n";

    eval {
      $self->dbobj->get_CloneAdaptor->store($clone);
    };
#    if (!$@) {
#      foreach my $contig ($clone->get_all_Contigs) {
#	$sic->store_input_id_class_analysis($contig->id,'contig',$analysis[0]);
#      }
    if($@) {
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
  my $fileext;
  
  my $dir     = $self->mirror_dir;
  my $logfile = $self->logfile;

  $self->throw("Need to specify species: mm or hs")
   unless ($self->species eq 'mm' || $self->species eq 'hs');

  open(IN,"<$logfile") || $self->throw("Can't read $logfile file");

  my %oldfiles;

  while (<IN>) {
    chomp;
    $oldfiles{$_} = 1;
  }

  opendir(DIR,"$dir/embl_files") || $self->throw("Can't read directory $dir");

  my @allfiles = sort {$b cmp $a} readdir(DIR); # newest first

  closedir(DIR);

  my @newfiles;
  my $spec = $self->species;

  foreach my $file (@allfiles) {
#    print STDERR "Looking at file $file\n";
    if ($file !~ /^\./  && 
	! -d $file      &&
	$file =~ /\.$spec/ &&
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

  my $clonehash = $self->{'_clonehash'};

  my $dir = $self->mirror_dir;
  my $count = 1;
  my $finished = 0;
  my $ccfile;
  my ($spec,$ftype);
  my (%idlist);

  if (defined $self->ids) {
    my $file = $self->ids;
    open FILE, "< $file";
    while (<FILE>) {
      chomp;
      $idlist{$_} = 1;
    }
  }

  $self->throw("Need to specify species: mm or hs")
   unless ($self->species eq 'mm' || $self->species eq 'hs');
  if ($self->species eq 'mm') {
    $ftype = 'ROD';
  }
  else {
    $ftype = 'HUM';
  }
      
  print STDERR "\nReading fasta file " . $file . "\n\n";

  if ($file =~ /_f/) {
    $finished = 1;
  }
  else {
    $ccfile = $file;
    $spec = $self->species;
    $ccfile =~ s/$spec\_u/$spec\_u_cc/;
    if (! -e "$dir/embl_files/$ccfile") {
      $self->throw("No contig coordinate file [$dir/embl_files/$ccfile]");
    }
    open(CC,"gunzip -c $dir/embl_files/$ccfile |") || $self->throw("Can't unzip file [$ccfile]\n");
  }
  
  open(IN,"gunzip -c $dir/embl_files/$file   |") || $self->throw("Can't unzip file [$file]\n");
  
  my $seqio = new Bio::SeqIO(-fh     => \*IN,
			     -format => 'fasta');
  
  while (my $fasta = $seqio->next_seq) {
    my $id = $fasta->id;
    next if defined $self->ids && (!defined $idlist{$id});
    printf STDERR  "Sequence %20s %s\n", $id, $fasta->desc;
    if ($finished && $fasta->desc =~ /HTGS_PHASE/) {
      $self->warn("Found unfinished sequence amongst allegedly finished sequence; skipping "
      . $fasta->desc);
      next;
    }
    $id =~ s/\.(.*)//;
    
    $clonehash->{$id}{'seq'}     = $fasta;
    $clonehash->{$id}{'version'} = $1;
    $clonehash->{$id}{'file'}    = $file;
    my $accver = "$id." . $clonehash->{$id}{'version'};
    # $clonehash->{$id}{chr}     = $self->{'_chromosome'}{$accver}
    #  if (defined $self->{'_chromosome'}{$accver});
    
    my $desc = $fasta->desc;
    my ($tmp,$type,$phase) = split(' ',$desc);
    $clonehash->{$id}{'id'} = $id;
    if ($type eq $ftype) {
      $clonehash->{$id}{'phase'} = 4;
    } else {
      $phase =~ s/HTGS_PHASE//;
      $clonehash->{$id}{'phase'} = $phase;
    }
  }
  
  close(IN);
  
  my $accession;
  
  if (! $finished) {
    while (<CC>) {
      chomp;
      if ($_ =~ /^SV\s+(.*)\.(.*)/) {
	$accession = $1;
      
	$clonehash->{$accession}{'version'} = $2;
	$clonehash->{$accession}{'contigs'} = [];
      } elsif ($_ =~ /(\d+)\s+(\d+):? contig of/) {
      
        push(@{$clonehash->{$accession}{'contigs'}},{ 'start' => $1, 
						      'end'   => $2});
      
      } elsif ($_ =~ /Contig (\d+):\s+(\d+)\-(\d+)/) {
        push(@{$clonehash->{$accession}{'contigs'}},{ 'start' => $2,
						      'end'   => $3} );
      }
    }
    close(CC);
  }

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
  my ($self) = @_;

  my $dir = $self->mirror_dir;
  my $chr_file;

  if ($self->species eq 'mm') {
    $chr_file = 'mouse_chr.lis';
  }
  else {
    $chr_file = 'human_chr.lis';
  }

  open(CHR,"<$dir/chrlist/$chr_file") || $self->throw("Couldn't open chromosome file [$dir/chrlist/$chr_file");

  while (<CHR>) {
    chomp;
    
    my ($acc,$ver,$id,$chr,@dum) = split(' ',$_);
    my $accver = "$acc.$ver";
    $self->{'_chromosome'}{$accver} = $chr;
    
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
    
    if (defined($clones->{$acc}{'seq'})) {
      print "\nProcessing $acc\n\n";
      
      my $clone     = new Bio::EnsEMBL::Clone;      
      
      $self->clones($clone);

      my $seq   = $clones->{$acc}{'seq'};
      my $ver   = $clones->{$acc}{'version'};
      my $id    = $clones->{$acc}{'id'};
      my $chr   = $clones->{$acc}{'chr'};
      my $phase = $clones->{$acc}{'phase'};

      $clone->htg_phase   ($phase);
      $clone->id          ($acc);
      $clone->embl_version($ver);
      $clone->embl_id     ($id);
      $clone->version     (1);

      printf (STDERR "\tFound %10s seq\n",    $clones->{$acc}{'seq'}->length);
      printf (STDERR "\tFound %10s version\n",$clones->{$acc}{'version'});
      printf (STDERR "\tFound %10s id\n",     $clones->{$acc}{'id'});
      printf (STDERR "\tFound %10s chr\n",    $clones->{$acc}{'chr'});
      printf (STDERR "\tFound %10s phase\n",  $clones->{$acc}{'phase'});
      
      if ($phase != 4) {
	my @contigs = @{$clones->{$acc}{'contigs'}};

	# if we haven't got a list of clones by reading the cc file
	# need to divide clone into contigs the dirty way ...
	if (scalar @contigs == 0) {
	  my @limits = $self->scanClone($seq->seq);
	  foreach my $l (@limits) {
	    push @contigs, {
	      'start' => $l->[0],
	      'end'   => $l->[1]
	    };
	  }
	}
	
	my $count = 1;
	print STDERR "\nContigs : " . scalar(@contigs) . "\n";
	foreach my $contig (@contigs) {
	  my $seqstr    = $seq->subseq($contig->{'start'},$contig->{'end'});
	  my $offset    = $contig->{'start'};
	  
	  my $newcontig    = new Bio::EnsEMBL::RawContig;
	  
	  my $contigid  = "$acc.$ver.$offset." . $contig->{'end'};
	  
	  $newcontig->name          ($contigid);
	  $newcontig->seq         ($seqstr);
	  $newcontig->embl_offset ($offset);
	  # $newcontig->chromosome  ($chr);

#	  print (STDERR "\tContig " . $newcontig->id . "\t : " . $newcontig->embl_offset . "\t" . ($newcontig->offset+$newcontig->length-1) . "\n");
	  
	  $clone->add_Contig   ($newcontig);
	  $count++;
	}
      } else {
	  my $newcontig = new Bio::EnsEMBL::RawContig;

	  my $contigid  = "$acc.$ver.1." . $seq->length;

	  $newcontig->name        ($contigid);
	  $newcontig->seq         ($seq->seq);    
	  $newcontig->embl_offset (1);
	  # $newcontig->chromosome  ($chr);

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
    if (!($arg->isa("Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor"))) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{'_dbobj'} = $arg;
  }

  return $self->{'_dbobj'};
}


=head2 species

    Title   :   species
    Usage   :   $species = $obj->species
    Function:   Get/set method for species (mm/hs)
    Returns :   string
    Args    :   string

=cut

sub species {
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
    if ($arg ne 'hs' && $arg ne 'mm') {
      $self->throw("Only mm [mouse] and hs [human] options supported");
    }
    $self->{'_species'} = $arg;
  }

  return $self->{'_species'};
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
    $self->{'_mirror_dir'} = $arg;
  }

  return $self->{'_mirror_dir'};
}

=head2 clones

    Title   :   clones
    Usage   :   @clones = $obj->clones
    Function:   Get/set method for adding a clone and returning all clones
    Returns :   Bio::EnsEMBL::Clone
    Args    :   Bio::EnsEMBL::Clone

=cut

sub clones {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    push(@{$self->{'_clones'}},$arg);
  }

  return @{$self->{'_clones'}};
}


sub logfile {
  my ($self) = @_;

  $self->throw("No mirror directory defined") unless defined($self->mirror_dir);
  $self->throw("No species defined") unless defined($self->species);

  return $self->mirror_dir . "import_" . $self->species . ".log";
}

=head2 retain_old_versions

    Title   :   retain_old_versions
    Usage   :   $dir = $obj->retain_old_versions
    Function:   Whether to delete or retain old versions of clones
		true/flase to delete/retain previous versions
    Returns :   boolean
    Args    :   boolean

=cut

sub retain_old_versions {
  my ($self,$arg) = @_;

  if (defined $arg) {
    if ($arg) {
      $self->{'_retain_old_versions'} = 1;
    }
    else {
      $self->{'_retain_old_versions'} = 0;
    }
  }
  return $self->{'_retain_old_versions'};
}


=head2 scanClone

    Title   :   scanClone
    Usage   :   @contigs = $obj->scanClone($seq)
    Function:   Scans the clone sequence to find positions of contigs
		by assuming at least x n's  between contigs
    Returns :   contig positions: list of lists (start, end)
    Args    :   seq: string

=cut
#'

sub scanClone {
  my($self, $seq) = @_;
  my(@gaps, @contig);
  my($start, $gap);

  # get a list of gaps - at least 50 bp
  my $pos = 0;
  while ($pos < length $seq) {
    my $unused = substr $seq, $pos;
    ($gap) = $unused =~ /(n{50,})/;
    last unless $gap;
    $start = 1 + index $seq, $gap, $pos;
    push @gaps, [ $start, $start + length($gap) - 1 ];
    $pos = $start + length $gap;
  }

  # calc coords of contigs

  if (@gaps){
    # 1st contig before 1st gap unless the sequence starts off with a gap
    push @contig, [1, $gaps[0]->[0] - 1] unless $gaps[0]->[0] == 1;

    # contigs other than 1st and last are between gaps
    foreach my $i (0 .. $#gaps - 1) {
      push @contig, [$gaps[$i]->[1] + 1, $gaps[$i + 1]->[0] - 1];
    }

    # last contig after last gap unless the sequence ends with a gap
    push @contig, [$gaps[$#gaps]->[1] + 1, length($seq)]
     unless $gaps[$#gaps]->[1] == length($seq);
  }
  else {
    # no gaps
    push @contig, [1, length($seq)];
  }

  return @contig;
}

=head2 files

    Title   :   files
    Usage   :   $obj->files([...])
    Function:   Get/set for a specific list of files to read
    Returns :   array ref
    Args    :   array ref

=cut

sub files {
  my($self, $files) = @_;

  if (defined $files) {
    $self->{'_files'} = $files;
  }
  return $self->{'_files'};
}

=head2 ids

    Title   :   ids
    Usage   :   $obj->ids('ids.txt')
    Function:   Get/set for a specific list of ids
    Returns :   filename
    Args    :   filename

=cut

sub ids {
  my($self, $ids) = @_;

  if (defined $ids) {
    $self->{'_ids'} = $ids;
  }
  return $self->{'_ids'};
}

1;
