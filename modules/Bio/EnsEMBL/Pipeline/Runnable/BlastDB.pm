
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BlastDB

=head1 SYNOPSIS

=head2 Methods:

=head1 DESCRIPTION

Indexing methods available:

wu_new - uses xdformat.
wu_old - uses setdb or pressdb, depending on you input sequence type.
ncbi   - uses formatdb

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BlastDB;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

=cut

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($sequences,
      $dbfile, 
      $type,
      $molecule_type,
      $workdir,
      $copy,
      $index_type,
      $make_fetchable_index) = $self->_rearrange([qw(SEQUENCES
						     DBFILE
						     TYPE
						     MOLECULE_TYPE
						     WORKDIR
						     COPY
						     INDEX_TYPE
						     MAKE_FETCHABLE_INDEX
						    )], @args);

  if (!defined($dbfile) && !defined($sequences)) {
    $self->throw("No dbfile or sequence array ref entered for indexing");
  }
  if (!defined($type))  { 
    $self->throw("Not database type entered");
  }

  $self->dbfile($dbfile)               if defined($dbfile);
  $self->sequences($sequences)         if defined($sequences);
  $self->molecule_type($type)          if defined($type);
  $self->molecule_type($molecule_type) if defined($molecule_type);
  $self->workdir($workdir)             if defined($workdir);
  $self->index_type($index_type);
  $self->make_fetchable_index($make_fetchable_index) 
    if defined($make_fetchable_index);

  $self->copy_db($dbfile) if (defined($dbfile) && $copy);

  return $self;
}

sub DESTROY {
  my $self = shift;

  $self->remove_index_files;
}

sub copy_db {
  my ($self, $dbfile) = @_;

  system("cp $dbfile " . $self->blastdb_dir);

  $dbfile =~ /([^\/]+)$/;

print "Database name is " . $1 . "\n";

  $self->throw("Didnt parse database names successfully.") unless $1;

  $self->dbfile($self->blastdb_dir."\/".$1);
  $self->copied_dbfile(1);
}

sub sequences {
  my ($self,$sequences) = @_;

  if (!defined($self->{_sequences})) {
    $self->{_sequences} = [];
  }
  if (defined($sequences)) {
    if (ref($sequences) eq "ARRAY") {
      push(@{$self->{_sequences}},@$sequences);
    } else {
      $self->throw("Argument must be an array ref . Currently is [$sequences]");
    }
  }
  return @{$self->{_sequences}};
}


sub run {
  my ($self) = @_;
  
  if (!defined($self->dbfile)) {
    my $seqfile = $self->make_seqfile;
    $self->dbfile($seqfile);
  }
  
  my $seqfile = $self->dbfile;
  
  if (! -e $seqfile) {
    $self->throw("Database [$seqfile] doesn't exist. Can't index");
  }
 
  my $command = $self->format_command . " " . $seqfile;
  print STDERR "Doing $command\n";
  my $exit_status = system($command);

  if ($exit_status) {
    return 0 
  } else {
    $self->db_formatted(1);
    return 1
  }
}


sub make_seqfile {
  my ($self) = @_;
  
  my $tmpdir = $self->blastdb_dir;
  
  my $blastfile = $self->get_tmp_file($tmpdir,'blast','fa');
  
  my $seqio = Bio::SeqIO->new(-format => 'Fasta',
			      -file   => ">$blastfile");
  
  foreach my $seq ($self->sequences) {
#    $seq->description('');
    $seqio->write_seq($seq);
  }
  
  close($seqio->_filehandle);
  
  return $blastfile;
}


sub dbname {
  my ($self) = @_;

  if (!$self->dbfile) {
    $self->throw("No database file defined.");
  }

  my $dbname = $self->dbfile;
  $dbname =~ /([^\/]+)$/;

  return $1;
}

sub remove_index_files {
  my ($self) = @_;

  #print "Removing index files created by $self.\n";
  
  if (!defined($self->dbfile)) {
    $self->throw("No database file defined - can't remove index files.");
  }

  opendir (BLASTDB_DIR, $self->blastdb_dir) 
    or $self->throw("Something has happened to the blast database"
		    ." directory.  Was it ever created?");

  my @dbfiles = readdir BLASTDB_DIR;
  
  closedir BLASTDB_DIR;

  foreach my $file (@dbfiles){
    next if $file =~ /^\./;
    unlink $self->blastdb_dir . "/" . $file;
  }
  
  if ($self->copied_dbfile) {
    unlink $self->dbfile;
  }

  rmdir $self->blastdb_dir
}

#################
# get/set methods 
#################

=head2 dbfile

 Title   : dbfile
 Usage   : $obj->dbfile($newval)
 Function: 
 Example : 
 Returns : value of dbfile
 Args    : newvalue (optional)


=cut

sub dbfile {
  my ($self,$value) = @_;

  if(defined $value) {
    $self->{_dbfile} = $value;
  }

  return $self->{_dbfile};  
}

=head2 db

 Title   : db
 Usage   : $obj->db
 Function: DEPRECATED
 Example : 
 Returns : DEPRECATED returns $self->dbfile
 Args    : none


=cut

sub db {
  my $self = shift;

  return $self->dbfile
}


=head2 molecule_type

 Title   : molecule_type
 Usage   : $self->molecule_type($newval)
 Function: Stores the sequence type - PROTEIN or DNA
 Example : 
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub molecule_type{
  my $self = shift;

  if(@_) {
    $self->{_molecule_type} = shift;
  }

  $self->warn("Molecule type not set.") 
    unless $self->{_molecule_type};

  return $self->{_molecule_type}
}

=head2 type

 Title   : type
 Usage   : $self->type($newval)
 Function: Maintains backwards compatibility.  Replaced by molecule_type.
 Example :
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub type {
  my $self = shift;

  if (@_) {
    my $type = shift;
    $self->molecule_type($type);
  }

  return $self->molecule_type;
}


sub db_formatted {
  my $self = shift;

  if (@_){
    $self->{_db_formatted} = shift
  }

  return $self->{_db_formatted}
}


sub copied_dbfile {
  my $self = shift;

  if (@_) {
    $self->{_copy} = shift;
  }

  return $self->{_copy}
}


sub index_type {
  my $self = shift;

  if (@_ || 
      (($self->{_index_type} !~ /ncbi/i)
       &&($self->{_index_type} !~ /wu/i))
     ) {

    my $value = shift if @_;

    # Default to new wu blast.
    $value = 'new wu' 
      unless $value;

    my $seqfetch_cmd_ok;

    if ($value =~ /wu/i) {

      if (($value =~ /new/i)&&($self->molecule_type eq 'DNA')) {
	$self->format_command('xdformat -n -I');
	$self->seqfetch_command(['xdget -n' , 'blastdb' , 'seqid']);
	$seqfetch_cmd_ok = 1;
	$self->{_index_type} = 'new_wu';
	
      } elsif (($value =~ /new/i)&&($self->molecule_type eq 'PROTEIN')) {
	$self->format_command('xdformat -p -I');
	$self->seqfetch_command(['xdget -p' , 'blastdb' , 'seqid']);
	$seqfetch_cmd_ok = 1;
	$self->{_index_type} = 'new_wu';
	
      } elsif ($value =~ /old/i) {
	
	$self->{_index_type} = 'old_wu';
	
	if ($self->molecule_type eq 'DNA') {
	  $self->format_command('pressdb');
### Can this be fixed?
	  $self->seqfetch_command('throw');
	  $seqfetch_cmd_ok = 1;

	} elsif ($self->molecule_type eq 'PROTEIN') {
	  $self->format_command('setdb');
### Can this be fixed?
	  $self->seqfetch_command('throw');
	  $seqfetch_cmd_ok = 1;
	}
      }
    } elsif ($value =~ /ncbi/i) {
      my $format_command = 'formatdb -o T ';;
      $format_command .= '-p F ' 
	if $self->molecule_type eq 'DNA';
      $format_command .= '-i ';

      $self->format_command($format_command);

      $self->seqfetch_command(['fastacmd -d' , 'blastdb' , '-s', 'seqid']);
      $seqfetch_cmd_ok = 1;

      $self->{_index_type} = 'ncbi';
    }

    $self->throw("Database indexing method [$value] not recognised.")
      unless ($self->format_command
	      && $seqfetch_cmd_ok);
  }

  return $self->{_index_type}
}

sub format_command {
  my $self = shift;

  if (@_){
    $self->{_format_command} = shift
  }

  $self->throw("Unable to return unset format command")
    unless $self->{_format_command};

  return $self->{_format_command}
}

sub seqfetch_command {
  my $self = shift;

  if(@_){
    $self->{_seqfetch_command} = shift;
    return
  }

  $self->throw("Seqfetch command not yet set.")
    unless $self->{_seqfetch_command};

  $self->throw("Blast database is not configured for sequence fetching (yet).") 
    if $self->{_seqfetch_command} eq 'throw';

  $self->throw("Blast database has not been set and/or formatted.\n"
	       . "You may need to use the \'run\' method on your object first.")
    unless ($self->dbfile && $self->db_formatted);

  my $command;

  foreach my $bit (@{$self->{_seqfetch_command}}){
    if ($bit eq 'blastdb'){
      $command .= $self->dbfile . ' ';
      next
    }
    if ($bit eq 'seqid') {
      # Unimplemented.  Could accept ids directly.
      next
    }
    
    $command .= $bit . ' ';
  }

  return $command
}


sub make_fetchable_index {
  my $self = shift;

  if (@_){
    $self->{make_fetchable_index} = shift;
  }

  return $self->{make_fetchable_index}
}

### Working directories ###

sub work_dir {
  my $self = shift;

  if (@_) {
    $self->{_work_dir} = shift;
    $self->{_work_dir} .= '/' 
      unless $self->{_work_dir} =~ /\/$/;
    return 1
  } else {
    $self->{_work_dir} = '/tmp/' if !defined $self->{_work_dir};
    return $self->{_work_dir}
  }
}

sub blastdb_dir {
  my $self = shift;

  unless ($self->{_blastdb_dir} && -d $self->{_blastdb_dir}) {

    my $blastdb_dir = $self->work_dir . "\/tempblast." . $$ . '/';

    mkdir $blastdb_dir;

    $self->throw("Failed to create a working directory for the blast database.") 
      unless (-d $blastdb_dir);

    $self->{_blastdb_dir} = $blastdb_dir;
  }

  return $self->{_blastdb_dir};
}

1;
