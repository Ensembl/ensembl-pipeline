# Bio::EnsEMBL::Pipeline::IDFactory::File
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright (c) EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::IDFactory::File

=head1 SYNOPSIS

  ...
  my $idfactory = Bio::EnsEMBL::Pipeline::IDFactory::File->new(
                               -NAME => $name,
                               -CONFIG   => $config);

  while(my $id_set = $idfactory->next(20)) {
    #create some jobs with input ids $id_set
    ...
  }

=head1 DESCRIPTION

This is an implementation of the IDFactory interface.  It provides a means to 
retrieve identifiers stored in a set of files.  The files must be defined in 
the configuration that is passed to the constructor and must contain a single
id per line.

=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods have an underscore prefix.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::IDFactory::File;
 
use vars qw(@ISA);

use Cwd qw(cwd);
use IO::File;
use Bio::EnsEMBL::Pipeline::IDFactory;
use Bio::EnsEMBL::Pipeline::IDSet;

@ISA = qw(Bio::EnsEMBL::Pipeline::IDFactory);



=head2 new

  Arg [-CONFIG]: Bio::EnsEMBL::Pipeline::Config
               The configuration for the running pipeline
  Arg [-NAME]: string
                The name of this ID Factory in the configuration file
  Example    : $idf = Bio::EnsEMBL::Pipeline::IDFactory::File->new(
                                           -NAME => $name,
                                           -CONFIG => $conf);
  Description: Instantiates an IDFactory class.  The confuration object
               that is passed to this constructor must define a list of
               files under the heading of the id factory name provided.  The
               filenames should be ';' delimited and may be relative
               pathnames or fully qualified paths.  If the 'dir' parameter
               is specified  this value will be used as the base for relative 
               path names, otherwise the current working directory will be 
               used.  The files are expected to contain a single identifier 
               per line, all-whitespace lines are ignored.  An example of the 
               file and dir parameters in a configuration file are as follows:
               ------
               [FILE_ID_FACTORY]
               module=Bio::EnsEMBL::Pipeline::IDFactory::File
               files=file1;/nfs/acari/ensembl/data/file2;file3
               dir=/nfs/acari/mcvicker/data/

               [REPEAT_MASKER_TASK]
               id_factory = file_id_factory
               ------
               Files '1' and '3' will be looked for in the directory
               /nfs/acari/ensembl/data directory. File '2' will be looked
               for in the '/nfs/acari/ensembl/data/' directory.

               The files are not opened or verified for existance until the
               next() method is called for the first time.  All whitespace
               filenames are ignored. 

  Returntype : Bio::EnsEMBL::Pipeline::IDFactory
  Exceptions : thrown on incorrect arguments
  Caller     : none

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = bless({}, $class);

  #call superclass constructor
  $self = $self->SUPER::new(@_);

  my $name = $self->name();

  my $filenames = $self->config->get_parameter($name, 'files');


  #filter out any empty filenames
  my @files = grep {$_ =~ /\S/} split(/;/, $filenames);

  if(!@files) {
    $self->throw("[$name] configuration must define a value for" .
                 " the parameter 'files' in order to use the File IDFactory");
  }

  #use a specified directory if ther was one
  my $dir = $self->config->get_parameter($name, 'dir') || '';

  #if the dir starts with '/' it is fully qualified, otherwise it isn't
  if(!$dir || substr($dir,0,1) ne '/') {
    $dir = cwd() . $dir;
  }

  #make sure the directory ends with a '/'
  $dir = "$dir/" if(substr($dir, length($dir)-1, 1) ne '/');

  #construct a full path to each of the files
  foreach my $file (@files) {
    #if the file starts with a '/' it is already a full path
    if(substr($file, 0, 1) ne '/') {
      $file = join('',$dir,$file);
    }
  }

  $self->filenames(\@files);

  return $self;
}



=head2 filenames

  Arg [1]    : (optional) arrayref $filenames
  Example    : $self->filenames(\@filenames);
  Description: Getter/Setter for the fully qualified list of filenames used
               to obtain the inputs for this sequence fetcher
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : general

=cut

sub filenames {
  my $self = shift;

  $self->{'filenames'} = shift if(@_);

  return $self->{'filenames'};
}


=head2 next

  Arg [1]    : (optional) int $num_ids
               The number of ids to retrieve
  Example    : while($idset = $id_factory->next(20)) { do something; }; 
  Description: Retrieves an idset containing the next $num_ids from the files
               that this id factory is using. If there are less than 
               $num_ids left in the files, then the remaining ids will
               be returned as part of the idset.  If there are no
               ids left in the files undef will be returned instead.  If 
               the $num_ids argument was not defined, all remaining ids will
               be returned.
  Returntype : Bio::EnsEMBL::Pipeline::IDSet or undef
  Exceptions : none
  Caller     : Tasks

=cut

sub next {
  my $self = shift;
  my $num_ids = shift;

  my $cur_file = $self->{'current_file_num'} || 0;
  my @files = @{$self->filenames};

  #return undef if we are past last file in list of files
  if($cur_file == scalar(@files)) {
    undef;
  }

  my @fhs = @{$self->{'filehandles'} || []};

  #
  # If the files have not been opened yet, open them all at once.
  # This at least means that this will fail right away rather then
  # failing potentially much later if a single file is unreadable
  #
  if(!@fhs) {
    foreach my $file (@files) {
      my $fh = IO::File->new();
      $fh->open($file) or $self->throw("Could not open file $file: $!\n");
      push @fhs, $fh;
    }
    $self->{'filehandles'} = \@fhs;
  }

  my $fh = $fhs[$cur_file];

  #
  # Retrieve ids from files until we have retrieved the number requested
  # or we have run out of files
  #
  my @ids;
  while((!defined($num_ids) || (scalar(@ids) < $num_ids)) 
        && ($cur_file < scalar(@files))) {
    my $fh = $fhs[$cur_file];
    my $id = <$fh>;

    #if at EOF move onto the next file
    if(!defined($id)) {
      #close filehandles we no longer need to be open
      $fhs[$cur_file]->close;
      $fhs[$cur_file] = undef;
      $cur_file++;
      next;
    }

    chomp($id);

    #skip all whitespace lines
    next if($id =~ /^\s*$/);

    push(@ids, $id); 
  }

  #
  # save the current file so we can pick up where we left of next time 
  # this function is called
  #
  $self->{'current_file_num'} = $cur_file;

  #return undef if there were in fact no ids left
  return undef if((!defined($num_ids) || $num_ids != 0) && scalar(@ids) == 0);

  return Bio::EnsEMBL::Pipeline::IDSet->new(-ID_LIST => \@ids);
}


sub DESTROY {
  my $self = shift;

  #close any filehandles which are still open
  foreach my $fh (@{$self->{'filehandles'}}) {
    $fh->close() if($fh);
  }
}


1;
