# Bio::EnsEMBL::Pipeline::IDFactory::Filename
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright (c) EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::IDFactory::Filename

=head1 SYNOPSIS

  ...
  my $idfactory = Bio::EnsEMBL::Pipeline::IDFactory::Filename->new(
                               -NAME => $name,
                               -CONFIG   => $config);

  while(my $id_set = $idfactory->next(20)) {
    #create some jobs with input ids $id_set
    ...
  }

=head1 DESCRIPTION

This is an implementation of the IDFactory interface.  It provides a means to
retrieve identifiers which are a list of filenames.  The filenames are
matched using a directory and regular expression from the configuration.

=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::IDFactory::Filename;

use Bio::EnsEMBL::Pipeline::IDFactory;
use Bio::EnsEMBL::Pipeline::IDSet;
use Cwd qw(cwd);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::IDFactory);



=head2 new

  Arg [-CONFIG]: Bio::EnsEMBL::Pipeline::Config
               The configuration for the running pipeline
  Arg [-NAME]: string
                The name of this ID Factory in the configuration file
  Example    : $idf = Bio::EnsEMBL::Pipeline::IDFactory::File->new(
                                           -NAME => $name,
                                           -CONFIG => $conf);
  Description: Instantiates an IDFactory class.  The configuration object
               that is passed to this constructor must define a directory to
               retrieve files from and a regular expression to be used to
               match the names of files.  This works in a similar way to the
               program 'find' in that all subdirectories are also searched
               recursively for matching filenames.
               If no directory is defined in the configuration the current
               working directory will be used instead.

               Filenames and directories beginning with '.' are ignored.
               The matching filenames are internally stored when the next
               method is called for the first time.  Any files added after this
                will not be included subsequent calls to next().

               ------
               [PEPSET_ID_FACTORY]
               module=Bio::EnsEMBL::Pipeline::IDFactory::Filename
               dir=/nfs/acari/mcvicker/data/
               regexp=^PeptideSet\.\d+$

               [BLASTP_TASK]
               id_factory = pepset_id_factory
               ------

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

  my $regexp = $self->config->get_parameter($name, 'regexp');

  if(!$regexp) {
    $self->throw("[$name] configuration must define a regular expression for" .
                 " the parameter regexp that can be used to match filenames");
  }

  $self->regexp($regexp);

  #use a specified directory if ther was one
  my $dir = $self->config->get_parameter($name, 'dir') || '';

  #if the dir starts with '/' it is fully qualified, otherwise it isn't
  if(!$dir || substr($dir,0,1) ne '/') {
    $dir = cwd() . $dir;
  }

  $self->dir($dir);

  return $self;
}


sub regexp {
  my $self = shift;
  $self->{'regexp'} = shift if(@_);
  return $self->{'regexp'};
}


sub dir {
  my $self = shift;
  $self->{'dir'} = shift if(@_);
  return $self->{'dir'};
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
  my $num = shift;

  my $files;
  #search for matching files the first time next is called
  if(!exists $self->{'files'}) {
    $self->{'files'} = $self->_matching_files($self->dir(), $self->regexp());
  }

  $files = $self->{'files'};

  my @ids;
  while(@$files && scalar(@ids) < $num) {
    push @ids, shift(@$files);
  }

  #if there are no ids left and some were requested, return undef instead
  return undef if(!@ids && (!defined($num) || $num > 0));

  return Bio::EnsEMBL::Pipeline::IDSet->new(-ID_LIST => \@ids);
}


sub _matching_files {
  my $self = shift;
  my $dir = shift;
  my $regexp = shift;

  #keep track of dirs that have been entered to avoid infinite recursion
  #when symlinks have created circular directory structure
  $self->{'seen_dir'} ||= {};
  return [] if($self->{'seen_dir'}->{$dir});
  $self->{'seen_dir'}->{$dir} = 1;

  if(!opendir(DIR, $dir)) {
    $self->warn("could not open dir '$dir'");
    return [];
  }

  my @files = grep !/^\./, readdir(DIR); #ignore files & dirs that start with .

  my @out;
  foreach my $file (@files) {
    #recurse through subdirectories
    if(-d $file) {
      push @out, @{$self->_matching_files("$dir/$file", $regexp)};
    }

    #take only files which match the regex
    if($file =~ /$regexp/) {
      push @out, "$dir/$file";
    }
  }

  closedir(DIR);

  return \@out;
}

1;
