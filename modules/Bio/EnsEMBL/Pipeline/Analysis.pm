#
# Object for storing sequence analysis details
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

Bio::EnsEMBL::Pipeline::Analysis - Stores details of an analysis run

=head1 SYNOPSIS

 my $obj    = new Bio::EnsEMBL::Analysis::Analysis
    ('-id'              => $id,
     '-logical_name'    => 'SWIRBlast',
     '-db'              => $db,
     '-db_version'      => $db_version,
     '-db_file'         => $db_file,
     '-program'         => $program,
     '-program_version' => $program_version,
     '-program_file'    => $program_file,
     '-gff_source'      => $gff_source,
     '-gff_feature'     => $gff_feature,
     '-module'          => $module,
     '-module_version'  => $module_version,
     '-parameters'      => $parameters,
     '-created'         => $created
     );

=head1 DESCRIPTION

Object to store details of an analysis run.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Analysis;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::AnalysisI;
use Bio::Root::RootI;

# Inherits from the base bioperl object
@ISA = qw(Bio::EnsEMBL::AnalysisI Bio::Root::RootI);

=head2 new

 Title   : new
 Usage   : my $analysis = new Bio::EnsEMBL::Pipeline::Analysis(%params);
 Function: Creates a new Bio::EnsEMBL::Pipeline::Analysis object
 Returns : Bio::EnsEMBL::Pipeline::Analysis
 Args    : -id:              Analysis ID
           -db:              DB Name
           -db_version:      DB Version
           -db_file:         DB File
           -program:         Program name
           -program_version: Program version
           -program_file:    Program file
           -gff_source:      GFF source
           -gff_feature:     GFF feature
           -module:          Module 
           -module_version:  Module version
           -parameters:      Parameters used
           -created:         Date created
           -logic_name:      Logical name

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($id,$db,$db_version,$db_file,$program,$program_version,$program_file,
      $gff_source,$gff_feature,$module,$module_version,$parameters,$created,
      $logic_name ) =
	  $self->_rearrange([qw(ID
				DB
				DB_VERSION
				DB_FILE
				PROGRAM
				PROGRAM_VERSION
				PROGRAM_FILE
				GFF_SOURCE
				GFF_FEATURE
				MODULE
				MODULE_VERSION
				PARAMETERS
				CREATED
				LOGIC_NAME
				)],@args);

  $self->dbID           ($id);
  $self->db             ($db);
  $self->db_version     ($db_version);
  $self->db_file        ($db_file);
  $self->program        ($program);
  $self->program_version($program_version);
  $self->program_file   ($program_file);
  $self->module         ($module);
  $self->module_version ($module_version);
  $self->gff_source     ($gff_source);
  $self->gff_feature    ($gff_feature);
  $self->parameters     ($parameters);
  $self->created        ($created);
  $self->logic_name     ( $logic_name );

  return $self; 
}


=head2 dbID

  Title   : dbID
  Usage   : $self->dbID
  Function: Get/set method for the id
  Returns : int
  Args    : int

=cut

sub dbID {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_dbid'} = $arg;
    }
    return $self->{'_dbid'};
}


=head2 db

  Title   : db
  Usage   : $self->db
  Function: Get/set method for database
  Returns : String
  Args    : String

=cut

sub db {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_db'} = $arg;
    }

    return $self->{'_db'};
}


=head2 db_version

  Title   : db_version
  Usage   : $self->db_version
  Function: Get/set method for the database version number
  Returns : int
  Args    : int

=cut

sub db_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_db_version'} = $arg;
    }

    return $self->{'_db_version'};
}


=head2 db_file

  Title   : db_file
  Usage   : $self->db_file
  Function: Get/set method for the database file
  Returns : string
  Args    : string

=cut

sub db_file {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_db_file'} = $arg;
    }

    return $self->{'_db_file'};
}


=head2 program

  Title   : program
  Usage   : $self->program
  Function: Get/set method for the program name
  Returns : String
  Args    : String

=cut

sub program {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_program'} = $arg;
    }

    return $self->{'_program'};
}


=head2 program_version

  Title   : program_version
  Usage   : $self->program_version
  Function: Get/set method for the program version number
  Returns : int
  Args    : int

=cut

sub program_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_program_version'} = $arg;
    }

    return $self->{'_program_version'};
}

=head2 program_file

  Title   : program_file
  Usage   : $self->program_file
  Function: Get/set method for the program file
  Returns : string
  Args    : string

=cut

sub program_file {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_program_file'} = $arg;
    }

    return $self->{'_program_file'};
}


=head2 module

  Title   : module
  Usage   : $self->module
  Function: Get/set method for the module name
  Returns : String
  Args    : String

=cut

sub module {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_module'} = $arg;
    }

    return $self->{'_module'};
}


=head2 module_version

  Title   : module_version
  Usage   : $self->module_version
  Function: Get/set method for the module version number
  Returns : string
  Args    : string

=cut

sub module_version {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_module_version'} = $arg;
    }

    return $self->{'_module_version'};
}

=head2 gff_source

  Title   : gff_source
  Usage   : $self->gff_source
  Function: Get/set method for the gff_source tag
  Returns : String
  Args    : String

=cut

sub gff_source {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_gff_source'} = $arg;
    }

    return $self->{'_gff_source'};
}

=head2 gff_feature

  Title   : gff_feature
  Usage   : $self->gff_feature
  Function: Get/set method for the gff_feature tag
  Returns : String
  Args    : String

=cut

sub gff_feature {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_gff_feature'} = $arg;
    }

    return $self->{'_gff_feature'};
}

=head2 parameters

  Title   : parameters
  Usage   : $self->parameters
  Function: Get/set method for the parameter string
  Returns : String
  Args    : String

=cut

sub parameters {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_parameters'} = $arg;
    }

    return $self->{'_parameters'};
}

=head2 created

  Title   : created
  Usage   : $self->created
  Function: Get/set method for the created time
  Returns : String
  Args    : String

=cut

sub created {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{'_created'} = $arg;
    }

    return $self->{'_created'};
}

=head2 logic_name

  Title   : logic_name
  Usage   : $self->logic_name
  Function: Get/set method for the logic_name, the name under
            which this typical analysis is known.
  Returns : String
  Args    : String

=cut


sub logic_name {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ($self->{'_logic_name'} = $arg);
  $self->{'_logic_name'};
}

1;
















