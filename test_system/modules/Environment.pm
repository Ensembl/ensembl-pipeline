package Environment;

use lib './config';
use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Root;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);





sub new{
  my ($class) = @_;
  my $self = bless {}, $class;
  return $self;
}


#containers
sub old_perl5lib{
  my $self = shift;
  $self->{'old_perl5lib'} = shift if(@_);
  return $self->{'old_perl5lib'};
}

sub old_blastdb{
  my $self = shift;
  $self->{'old_blastdb'} = shift if(@_);
  return $self->{'old_blastdb'};
}



sub add_to_perl5lib{
  my ($self, $addition) = @_;
  my $perl5lib = $ENV{'PERL5LIB'};
  $self->old_perl5lib($perl5lib);
  $addition .= ':'.$perl5lib;
  $ENV{'PERL5LIB'} = $addition;
  return $addition;
}

sub change_blastdb{
  my ($self, $blastdb) = @_;
  $self->old_blastdb($ENV{'BLASTDB'});
  $ENV{'BLASTDB'} = $blastdb;
  return $blastdb;
}


sub setup_paths{
  my ($self, $logic_name, $species) = @_;
  my $method = "setup_".$logic_name."_paths";
  if($self->can($method)){
    $self->$method($species);
  }
}


sub setup_BestPmatch_paths{
  my ($self, $species) = @_;
  $self->setup_Pmatch_paths($species);
}

sub setup_Pmatch_paths{
  my ($self, $species) = @_;
  my $pwd = $ENV{PWD};
  my $pmatch_fasta_path = $pwd."/".$species."/data/pmatch_proteins.fa";
  if(! -e $pmatch_fasta_path){
    throw("Can't setup environment for pmatch config ".$pmatch_fasta_path.
          " doesn't exist");
  }else{
    $ENV{PMATCH_FASTA} = $pmatch_fasta_path;
  }
}

sub setup_all_paths{
  my ($self, $species) = @_;
  $self->setup_Pmatch_paths($species);
}

sub return_environment{
  my ($self) = @_;

  $ENV{'PERL5LIB'} = $self->old_perl5lib;
  $ENV{'BLASTDB'} = $self->old_blastdb;
}


sub DESTROY{
  my ($self) = @_;
  $self->return_environment;
}

1;
