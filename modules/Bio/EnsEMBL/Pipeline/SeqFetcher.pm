#!/usr/local/bin/perl

#
#
# Cared for by Val Curwen  <vac@sanger.ac.uk>
#
# Copyright Val Curwen
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::SeqFetcher

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher->new(
                                             );
    $obj->pfetch('/path/to/pfetch');
    my $seq = $obj->run_pfetch('z87703');
    $obj->getz('/path/to/getz');
    my $seq2 = $obj->run_getz('z87703','embl emblnew');

=head1 DESCRIPTION

Object to perform various sequence retrieval functions

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher;

use strict;
use Bio::Root::RootI;
use Bio::Seq;
use Bio::SeqIO;

use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  
  return $self; # success - we hope!
}

=head2 pfetch

    Title   :   pfetch
    Usage   :   $self->pfetch('/usr/local/ensembl/bin/pfetch')
    Function:   Get/set for pfetch executable
    Returns :   
    Args    :   path to pfetch executable

=cut

sub pfetch {
  
  my ($self, $pfetch) = @_;
  if ($pfetch)
    {
      $self->{'_pfetch'} = $pfetch;
    }
  return $self->{'_pfetch'};  
  
}

=head2 efetch

    Title   :   efetch
    Usage   :   $self->efetch('/usr/local/ensembl/bin/efetch')
    Function:   Get/set for efetch executable
    Returns :   
    Args    :   path to efetch executable

=cut

sub efetch {
  
  my ($self, $efetch) = @_;
  if ($efetch)
    {
      $self->{'_efetch'} = $efetch;
    }
  return $self->{'_efetch'};  
  
}


=head2 getz

    Title   :   getz
    Usage   :   $self->getz('/usr/local/ensembl/bin/getz')
    Function:   Get/set for getz executable
    Returns :   
    Args    :   path to getz executable

=cut

sub getz {
  
  my ($self, $getz) = @_;
  if ($getz)
    {
      $self->{'_getz'} = $getz;
    }
  return $self->{'_getz'};  
  
}


=head2 run_pfetch

    Title   :   run_pfetch
    Usage   :   $self->run_pfetch($id)
    Function:   Retrieves a sequence using pfetch
    Returns :   Bio::Seq, or undef
    Args    :   id of sequence to be retrieved

=cut

sub run_pfetch {

  my ($self,$id) = @_;

  if (!defined($id)) {
    $self->throw("No id input to run_pfetch");
  }  

  my $seqstr;
  my $seq;
  my $newid      = $self->parse_header($id);
  my $pfetch     = $self->pfetch;

  # if pfetch path not explicitly set, assume it's in $PATH
  $pfetch        = 'pfetch' unless defined($pfetch); 

  open(IN,"$pfetch -q $newid |") || $self->throw("Error running pfetch for id [$newid]: $pfetch");
  $seqstr = <IN>;
  close IN;

  chomp($seqstr);
  if(defined $seqstr && $seqstr ne "no match") {
    $seq = new Bio::Seq(-seq => $seqstr,
			-id  => $newid);
  }
  
  return $seq;

}

=head2 run_efetch

    Title   :   run_efetch
    Usage   :   $self->run_efetch($id)
    Function:   Retrieves a sequence using efetch
    Returns :   Bio::Seq, or undef
    Args    :   id of sequence to be retrieved

=cut

sub run_efetch {

  my ($self,$id) = @_;

  if (!defined($id)) {
    $self->throw("No id input to run_efetch");
  }  

  my $seqstr;
  my $seq;
  my $newid      = $self->parse_header($id);
  my $efetch     = $self->efetch;

  # if efetch path not explicitly set, assume it's in $PATH
  $efetch        = 'efetch' unless defined($efetch); 

  open(IN,"$efetch -q $newid |") || $self->throw("Error running efetch for id [$newid]: $efetch");
  $seqstr = <IN>;
  close IN;

#  chomp($seqstr);
  $seq = new Bio::Seq(-seq => $seqstr, -id  => $newid)
    unless (!defined $seqstr || $seqstr =~ "not found");
  
  return $seq;

}


=head2 run_getz

    Title   :   run_getz
    Usage   :   $self->run_getz($id,$libs)
    Function:   Retrieves a sequence using getz from the specified libraries
    Returns :   Bio::Seq, or undef
    Args    :   id of sequence to be retrieved, string representing libraries to be searched

=cut

sub run_getz {

  my ($self,$id,$libs) = @_;

  if (!defined($id)) {
    $self->throw("No id input to run_getz");
  }

  if (!defined($libs)) {
    $self->throw("No libs input to run_getz");
  }

  my $seqstr;
  my $seq;
  my $newid      = $self->parse_header($id);
  $self->throw("Could not parse id [$id]") unless defined $newid;
  my $getz       = $self->getz;

  # if getz path not explicitly set, assume it's in $PATH
  $getz          = 'getz' unless defined($getz); 

  open(IN, "getz -e '[libs={$libs}-ID:$id] | [libs-AccNumber:$id]' |") 
    || $self->throw("Error running getz for id [$newid]: $getz");

  # hack just for rikens
  my $format = 'EMBL';
  if($libs eq 'mouseprot') { $format = 'Fasta'; }
  my $fh = Bio::SeqIO->new(-fh   => \*IN, "-format"=>$format);

  $seq = $fh->next_seq();
  close IN;

  $self->warn("Problem with getz for [$id]") unless defined $seq;

  return $seq;
  
}

=head2 parse_header

  Title   : parse_header
  Usage   : my $newid = $self->parse_header($id);
  Function: Parses different sequence headers
  Returns : string
  Args    :string to be parsed

=cut

sub parse_header {
    my ($self,$id) = @_;

    if (!defined($id)) {
	$self->throw("No id input to parse_header");
    }

    my $newid = $id;

    if ($id =~ /\/ug=(\S+)\s+/){
      $newid = $1;
    }
    
    elsif ($id =~ /^(.*)\|(.*)\|(.*)/) {
      if ($2 eq "UG") {
	$newid = $3;
      }
      else {
	$newid = $2;
      }
      $newid =~ s/(.*)\..*/$1/;
      
    }
 
    elsif ($id =~ /^..\:(.*)/) {
	$newid = $1;
    }

    $newid =~ s/ //g;

    print STDERR "newid: $newid\n";

    return $newid;
}

