
#
# Ensembl module for ScanProsite
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Profile - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Runnable::Protein::ScanProsite;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object



use Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation;
use Bio::Seq;
use Bio::SeqIO;


@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation);







sub multiprotein{
  my ($self) = @_;
  return 1;
}


=head2 run_analysis

    Title   :   run_analysis
    Usage   :   $obj->run_analysis
    Function:   Runs the blast query
    Returns :   nothing
    Args    :   none

=cut

sub run_analysis {
    my ($self) = @_;

    my $run = "/usr/local/bin/perl ".$self->program . ' -pattern ' .
      $self->database.' '.$self->parameters.' '.$self->filename . ' > ' 
        .$self->results;

    print STDERR "RUNNING: $run\n";
    
    $self->throw("Failed during ScanProsite run $!\n") unless 
      (system ($run) == 0) ;
}



=head2 parse_results

    Title   :  parse_results
    Usage   :   $obj->parse_results($filename)
    Function:   Parses cpg output to give a set of features
                parsefile can accept filenames, filehandles or pipes (\*STDIN)
    Returns :   none
    Args    :   optional filename

=cut
sub parse_results {
    my ($self) = @_;
    
    my $filehandle;
    my $resfile = $self->results();
    
    if (-e $resfile) {
        
      if (-z $self->results) {  
        print STDERR "pfscan didn't find any hits\n";
        return; }       
      else {
        open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
      }
    }
    my %printsac;
    my $line;
    
    while (<CPGOUT>) {
      $line = $_;
      chomp $line;
      #print STDERR "$line\n";
      my ($id,$hid,$name,$from,$to,$confirmed) = split (/\|/,$line);
      
      my $fp = $self->create_protein_feature($from, $to, $confirmed, $id,
                                             0, 0, $hid, $self->analysis,
                                             0, 0);
      $self->add_to_output($fp);
    }
    
     
}


1;
