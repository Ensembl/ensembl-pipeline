
#
# Ensembl module for Prints
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Prints - DESCRIPTION of Object

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


package Bio::EnsEMBL::Pipeline::Runnable::Protein::Prints;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object



use Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation;
use Bio::Seq;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation);




###########
# Analysis methods
##########



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

    
    # This routine expands the database name into $db-1 etc for
    # split databases

     my $command =  $self->program." ".$self->database." ".$self->filename.
       " "."-fjR  > ".$self->results, "\n";

    $self->throw("Failed during prints run $!\n") unless 
      (system ($command) == 0) ;
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
	    print STDERR "Printscan didn't find any hits\n";
	    return; }       
    else {
	    open (CPGOUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");#
    }
  }
  my %printsac;
  my $line;
  
  my $sequenceId;
  while (<CPGOUT>) {
    $line = $_;
    chomp $line;
    # Pattern match the Sn; field which should contain the SequenceId and Accession
    
    if ($line =~ s/^Sn;//) { # We have identified a Sn; line so there should be the following:
	    
	    #ENSP00000003603 Gene:ENSG00000000003 Query:AL035608 Contig:AL035608.00001 Chr:chrX basepair:97227305
	    ($sequenceId) = $line =~ /^\s*(\w+)/;
    }
    
    
    if ($line =~ s/^1TBH//) {
      my  ($id) = $line =~ /^\s*(\w+)/;
      my ($ac) = $line =~ /(PR\w+);?\s*$/;
      $printsac{$id} = $ac;
    }
    
    if ($line =~ s/^3TB//) {
      if ($line =~ s/^[HN]//) {
        my ($num,$temp1,$tot1) = "";
        # Grab these lines
        #       1433ZETA        1  of  6  88.19   1328    1.00e-16  ELTVEERNLLSVAYKNVIGARRASWRIITS                          30   35   36   48
        # split line on space, hence strip off all leading spaces first.
        $line =~ s/^\s+//;
        
        # Place all elements of list into an array      
        my @elements = split /\s+/, $line; 
        
        # Name each of the elements in the array
        my ($fingerprintName,$motifNumber,$temp,$tot,$percentageIdentity,$profileScore,$pvalue,$subsequence,$motifLength,$lowestMotifPosition,$matchPosition,$highestMotifPosition) = @elements;
        
        my $start = $matchPosition;
        #
        # If the match to the pattern lies at the end of the protein we might get padding of the subsequence with #'s, and the 
        # end position will be bigger than the actual end of the protein. So we'll strip the #'s off the end, adjust the 
        # motif length accordingly, and only then derive the match end.
        my $hash_substring;
        my $end;
        
        if($subsequence =~ /(\#+)$/){
          $hash_substring = $1;
          $end = $matchPosition + $motifLength - 1 - length($hash_substring);
        }else{
          $end = $matchPosition + $motifLength - 1;
        }
        
        my $print =  $printsac{$fingerprintName};
        
    
        my $fp = $self->create_protein_feature($start, $end, $profileScore,
                                               $sequenceId, 0, 0, $print, 
                                               $self->analysis, $pvalue, 
                                               $percentageIdentity);
        $self->add_to_output($fp);	
      }
    }
  }
}


1;

