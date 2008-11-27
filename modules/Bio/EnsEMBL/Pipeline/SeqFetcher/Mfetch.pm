#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch->new(
							      -executable => $exe
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using mfetch.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

Method Bio::EnsEMBL::Root::_rearrange is deprecated.
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
rearrange(order, list); #instead

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch;

use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::DB::RandomAccessI;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $options) = rearrange(['EXECUTABLE', 'OPTIONS'], @args);
  
  if (!defined $exe) {
    $exe = 'mfetch';
  }
  $self->executable($exe);

  if (defined $options) {
    $self->options($options);
  } 
  # caching of sequences 
  $self->{_seq_cache}={};    

  return $self; # success - we hope!
}

=head2 executable

  Title   : executable
  Usage   : $self->executable('/path/to/executable');
  Function: Get/set for the path to the executable being used by the module. If not set, the executable is looked for in $PATH.
  Returns : string
  Args    : string

=cut

sub executable {
  my ($self, $exe) = @_;
  if ($exe)
    {
      $self->{'_exe'} = $exe;
    }
  return $self->{'_exe'};  
}

=head2 options

  Title   : options
  Usage   : $self->options('tc');
  Function: Get/set for options to mfetch
  Returns : string
  Args    : string

=cut

sub options {

  my ($self, $options) = @_;
  if ($options)
    {
      $self->{'_options'} = $options;
    }
  return $self->{'_options'};  

}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_eq_by_acc($accession);
  Function: Does the sequence retrieval via mfetch
  Returns : Bio::Seq
  Args    : 

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;
  
  if (!defined($acc)) {
    throw("No accession input");
  }  
 
  if ( defined ( $self->{_seq_cache}{$acc})) {  
     return $self->{_seq_cache}{$acc};
  } 
  # seqence not cached, needs to be fetched 
  my $seqstr;
  my $seq;
  my $mfetch = $self->executable; 
  # the option for fetching a sequence only by mfetch is mfetch -v fasta  
  my $options = $self->options;
  if (defined($options)) { $options = '-' . $options  unless $options =~ /^-/; }

  my $command = "$mfetch -v fasta  ";
  if (defined $options){
    $command .= "$options ";
  }
  
  $command .= $acc;
  print STDERR "$command\n"; 
  open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc]: $mfetch");  
   
  while  (my $line=<IN>){  
     chomp($line) ; 
     next if $line=~m/>/;  
     $seqstr.=$line; 
  }  

  close IN or throw("Error running mfetch for accession [$acc]: $mfetch");
  chomp($seqstr);
  eval{
    if(defined $seqstr && $seqstr ne "no match") {
      $seq = new Bio::Seq('-seq'               => $seqstr,
			  '-accession_number'  => $acc,
			  '-display_id'        => $acc);
    }
  };

  if($@){
    print STDERR "$@\n";
  }
  
  throw("Could not mfetch sequence for [$acc]\n") unless defined $seq;

  $self->{_seq_cache}{$acc}=$seq;
  return $seq;
} 


=head2 get_Entry_fields

  Title   : get_Entry_Fields
  Usage   : $self->get_Entry_Fields(\@accessions,\@fields); 
          : $self->get_Entry_Fields("Q5RFX5", \@fields); 
  Arg [0] : $self 
  arg [1] : ACC as string ( Q5RFX5 ) OR arrayref to array of acc 
  arg [2] : arrayref to fields which you want to receive     
            \@field = qw(  pe taxon acc )  ;  
  Function: Does the retrieval of different files like Taxonomy_id or PE level via mfetch, 
            either for one acc or in batch mode.
  Returns : arrayref to array of hashes for each acc.?
  Args    : 

=cut

sub  get_Entry_Fields {
  my ($self, $acc,$fields) = @_;
 
  if ( ref($fields)=~m/ARRAY/ ) {  
     print "fields to get : " . join ( " " , @$fields )."\n"; 
  }  
  if (!defined($acc)) {
    throw("No accession input");
  }  

  my @acc_to_fetch ; 
  if ( ref($acc) =~m/ARRAY/ ) {   
    # batch mode - more than 1 acc to fetch 
     @acc_to_fetch = @$acc; 
  }else {  
    #have only single acc to fetch 
    push @acc_to_fetch,$acc; 
  }  

  my $seqstr;
  my $seq;
  my $mfetch = $self->executable;  

  my $options = $self->options;
  if (defined($options)) { $options = '-' . $options  unless $options =~ /^-/; }

  my $command = "$mfetch "; 
  if ( scalar(@$fields) > 0 ) {   
    my $f = join (" ",@$fields) ;  
    if ( $f=~m/acc/) {  
      $f =~s/acc//g; 
    }
    $f = "acc " .$f; 
    $command .= " -f \"$f \" "      ;
  }else {  
    $command .= " -v  full"  ;   
  }   

  if (defined $options){
    $command .= "$options ";
  }

  # we will hand over the acc to fetch as string.  
  # if in batch mode do not submit more than 50 entries / strings to mfetch at one time.  
  my @fetch_strings ; 
  if ( scalar(@acc_to_fetch) > 500 ) {   
    # more than 500 acc to fetch - make strings of 500 acc which is faster.
    while ( @acc_to_fetch ) {  
      my @tmp =  splice (@acc_to_fetch, 0, 300 ) ;   
      push @fetch_strings, join ( " ", @tmp) ;  
      #print $fetch_strings[-1]."\n" ; 
    }  
  } else {  
    push @fetch_strings, join ( " ", @acc_to_fetch ) ; 
  }   
  
  my $cmd_prefix = $command ;   

  my %all_entries; 
  my @entries_not_found;  

  for my $acc_string ( @fetch_strings ) {  
    $command = $cmd_prefix ." " .  $acc_string;
    my @nr_acc = split /\s+/, $acc_string ; 
    print "cmd: $command\n";  

    open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc_string]: $mfetch");  

    my $acc;  
    my %entry_hash ;   
    my $accPointer = -1 ; 
    my $last_field ="";   

    LINE: while  (my $line=<IN>){  
        chomp($line) ;    
        #print "pointer : $nr_acc[$accPointer] $accPointer\n" ; 
        if ( $line =~m/no match/ ) {  
          #print  "no entry found for $nr_acc[$accPointer] - exiting\n" ;  
          $accPointer++;   
          print "no -match - pointer val is : $accPointer\n" ;  
          push @entries_not_found, $nr_acc[$accPointer] ; 
          next LINE;  
        } 
        my @l = split /\s+/, $line;    
        my $field = shift @l ;    
        print "last_field : $last_field  VS $field\n" ;  
        # result can contain more than one line begining with same field identifier , ie  
        #   AC Q4589; 
        #   AC Q0999;
        #   not sure how this works if we only get AC's ....   but why would we do this anyway ? 
        #   
        
        if ($field =~m/AC/ ) {      
          if ( $field eq $last_field ) {  
             print "we got more than one line starting with AC ... wthis means it's not a new entry. \n" ; 
             # as we don't store/process the multipe acc returned we just ignore this line. this will caus 
             # trouble if we only have acc's returned.   
             $last_field = $field ; 
             next LINE; 
             
          } else {  
             print "NEW ENTRY : " . join ( " " , @l ) ;  
             $acc = $l[0];  
             $acc=~s/\;//; 

             $accPointer++;  
             print "increment pointer : $accPointer\n" ;  
             print "test : $acc vs $nr_acc[$accPointer]\n" ;  
   
             #if ( scalar(@l) > 1 ) {  
             #  warning ("more than one AC number returned : " . join(" ", @l) . "  - we ignore the others " . scalar(@l) . " objects\n") ; 
             #}  
             
             # this moves all already read information into the big hash before we 
             # read the new entry . 
             if ( $entry_hash{AC} ) {  
               # entry is defined so we already read all AC .... 
               my @acc_for_entry = @{$entry_hash{AC}} ; 
               for my $acc(@acc_for_entry) { 
                 for my $f ( keys %entry_hash ) { 
                    $all_entries{$acc}{$f}= $entry_hash{$f}; 
                 }    
               }
             } 
             undef %entry_hash;  
            
             # we now got a clean new entry_hash  for present entry  
             # $entry_hash{AC} = \@l;    
             print "pointer val is : $accPointer\n" ;  
             my $tmp_acc = $nr_acc[$accPointer] ;  
             unless ($tmp_acc =~m/$acc/ ) {  
               throw  " acc does not match - somethings wrong with the pointers to acc. $tmp_acc vs $acc\n" ;  
             } 
             #$entry_hash{AC} = [$acc];    
             $entry_hash{AC} = [$tmp_acc];    
             $last_field = $field ; 
              next LINE ; 
           }  
        }
        $entry_hash{$field}.= join (" ", @l);  
        $last_field = $field ; 
     }   
     close IN or throw("Error running mfetch for accession [$acc_string]: $mfetch");  
  } 
  if ( 0 ) { 
    for my $key ( keys %all_entries ) {  
       print "KEY $key\n";  
         my %tmp = %{$all_entries{$key}} ;
         for ( keys %tmp ) {   
            if ( /AC/ ) { 
              print "\t\t$_ --> " . join(" ", @{$tmp{$_}} ). "\n";  
            }else {  
              print "\t\t$_ --> $tmp{$_}\n";  
            }
         }  
         print "\n" ;
    }   
  }
  for ( @entries_not_found ) {  
     print "no entry found for : $_\n" ; 
  } 
  return [\%all_entries , \@entries_not_found ] ; 
} 

##
##
##=head2 batch_fetch
##
##  Title   : batch_fetch
##  Usage   : $self->batch_retrieval(@accession_list);
##  Function: Retrieves multiple sequences via mfetch
##  Returns : reference to a list of Bio::Seq objects
##  Args    : array of accession strings
##
##=cut
##
##sub  batch_fetch {
##  my $self = shift @_;
##
##  my @sequence_list;
##  
##  unless (scalar @_) {
##    throw("No accession input");
##  }  
##
##  my $accession_concatenation;
##  my @accession_list;
##
##  while (my $acc = shift @_) {
##    push (@accession_list, $acc);
##    $accession_concatenation .= $acc . ' ';
##  }
##  
##  my $mfetch = $self->executable;
##  my $options = $self->options;
##  if (defined($options)) { $options = '-' . $options  unless $options =~ /^-/; }
##
##  my $command = "$mfetch -q ";
##  if (defined $options){
##    $command .= "$options ";
##  }
##
##  $command .= $accession_concatenation;
##
##  #print STDERR "$command\n";
##
##  open(IN,"$command |") or throw("Error opening pipe to mfetch : $mfetch");
##
##  my $seq_placemarker = -1;
##
## SEQ:
##  while (my $seqstr = <IN>) {
##    
##    $seq_placemarker++;
##
##    chomp($seqstr);
##    
##    my $seq;
##    
##    unless (defined $seqstr && $seqstr eq 'no match') {
##
##      eval{
##	if(defined $seqstr && $seqstr ne "no match") {
##	  $seq = new Bio::Seq('-seq'               => $seqstr,
##			      '-accession_number'  => $accession_list[$seq_placemarker],
##			      '-display_id'        => $accession_list[$seq_placemarker]);
##	}
##      };
##      
##      if($@){
##	print STDERR "$@\n";
##      }
##    }
##
##    unless (defined $seq){
##      $self->warn("PFetch Error : Could not mfetch sequence for " . 
##		  $accession_list[$seq_placemarker] . "\n");
##      next SEQ;
##    }
##
##    push (@sequence_list, $seq);
##  }
##  
##  close IN or throw("Error running mfetch : $mfetch");
##  
##  return \@sequence_list;
##}
##
##
##
1;
