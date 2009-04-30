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
  $self->{verbose} = 0 ; 
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
  print STDERR "$command\n" if $self->{verbose}; 
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
  my ($self, $acc_to_fetch,$fields) = @_; 

   print "Fields to get : " . join ( " " , @$fields )."\n" if $self->{verbose}; 

  if ( ref($acc_to_fetch)=~m/ARRAY/ ) {   
    print "BatchFetch fields to get : " . join ( " " , @$fields )."\n" if $self->{verbose}; 
    return $self->get_Entry_Fields_BatchFetch($acc_to_fetch,$fields); 
  } 
 
  if (!defined($acc_to_fetch)) {
    throw("No accession input");
  }  

  print " try to fetch $acc_to_fetch\n" ; 

  my $command ;  
  my %all_entries; 
  my @entries_not_found;   

  my $cmd_prefix = $self->_make_mfetch_prefix_command($fields);  
  $command = $cmd_prefix ." " .  $acc_to_fetch; 

  # extract acc if you do a wildcard like mfetch -i  "AJFLD%"  

  my $acc_format = $acc_to_fetch;
  $acc_format=~s/acc://g;
  $acc_format=~s/\%//g;
  $acc_format=~s/\"//g; 

  print "cmd: $command\n"; 

  open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc_to_fetch]:  $command ");  
   my @lines = <IN> ; 
  close IN or throw("Error running mfetch for accession [$acc_format]: $command ");  
 
  my %entry;  

  for my $line ( @lines ) { 
        chomp($line) ;    

        print "LINE $line\n" ;  
        # we're just parsing one entry so if we get a no_match we just return  .... 
        
        if ( $line =~m/no match/ ) {    
           print  "no entry found for $acc_to_fetch with $command\n";
           push @entries_not_found, $acc_to_fetch; 
           return [\%all_entries , \@entries_not_found ] ; 
        }  

        my @l = split /\s+/, $line;     
        my $key_field = shift @l ;     

        # result can contain more than one line begining with same field identifier , ie  
        #   AC Q4589; 
        #   AC Q0999;
        #   not sure how this works if we only get AC's ....   but why would we do this anyway ?  
        
        if ( $key_field =~/AC/) {  
          if ( scalar(@l) > 1 ) {  
            warning ("more than one AC number returned : " . join(" ", @l) . "  - we ignore the others " . scalar(@l) . " objects\n") ; 
          }   
        }      

        $entry{$key_field}.= join(" ", @l);    
        #print "populating $key_field ... ".join(" ", @l) . "\n";
  }  

  # aim is to return $all_entries{Q4981}{AC}  -> string  
  # aim is to return $all_entries{Q4981}{org} -> string 
  # aim is to return $all_entries{ ACC }{ FIELD } -> string 
 
  # populate all_entry-hash 
  for my $field ( keys %entry ) {  
    $all_entries{$acc_format}{$field}=$entry{$field}; 
  }
        # my $hit=0; 
        # for my $result ( @l ) { 
            #print "testing : $result <-->  $tmp_acc .... \n" ;    
        #    if ( $result =~m/$tmp_acc/ ) {  
        #         $hit = 1 ;  
        #         last; 
        #    } 
        # } 


   for my $key ( keys %all_entries ) {  
       print "KEY $key\n";  
         my %tmp = %{$all_entries{$key}} ;
         for ( keys %tmp ) {   
            #if ( /AC/ ) { 
            #  print "\t\t$_ --> " . join(" ", @{$tmp{$_}} ). "\n";  
            #}else {  
              print "\t\t$_ --> $tmp{$_}\n";  
            #}
         }  
         print "\n" ;
   }    
  return [\%all_entries , \@entries_not_found ] ; 
} 


=head2 get_Entry_Fields_BatchFetch 

  Title   : batch_fetch
  Usage   : $self->batch_retrieval(@accession_list);
  Function: Retrieves multiple sequences via mfetch
  Returns : reference to a list of Bio::Seq objects
  Args    : array of accession strings

=cut

sub get_Entry_Fields_BatchFetch { 
  my ($self, $acc,$fields) = @_;
  print "using batchFetch\n" ; 
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
  my $last_acc; 
  my %last_entry ; 
  for my $acc_string ( @fetch_strings ) {  
    $command = $cmd_prefix ." " .  $acc_string;
    my @nr_acc = split /\s+/, $acc_string ;    
    my @tmp; 
    for my  $acc_format  ( @nr_acc) {    
         $acc_format=~s/acc://g;
         $acc_format=~s/\%//g;
         $acc_format=~s/\"//g;
         push @tmp, $acc_format ; 
     }
      
    print "cmd: $command\n" if $self->{verbose}; 

    open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc_string]: $mfetch");  
    @nr_acc = @tmp ; 
    my $acc;  
    my %entry_hash ;   
    my $accPointer = -1 ; 
    my $last_field ="";   

    LINE: while  (my $line=<IN>){  
        chomp($line) ;  
        #print "line $line\n";   
        #print "pointer : $nr_acc[$accPointer] $accPointer\n" ; 
        if ( $line =~m/no match/ ) {  
          #print  "no entry found for $nr_acc[$accPointer] \n" ;  
          $accPointer++;   
          #print "no -match - pointer val is : $accPointer\n" ;   
          #print "adding entry_not_found $nr_acc[$accPointer]\n"; 
          push @entries_not_found, $nr_acc[$accPointer] ; 
          next LINE;  
        } 
        my @l = split /\s+/, $line;    
        my $field = shift @l ;    
        #print "last_field : $last_field  VS $field\n" ;  
        # result can contain more than one line begining with same field identifier , ie  
        #   AC Q4589; 
        #   AC Q0999;
        #   not sure how this works if we only get AC's ....   but why would we do this anyway ? 
        #   
        
        if ($field =~m/AC/ ) {      
          if ( $field eq $last_field ) {  
             print "we got more than one line starting with AC ... wthis means it's not a new entry. \n" if $self->{verbose};
             # as we don't store/process the multipe acc returned we just ignore this line. this will caus 
             # trouble if we only have acc's returned.   
             $last_field = $field ; 
             next LINE; 
             
          } else {  
             print "NEW ENTRY : " . join ( " " , @l ) if $self->{verbose}; 
             $acc = $l[0];  
             $acc=~s/\;//; 

             $accPointer++ unless ( length(@nr_acc) == 1 ) ;  
   
             #if ( scalar(@l) > 1 ) {  
             #  warning ("more than one AC number returned : " . join(" ", @l) . "  - we ignore the others " . scalar(@l) . " objects\n") ; 
             #}  
             
             # this moves all already read information into the big hash before we 
             # read the new entry .  
             #
             
             if ( $entry_hash{AC} ) {  
               # entry is defined so we already read all AC information and we can add safely. 
               # problem if there's only one entry as it does not go here . ....  is the last entry always lost ? 
               my @acc_for_entry = @{$entry_hash{AC}} ; 
               for my $acc(@acc_for_entry) { 
                 for my $f ( keys %entry_hash ) { 
                    $all_entries{$acc}{$f}= $entry_hash{$f};  
                     #print "adding $acc $f $entry_hash{$f};\n"; 
                 }    
               }
             }  
             %last_entry = %entry_hash ; 
             undef %entry_hash;  
            
             # we now got a clean new entry_hash  for present entry  
             # $entry_hash{AC} = \@l;    
             #print "pointer val is : $accPointer\n" ;  
             my $tmp_acc = $nr_acc[$accPointer] ;  
             #print "test : $acc versus $tmp_acc \n" ;    

             $last_acc = $tmp_acc ;  

             my $hit=0; 
             for my $result ( @l ) { 
                #print "testing : $result <-->  $tmp_acc .... \n" ;    
                if ( $result =~m/$tmp_acc/ ) {  
                    $hit = 1 ;  
                    last; 
                } 
             } 

             unless ($tmp_acc =~m/$acc/ ) { 
                 if ( $hit == 0 ) { 
                   warning   " acc does not match - somethings wrong with the pointers to acc. $tmp_acc vs $acc\n" ;   
                 }
             }
             #print "entries match - $tmp_acc $acc " ;  
            # print "other_entries: " . join(" " , @l) . "\n"; 
               
             #$entry_hash{AC} = [$acc];     
            # print "populating $tmp_acc to entry_hash - NOT $acc\n"; 
             $entry_hash{AC} = [$tmp_acc];    
             $last_field = $field ; 
              next LINE ; 
           }  
        } 
       #print "populating $field ... ".join(" ", @l) . "\n";
        $entry_hash{$field}.= join (" ", @l);  
        $last_field = $field ; 
     }   
     %last_entry = %entry_hash ; 

     close IN or throw("Error running mfetch for accession [$acc_string]: $mfetch");  
  }  

  # add last entry to all_entries .
  unless ( $all_entries{$last_acc} ) {    
    # add content of last hash last_ENTRY HERE  - also populate last_entry as well beofre .  
    #print "need to add last_entry now \n" ;   

    if ( $last_entry{AC} ) { 
      my @acc_for_entry = @{$last_entry{AC}} ; 
       for my $acc(@acc_for_entry) { 
          for my $f ( keys %last_entry ) { 
             $all_entries{$acc}{$f}= $last_entry{$f};  
             #print "adding $acc $f $last_entry{$f};\n"; 
          }    
       }
     }  
   } else { 
     print "NO_ACC defined for last entry - skipping \n" if $self->{verbose} ; ;  
   } 


  if ( $self->{verbose}  ) { 
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
     print "no entry found for : $_\n" if $self->{verbose};  
  } 


  return [\%all_entries , \@entries_not_found ] ; 
} 








sub verbose { 
   my ($self, $arg ) = @_ ; 
   if ( $arg ) { 
      $self->{verbose}= $arg ;  
   }
}  


 
sub _make_mfetch_prefix_command {  
   my ($self, $f ) = @_ ; 
   my @fields = @$f ;  

   my $mfetch = $self->executable;  
   my $options = $self->options;
   if (defined($options)) { 
     unless ($options =~ /^-/ ) {
       $options = '-' . $options; 
     }
   }
  my $command = "$mfetch "; 
  if ( scalar(@fields) > 0 ) {   
    my $f = join (" ",@fields) ;  
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
  print "PREFIX COMMAND $command\n" ;  
  return $command ; 
}  

1;
