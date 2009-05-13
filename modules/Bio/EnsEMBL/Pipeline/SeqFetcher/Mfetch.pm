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

  print " try to fetch $acc_to_fetch\n"  if $self->{verbose}; 
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

  print "cmd: $command\n" if $self->{verbose}; 
  open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc_to_fetch]:  $command ");  
   my @lines = <IN> ; 
  close IN or throw("Error running mfetch for accession [$acc_format]: $command ");  
 
  my %entry;  

  for my $line ( @lines ) { 
        chomp($line) ;    

        print "LINE $line\n" if $self->{verbose};
        # we're just parsing one entry so if we get a no_match we just return  .... 
        
        if ( $line =~m/no match/ ) {    
           print  "no entry found for $acc_to_fetch with $command\n" if $self->{verbose};
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
        print "populating $key_field ... ".join(" ", @l) . "\n" if $self->{verbose};
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


sub build_acc_index { 
   my ($acc) = @_ ;   

   my %tmp ;  
   my $position = 0 ; 
   for my $ac ( @$acc ) {   
     $tmp{$position} = $ac; 
     $position++;
   }  
   return \%tmp ; 
}

sub get_Entry_Fields_BatchFetch { 
  my ($self, $acc,$fields) = @_;

  unless ( ref($acc) =~m/ARRAY/ ) {   
    throw("if you use batchfetching you have to submit the acc-numbers in an array-reference\n"); 
  }  

  # NOTE : mfetch does currently not work for fetching with wildcards. use get_entryFields() 


  my $cmd_prefix = $self->_make_mfetch_prefix_command($fields);  

  my %acc_index = %{build_acc_index($acc)};
  # fetch in batches of 300   
  my @fetch_strings = @{make_fetch_strings($acc, 250 )};  

  my $command ;  
  my %all_entries; 
  my @entries_not_found;   
  my %all_entries; 
  my @entries_not_found;  
  my $last_acc; 
  my %last_entry ;  
  my $entry_number = 0;
  my @lines ;  

  # data fetching + parsing  
  
  my $last_field ="";     
  my ( %little_hash, %all_entries ) ; 
  my $new_entry = 0;  
  my $no_match_found = 0 ; 

  STRING: for my $acc_string ( @fetch_strings ) {   
    $command = $cmd_prefix ." " .  $acc_string;
    my @nr_acc = split /\s+/, $acc_string ;    
    print $entry_number . " / " . keys (%acc_index) . " entires fetched\n" ;  
    print "\n\n$command \n\n" if $self->{verbose} ; 
    ## data retrieval  
    

    open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc_string]: $command ");   
      my @lines  = <IN> ; 
    close IN or throw("Error running mfetch for accession [$acc_string]: $command");   
    
#    my ( %little_hash, %all_entries ) ;  

  #  my $new_entry = 0;  
  #  my $no_match_found = 0 ;  
  
    # data parsing  
   
    LINE: for my $line ( @lines ) {  
      chomp($line) ;   

      if ( $line =~m/no match/ ) {      
        $entry_number++;  
        print "no match for $entry_number -- $acc_index{$entry_number}\n" if $self->{verbose};; 
        $last_field = "";  
        $no_match_found++ ; 
        next LINE;  
      } 
      #     $hash{accession}{AC} 
      #     $hash{accession}{PE}    
      
      if ( $self->{verbose} ) { 
        print "Match $acc_index{$entry_number}  --> entry: " ;  
        print " $entry_number LINE : $line \n";  
      } 

      my @l = split /\s+/, $line;    
      my $field = shift @l ;     

      # parsing the start of the ENTRY with AC field  .... 
       
      if ($field =~m/AC/ ) {       
        if ( $last_field =~m/AC/)  {    
             $new_entry = 0 ;   
         } else { 
             print "\nnew entry found ...\n" if $self->{verbose} ; 
             $new_entry = 1; 
             $last_field = $field ;  
        } 
      }   

      if ( $new_entry == 1 ) {   
         if (scalar(keys %little_hash) > 0 ) { 
           if ( $field =~/AC/ ) {  
              print " NEW ENTRY STARTS - adding info to big hash ... $field ----- $last_field  -- $entry_number\n" if $self->{verbose}; 

              my $query_acc ; 
              if ( $no_match_found > 0 ) {  
                print "no match found ... $no_match_found \n" if $self->{verbose}; 
                $query_acc = $acc_index{($entry_number-$no_match_found) };  
                $no_match_found = 0; 
              } else { 
                $query_acc = $acc_index{$entry_number}; 
              }
              %all_entries = %{_add_information_to_big_hash(\%all_entries, \%little_hash,$query_acc )};   
              undef %little_hash;  
              $entry_number++; 
           } 
         }
      } 
        # add to little hash  
        $little_hash{$field}.=join (" " , @l); 
        $last_field = $field ;   
    } 
  } # fetch next round .. 
 
  # add last entry to all_entries . 
  
  if ( $self->{verbose}  ) { 
    for my $key ( keys %all_entries ) {  
       print "KEY $key\n" ;
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


sub _add_information_to_big_hash {
  my ($all, $little, $query_acc)  = @_ ;  

  my %all_entries = %$all; 
  my %little_hash = %$little; 

  # $little_hash{AC} =  Q7PYN8; Q01FQ6;
  # $little_hash{OC} =  Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota; 
  # $little_hash{PE} =  4: Predicted;
 
  my $acc_string = $little_hash{AC};
  $acc_string =~s/\;//g;
  my @accession_numbers = split /\s+/, $acc_string ;   

  # consistency check - the query acc which we used in the mfetch query should also be in the AC field of the entry returned ....  
  my $found ; 
  for my $returned_acc ( @accession_numbers ) {  
     if ($query_acc =~m/$returned_acc/) {  
        $found = 1; 
     } 
  } 
  unless ( $found ) {  
     throw( " the acc $query_acc can't be found in the query returned by mfetch [ " . join (" ", @accession_numbers) . " ]\n");  
  } 


    unless ( exists $all_entries{$query_acc} ) { 
      $all_entries{$query_acc} = \%little_hash; 
    }else { 
      # we already have an entry for this ACC. - check if the entries are the same ...  
      print "Warning ! The acc. you like to add has already been indexed ...\n" ; 
    
      # check if the entries are the same
      my %new_entry_to_add = %little_hash ; 
      my %stored_data = %{$all_entries{$query_acc}}; 

      for my $lk ( keys %little_hash ) {    
          if (  $little_hash{$lk}=~/$stored_data{$lk}/ ) { 
          } else {  
            print "DATA INCONSITENCY !!!\n" ; 
            print "NEW : $lk  --> $little_hash{$lk}\n" ; 
            print "OLD : $lk  --> $stored_data{$lk}\n" ;  
          } 
      } 
      print "\n\n\n" ; 
    }  
  return \%all_entries; 
}       


sub make_fetch_strings {  
   my ($acc_array, $size ) = @_;  

   my @acc_to_fetch = @$acc_array;  
   my @fetch_strings ; 

   if ( scalar(@acc_to_fetch) > $size ) {   
    while ( @acc_to_fetch ) {  
      my @tmp =  splice (@acc_to_fetch, 0, $size ) ;   
      push @fetch_strings, join ( " ", @tmp) ;  
      #print $fetch_strings[-1]."\n" ; 
    }  
   } else {  
     push @fetch_strings, join ( " ", @acc_to_fetch ) ; 
   }    
   return \@fetch_strings; 
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
