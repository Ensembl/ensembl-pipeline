=head1 LICENSE

 Copyright [1999-2013] Genome Research Ltd. and the EMBL-European Bioinformatics Institute
 
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 
      http://www.apache.org/licenses/LICENSE-2.0
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch - 

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Mfetch->new(
							      -executable => $exe
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using mfetch.

=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

Method Bio::EnsEMBL::Root::_rearrange is deprecated.
use warnings ;
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
use Time::HiRes qw ( gettimeofday ) ; 
use vars qw(@ISA);
$|=0;
@ISA = qw(Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $options ) = rearrange(['EXECUTABLE', 'OPTIONS'], @args);
 
  if (!defined $exe) {
    $exe = 'mfetch';
  }
  $self->executable($exe); 

  if (defined $options) { 
    $options=~s/^\s+//g; 
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
  if (defined $exe)
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
  if ($options) {
      if ( $self->{_options} ) {
         $self->{'_options'} = $self->{_options} . " " . $options;
      } else {
        $self->{'_options'} = $options;
     }
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

sub get_Seq_by_acc {
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
  close IN or throw("Error running mfetch for accession [$acc]: \n$mfetch");
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


=head2 get_Entry_Fields_no_mfetch

  Title   : get_Entry_Fields_no_mfetch
  Usage   : $self->get_Entry_Fields_no_mfetch(\@accessions,\@fields,"/nfs/ensembl/amonida/uniprot_sprot_summary.dat","/nfs/ensembl/amonida/uniprot_trembl_summary.dat"); 
          : $self->get_Entry_Fields_no_mfetch("Q5RFX5", \@fields,"/nfs/ensembl/amonida/uniprot_sprot_summary.dat","/nfs/ensembl/amonida/uniprot_trembl_summary.dat"); 
  Arg [0] : $self
  arg [1] : ACC as string ( Q5RFX5 ) OR arrayref to array of acc 
  arg [2] : arrayref to fields which you want to receive
            \@field = qw(  pe taxon acc )  ;
  arg [3] : Filepath to uniprot dumped file containing Swiss-Prot protein fields.
  arg [4] : Filepath to uniprot dumped file containing TrEMBL protein fields.
  Function: Does the retrieval of different files like Taxonomy_id or PE level via uniprot dumped files.
  Returns : arrayref to array of hashes for each acc.?
  Args    : 

=cut

sub get_Entry_Fields_no_mfetch {
  my ($self,$sprot_file,$trembl_file,$acc_to_fetch,$fields) = @_; 

   print "Fields to get : " . join ( " " , @$fields )."\n" if $self->{verbose}; 

  if ( ref($acc_to_fetch)=~m/ARRAY/ ) {
    print "BatchFetch fields to get : " . join ( " " , @$fields )."\n" if $self->{verbose}; 

    my @acc_to_fetch_no_version = ();
    foreach my $an_acc (@$acc_to_fetch)
    {
      my $an_acc_no_version = $an_acc;
      $an_acc_no_version =~ s/\..*//;
      push(@acc_to_fetch_no_version,$an_acc_no_version);
    }
    return $self->get_Entry_Fields_BatchFetch_no_mfetch($sprot_file,$trembl_file,\@acc_to_fetch_no_version,$fields);
  }
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

sub get_Entry_Fields {
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
  my @lines = @{$self->_mfetch_command($command)} ; 
  #open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc_to_fetch]:  $command ");  
  # my @lines = <IN> ; 
  #close IN or throw("Error running mfetch for accession [$acc_format]: $command ");  
 
  my %entry;  

  for my $line ( @lines ) { 
        chomp($line) ;    
        #print "LINE $line" if $self->{verbose};
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
  # populate all_entry-hash 
  for my $field ( keys %entry ) {  
    $all_entries{$acc_format}{$field}=$entry{$field}; 
  }

   if ( 0 )  { 
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
  }
  return [\%all_entries , \@entries_not_found ] ; 
} 


sub _mfetch_command {
  my  ($self, $command)= @_;

  open(IN,"$command |") or throw("Error opening pipe to mfetch for command  $command ");
  my @lines = <IN> ;
  close IN or throw("Error running mfetch for command $command ");
  return \@lines ;
}


=head2 get_Seq_by_acc_wildcard

  Title   : get_Seq_by_acc_wildcard
  Usage   : $self->get_eq_by_acc($accession);
  Function: Does the sequence retrieval via mfetch
  Returns : string  
  Args    : 

=cut

sub get_Seq_by_acc_wildcard {
  my ($self, $acc) = @_;

  $acc=~s/\..*//g; # chop version off 

  my $options = $self->options ;
  unless ( $options=~m/-v fasta/ ) {
    $self->options("-v fasta");
  }
  my $cmd_prefix = $self->_make_mfetch_prefix_command(\[]) ;
  $cmd_prefix .= " -i \"acc:$acc\%\""  ;

  my @entry;
  my $not_found;

  my @lines = @{ $self->_mfetch_command($cmd_prefix)};
  chomp(@lines) ;
  for ( @lines ) {
    if (/no match/) {
      return \@entry;
    }
  }
 return \@lines ;
}

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

sub get_Entry_Fields_file {
  my ($filename,@acc) = @_;

  my @lines = ();
  my $acc_found = 0;

  my %acc_hash = map { $_ => 1 } @acc;

  print("Starting get_Entry_Fields_file $filename, num acc=".scalar(@acc)."\n");

  open(FILE,$filename) || die "Could not open file $filename";

  while (my $line=<FILE>) {
    if ($acc_found == 1) {
      if ($line =~ /^PE/) { # end of fields for the accession found, last line to store for the current accession
        # adding the DT field artificially from the filename
        if ($filename =~ /uniprot_sprot/) {
          push(@lines,"DT   00-AAA-0000, integrated into UniProtKB/Swiss-Prot.");
        } else {
          push(@lines,"DT   00-AAA-0000, integrated into UniProtKB/TrEMBL.");
        }
        $acc_found = 0;
      }
      push(@lines,$line);
    } else {
      if ($line =~ /^AC/) {
        my $first_AC_line = $line;
        my @acc_line_fields = split /\s+/,$line;
        shift(@acc_line_fields); # remove AC field

        my $line_after_AC_line;

        while ($line=<FILE>) {
          if ($line !~ /^AC/) {
            $line_after_AC_line = $line;
            last;
          } else { # add more accessions from additional AC lines to our AC line
            my @more_acc_line_fields = split /\s+/,$line;
            shift(@more_acc_line_fields); # remove AC field
            push(@acc_line_fields,@more_acc_line_fields);
          }
        }

	# we only want to print the accs that we can find in our accs hash, not all of the line where at least we have 1 match
        my @acc_found = ();
        foreach my $an_acc (@acc_line_fields) {
	  $an_acc =~ s/;//; # remove last character ;
          if (exists $acc_hash{$an_acc}) {
            push(@acc_found,$an_acc);
            $acc_found = 1;
            #print "found an_acc: $an_acc\n";
          }
          #last if $acc_found == 1;
        }
        if ($acc_found == 1) {
          my $acc_str = join(" ",@acc_found);
          push(@lines,"AC $acc_str\n");
          push(@lines,$line_after_AC_line);
        }
      }
    }
  }
  close(FILE);
  return @lines;
}

sub get_Entry_Fields_BatchFetch_no_mfetch {
  my ($self,$sprot_file,$trembl_file,$acc,$fields) = @_;

  unless ( ref($acc) =~m/ARRAY/ ) {   
    throw("if you use batchfetching you have to submit the acc-numbers in an array-reference\n"); 
  }

  my $cmd_prefix = $self->_make_mfetch_prefix_command($fields);   
  my %acc_index = %{build_acc_index($acc)};

  my @acc_to_fetch = @$acc;

  my @fetch_strings;
  push @fetch_strings,join( " ",@acc_to_fetch); # everything is fetched

  print "got " . scalar(@fetch_strings) . " jobs to run \n"; 
  my $command ;  
  my @entries_not_found;   
  my $last_acc; 
  my %last_entry ;  
  my $entry_number = 0;
  my @lines ;  
  # data fetching + parsing  
  my %no_match_index;  
  my $last_field ="";     
  #my ( %little_hash, %all_entries ) ;  #c 
  my ( $little_hash, $all_entries ) ; 
  my $new_entry = 0;  
  #my $no_match_found = 0 ; 

  STRING: for my $acc_string ( @fetch_strings ) {    
    $command = $cmd_prefix ." " .  $acc_string;
    my @nr_acc = split /\s+/, $acc_string ;      

    print $entry_number . " / " . keys (%acc_index) . " entries fetched\n" ;  

    my $t0 = gettimeofday() ;  
    print "starting no_mfetch\n" ;  

    my @lines_sprot = get_Entry_Fields_file($sprot_file,@nr_acc);
    print ("sprot finished, num of lines_sprot:".scalar(@lines_sprot)."\n");
    my @lines_trembl = get_Entry_Fields_file($trembl_file,@nr_acc);
    print ("trembl finished, num of lines_trembl:".scalar(@lines_trembl)."\n");

    my @unsorted_lines = ();
    push(@unsorted_lines,@lines_sprot);
    push(@unsorted_lines,@lines_trembl);

    print("SORTING NOW!\n");

    my $is_found = 0;
    # have to sort lines so that the accessions are in the same order as in acc_string
    ONE_ACC: foreach my $one_acc (@nr_acc) {
           foreach my $one_line (@unsorted_lines) {
             if ($is_found == 0) {
               if ( ($one_line =~ /^AC/) and ($one_line =~ /$one_acc/) ) {
                 push(@lines,"AC $one_acc");
                 $is_found = 1;
               } else {
                 ;
               }
             } else {
               push(@lines,$one_line);
               if ($one_line =~ /^PE/) { # end of fields for the accession found, last line to store for the current accession
                 $is_found = 0;
                 next ONE_ACC;
               }
             }
           } # end foreach one_line
           push(@lines,"no match\n");
         } # end foreach one_acc

    my $t1 = gettimeofday() ; 
    my $delta_t = $t1 - $t0 ;  
    print "time for no_mfetch : $delta_t\n" ; 

    # data parsing  

    LINE: for my $line ( @lines ) {
      print "PARSING : $line\n" if $self->verbose(); 
      chomp($line) ;

      if ( $line =~m/no match/ ) {

        print "line contains \"no match\"  \n" if $self->{verbose} ;  
        $last_field = "";

        # if we have an entry in memory store it
         if (scalar(keys %$little_hash) > 0 ) {   # if we have read a full entry before which has not been stored 
              print " have no_match now, but a full entry in memory for ".  $acc_index{$entry_number} . "-- so let's try and store it.\n" if $self->{verbose}; 
              my  $query_acc = $acc_index{$entry_number}; 
              $all_entries = _add_information_to_big_hash($all_entries, $little_hash,$query_acc );   
              undef $little_hash;  
              $entry_number++ ;
              print "stored and entry incremented : $entry_number   $acc_index{$entry_number} NO_MATCH \n" ; 
              print "NEW adding $acc_index{$entry_number} $entry_number to the list of no_matches \n" if $self->{verbose};
              push @entries_not_found,  $acc_index{$entry_number};  
              $no_match_index{$entry_number}=1; 
              next LINE;  
         } 
          print "no match for $entry_number -- $acc_index{$entry_number}\n" if $self->{verbose};;   
          if ( exists $no_match_index{$entry_number} ) {  
            $entry_number++ ; 
          } 
          $no_match_index{$entry_number}=1;  
          print "adding $acc_index{$entry_number} $entry_number to the list of no_matches \n" if $self->{verbose}; 
          push @entries_not_found,  $acc_index{$entry_number};  
          $entry_number++ ; 
          next LINE;  
      }  

      my @l = split /\s+/, $line;    
      my $field = shift @l ;     

      # parsing the start of the ENTRY with AC field  .... 

      if ($field =~m/AC/ ) {
        if ( $last_field =~m/AC/)  {
             # we have multiple AC fields in the entry which follow each other  AC xxx AC xxx 
             $new_entry = 0 ;
         } else {
             print "\nnew entry found ...\n" if $self->{verbose} ;
             $new_entry = 1;
             $last_field = $field ;
        }
      }

      if ( $new_entry == 1 ) {    
         if (scalar(keys %$little_hash) > 0 ) {   # if we have read a full entry before which has not been stored 
           if ( $field =~/AC/ ) {                # if we NOW READ a  new entry we need to store the last one ...  

              print " NEW ENTRY STARTS\n" if $self->{verbose} ; 
              # determine correct entry index 
              while (  exists $no_match_index{$entry_number} ) { 
                print $acc_index{$entry_number} . " is recorded as NO MATCH - checking next ...\n" if $self->{verbose} ; 
                $entry_number++ ; 
              } 
              my $query_acc ; 
              #if ( $no_match_found > 0 ) {  
              #  print "no matches found ... $no_match_found \n" if $self->{verbose}; 
              #  $query_acc = $acc_index{($entry_number-$no_match_found) };  
              #  $no_match_found = 0; 
              #} else { 
                $query_acc = $acc_index{$entry_number}; 
              #}
              $all_entries = _add_information_to_big_hash($all_entries, $little_hash,$query_acc );   


              undef $little_hash;
              $entry_number++;
           }
         }
        elsif ( exists $no_match_index{$entry_number} ) {  
            warning("entry with number $entry_number  ( $acc_index{$entry_number} ) was recorded as no_match  -incrementing entry ... \n");
            $entry_number++; 
        }
      }
      # add to little hash  
      $$little_hash{$field}.=join (" " , @l); 
      $last_field = $field ;   
    }  # next LINE 
  } # next STRING - fetch next round 
 
  # add last entry to all_entries . 
  $all_entries = _add_information_to_big_hash($all_entries, $little_hash,$acc_index{$entry_number} );   
  
  if ( $self->{verbose}  ) { 
    for my $key ( keys %{$all_entries} ) {  
       print "KEY $key\n" ;
         my %tmp = %{$$all_entries{$key}} ;
         for ( keys %tmp ) {   
           print "\t\t$_ --> $tmp{$_}\n";  
         }
         print "\n" ;
    }
  }
  # combine both  
  return [$all_entries , \@entries_not_found ] ; 
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

  unless ( ref($acc) =~m/ARRAY/ ) {   
    throw("if you use batchfetching you have to submit the acc-numbers in an array-reference\n"); 
  }  

  # NOTE : mfetch does currently not work for fetching with wildcards. use get_entryFields() 


  my $cmd_prefix = $self->_make_mfetch_prefix_command($fields);   
  my %acc_index = %{build_acc_index($acc)};
  # fetch in batches of 300   
  my @fetch_strings = @{make_fetch_strings($acc, 1000 )};   
  print "got " . scalar(@fetch_strings) . " jobs to run \n"; 
  my $command ;  
  my @entries_not_found;   
  my $last_acc; 
  my %last_entry ;  
  my $entry_number = 0;
  my @lines ;  
  # data fetching + parsing  
  my %no_match_index;  
  my $last_field ="";     
  #my ( %little_hash, %all_entries ) ;  #c 
  my ( $little_hash, $all_entries ) ; 
  my $new_entry = 0;  
  #my $no_match_found = 0 ; 

  STRING: for my $acc_string ( @fetch_strings ) {    
    $command = $cmd_prefix ." " .  $acc_string;
    my @nr_acc = split /\s+/, $acc_string ;      

    print $entry_number . " / " . keys (%acc_index) . " entries fetched\n" ;  

    if ( $self->{verbose} ) { 
      print "\n\n$command \n\n" ;
    }
    
    my $t0 = gettimeofday() ;  
    print "starting mfetch\n" ;  
    my @lines = @{$self->_mfetch_command($command)} ; 
    my $t1 = gettimeofday() ; 
    my $delta_t = $t1 - $t0 ;  
    print "time for mfetch : $delta_t\n" ; 
    
    # data parsing  
   
    LINE: for my $line ( @lines ) {  
      print "PARSING : $line\n" if $self->verbose(); 
      chomp($line) ;   

      if ( $line =~m/no match/ ) {        
# print "line $line\n"; 
        print "line contains \"no match\"  \n" if $self->{verbose} ;  
        $last_field = "";    

        # if we have an entry in memory store it  
        
         if (scalar(keys %$little_hash) > 0 ) {   # if we have read a full entry before which has not been stored 
              print " have no_match now, but a full entry in memory for ".  $acc_index{$entry_number} . "-- so let's try and store it.\n" if $self->{verbose}; 
              my  $query_acc = $acc_index{$entry_number}; 
              $all_entries = _add_information_to_big_hash($all_entries, $little_hash,$query_acc );   
              undef $little_hash;  
              $entry_number++ ;
              print "stored and entry incremented : $entry_number   $acc_index{$entry_number} NO_MATCH \n" ; 
              print "NEW adding $acc_index{$entry_number} $entry_number to the list of no_matches \n" if $self->{verbose};
              push @entries_not_found,  $acc_index{$entry_number};  
              $no_match_index{$entry_number}=1; 
              next LINE;  
         } 
          print "no match for $entry_number -- $acc_index{$entry_number}\n" if $self->{verbose};;   
          if ( exists $no_match_index{$entry_number} ) {  
            $entry_number++ ; 
          } 
          $no_match_index{$entry_number}=1;  
          print "adding $acc_index{$entry_number} $entry_number to the list of no_matches \n" if $self->{verbose}; 
          push @entries_not_found,  $acc_index{$entry_number};  
          $entry_number++ ; 
          next LINE;  
      }  

      my @l = split /\s+/, $line;    
      my $field = shift @l ;     

      # parsing the start of the ENTRY with AC field  .... 
       
      if ($field =~m/AC/ ) {       
        if ( $last_field =~m/AC/)  {     
             # we have multiple AC fields in the entry which follow each other  AC xxx AC xxx 
             $new_entry = 0 ;   
         } else { 
             print "\nnew entry found ...\n" if $self->{verbose} ; 
             $new_entry = 1; 
             $last_field = $field ;  
        } 
      }   

      if ( $new_entry == 1 ) {    
         if (scalar(keys %$little_hash) > 0 ) {   # if we have read a full entry before which has not been stored 
           if ( $field =~/AC/ ) {                # if we NOW READ a  new entry we need to store the last one ...  
              
              print " NEW ENTRY STARTS\n" if $self->{verbose} ; 
              # determine correct entry index 
              while (  exists $no_match_index{$entry_number} ) { 
                print $acc_index{$entry_number} . " is recorded as NO MATCH - checking next ...\n" if $self->{verbose} ; 
                $entry_number++ ; 
              } 
              my $query_acc ; 
              #if ( $no_match_found > 0 ) {  
              #  print "no matches found ... $no_match_found \n" if $self->{verbose}; 
              #  $query_acc = $acc_index{($entry_number-$no_match_found) };  
              #  $no_match_found = 0; 
              #} else { 
                $query_acc = $acc_index{$entry_number}; 
              #}
              $all_entries = _add_information_to_big_hash($all_entries, $little_hash,$query_acc );   


              undef $little_hash;  
              $entry_number++; 
           }  
         } elsif ( exists $no_match_index{$entry_number} ) {  
             warning("entry with number $entry_number  ( $acc_index{$entry_number} ) was recorded as no_match  -incrementing entry ... \n");
             $entry_number++; 
         } 
      }
      # add to little hash  
      $$little_hash{$field}.=join (" " , @l); 
      $last_field = $field ;   
    }  # next LINE 
  } # next STRING - fetch next round 
 
  # add last entry to all_entries . 
  $all_entries = _add_information_to_big_hash($all_entries, $little_hash,$acc_index{$entry_number} );   
  
  if ( $self->{verbose}  ) { 
    for my $key ( keys %{$all_entries} ) {  
       print "KEY $key\n" ;
         my %tmp = %{$$all_entries{$key}} ;
         for ( keys %tmp ) {   
           # if ( /AC/ ) { 
           #   print "\t\t$_ --> " . join(" ", @{$tmp{$_}} ). "\n"; 
           # }else {  
              print "\t\t$_ --> $tmp{$_}\n";  
           # }
         }  
         print "\n" ;
    }   
  } 

  # combine both  

  return [$all_entries , \@entries_not_found ] ; 
} 


sub get_Seq_BatchFetch {
  my ($self, $acc ) = @_;

  unless ( ref($acc) =~m/ARRAY/ ) {
    throw("if you use batchfetching you have to submit the acc-numbers in an array-reference\n");
  }

  # NOTE : mfetch does currently not work for fetching with wildcards. use get_entryFields() 

  my $options = $self->options ;
  unless ( $options=~m/-v fasta/ ) {
    $self->options("-v fasta") ;
  }
  my $cmd_prefix = $self->_make_mfetch_prefix_command([]) ;

  my %acc_index = %{build_acc_index($acc)};

  # fetch in batches of 500   
  my @fetch_strings = @{make_fetch_strings($acc, 500 )};


  my (@clean , @entries_not_found, @lines ) ;
  my ($command) ;
  my $entry_number = 0;

  # data fetching + parsing  
  STRING: for my $acc_string ( @fetch_strings ) {
    $command = $cmd_prefix ." " .  $acc_string;
    my @nr_acc = split /\s+/, $acc_string ;

    # data retrieval   
    #print $command ;  
    my @lines = @{$self->_mfetch_command($command)} ; 
    #open(IN,"$command |") or throw("Error opening pipe to mfetch for accession [$acc_string]: $command ");
    #my @lines  = <IN> ;
    #close IN or throw("Error running mfetch for accession [$acc_string]: $command");


    LINE: for my $line ( @lines ) {
      chomp($line) ;

      if ( $line =~m/no match/ ) {
        $entry_number++;
        print "no match for $acc_index{$entry_number}\n" ;
        push @entries_not_found , $acc_index{$entry_number};
        next LINE;
      } else {
        push @clean, $line ;
      }
    }
  } # fetch next round .. 

  for ( @entries_not_found ) {
    print "no entry found for : $_\n" if $self->{verbose};
  }
  return  [\@clean, \@entries_not_found ] ;
}











sub _add_information_to_big_hash {
  my ($all, $little, $query_acc)  = @_ ;  

  #my %all_entries = %$all;  #c 
  #my %little_hash = %$little; #c  
  #my %little_hash = %$little; #c  

  # $little_hash{AC} =  Q7PYN8; Q01FQ6;
  # $little_hash{OC} =  Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota; 
  # $little_hash{PE} =  4: Predicted;
 
  my $acc_string = $$little{AC};
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


    unless ( exists $$all{$query_acc} ) { 
      $$all{$query_acc} = $little;
    }else { 
      # we already have an entry for this ACC. - check if the entries are the same ...  
      print "Warning ! The acc. you like to add has already been indexed ...\n" ; 
    
      # check if the entries are the same

      for my $lk ( keys %$little) {    
          if (  $$little{$lk}=~/$$all{$query_acc}{$lk}/ ) { 
          } else {  
            warning( "POSSIBLE DATA INCONSITENCY !!!\n" );
            print "NEW : $lk  --> $$little{$lk}\n" ; 
            print "OLD : $lk  --> $$all{$query_acc}{$lk}\n";
          } 
      } 
      print "\n\n\n" ; 
    }  
  return $all; 
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
   if (defined $arg ) { 
      $self->{verbose}= $arg ;  
   }
}  


 
sub _make_mfetch_prefix_command {  
   my ($self, $f ) = @_ ; 

   my $mfetch = $self->executable;  
   my $options = $self->options;
   if (defined($options)) { 
     unless ($options =~ /^-/ ) {
       $options = '-' . $options; 
     }
   }  

   my $command = "$mfetch "; 

  # case 1 : we want to fetch entries with -f flag, ie 
  # mfetch -f "acc Taxon " ABCDE890.1    
  # user should only submit field names,NOT " -f Taxon " ... 

  if ( defined $f && ref($f)=~m/ARRAY/) {
    my @fields = @$f ;   
    if ( scalar(@fields) > 0 ) {   
      my $f = join (" ",@fields) ;    
      # remove -f flag if user has submitted it 
      $f=~s/-f//g; 
      # put 'acc' field at the beginning of the string and remove redundancy 
      # and that it's there as well.

      if ( $f=~m/acc/) {  
        $f =~s/acc//g; 
      }
      $f = "acc " .$f; 
      $command .= " -f \"$f \" "      ;
    }
  } 
   

  if (defined $options){
    $command .= "$options ";
  }
  if ($self->verbose ) { 
    print "PREFIX COMMAND $command\n"  
  }
  return $command ; 
}  

1;
