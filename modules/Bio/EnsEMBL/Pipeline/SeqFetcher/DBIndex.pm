package Bio::EnsEMBL::Pipeline::SeqFetcher::DBIndex;

use strict;
use vars qw(@ISA);

use FileHandle;

use Bio::EnsEMBL::Root;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
    my($class, @args) = @_;

    my $self = $class->SUPER::new(@args);

    #Should we be able to search multiple databases?

    my( $database ) =  $self->_rearrange([qw(DATABASE)], @args);
    
    # Store any parameters passed
    $self->database($database);

    if ($self->database) {
      $self->read_config;
      $self->read_fileids;
      $self->read_header;
    }

    return $self;
}

sub read_config {
  my ($self) =  @_;

    if (!defined($self->database)) {
	$self->throw("Database not set in index. Can't read config.dat");
    }

  my $configfile = $self->database . "/config.dat";

  if (! -e $configfile) {
    $self->throw("No config file [$configfile]. Can't read namespace");
  }

  open(IN,"<$configfile") || $self->throw("Can't open configfile [$configfile]");

  while (<IN>) {
    chomp;

    my ($key,$value) = split(/\t/,$_);

    if ($key eq "namespace") {
      $self->namespace($value);
    } else {
      $self->throw("Unrecognised key value pair in config.dat [$key][$value]");
    }
  }
  close(IN);
}

sub read_fileids {
    my ($self) = @_;

    if (!defined($self->database)) {
	$self->throw("Database not set in index. Can't read fileids");
    }

    my $fileid_file = $self->database . "/fileids.dat";

    if (! -e $fileid_file) {
	$self->throw("Fileid file [$fileid_file] doesn't exist for database [" . 
		     $self->database . "\n");
    }

    open(FILEID,"<$fileid_file") || $self->throw("Can't open fileid file [$fileid_file");

    while (<FILEID>) {
	chomp;
	my ($id,$name,$size) = split(/\t/);

	if (!defined($self->{_fileid}{$id})) {
	    my $fh = new FileHandle("<$name");
	    $self->{_fileid}{$id}   = $fh;
	    $self->{_dbfile}{$name} = $id;
            $self->{_size}{$id}     = $size; 
	} else {
	    $self->throw("Something's wrong - there are two identical file ids $id");
	}
    }
    
    close(FILEID);
}



sub read_header {
    my ($self) = @_;

    my $file = $self->primary_index_file;;

    open (IN,"<$file") || $self->throw("Can't read index [$file]");

    my $fh = \*IN;

    my $num_fields;
    my $record_width;
   
    my @field_width;

    my $i = 0;
    
    sysread($fh,$record_width,4);
    sysread($fh,$num_fields,4);

    print "Width and fields $record_width:$num_fields:\n";

    $record_width =~ s/ //g;
    $num_fields   =~ s/ //g;

    $record_width = $record_width * 1;
    $num_fields   = $num_fields * 1;

    $self->{_num_fields} = $num_fields;
    $self->record_size($record_width);


    while ($i < $num_fields) {
      sysread($fh,$field_width[$i],5);

       print "Field width for $i is $field_width[$i]\n";
       $i++;
    }
    $self->{_start_pos}    = tell($fh);
    $self->{_field_widths} =  [];

    my $tmp;

    sysseek($fh,$self->{_start_pos},0);
    sysread($fh,$tmp,$field_width[0]);
    print "Tmp $tmp\n";

    sysseek($fh,$self->{_start_pos} + 1000*$self->record_size,0);
    sysread($fh,$tmp,$field_width[0]);
    print "Tmp $tmp\n";
    push(@{$self->{_field_widths}},@field_width);

  }

sub read_record {
  my ($self,$fh,$pos) = @_;

  sysseek($fh,$pos,0);

  my @field_width =  @{$self->{_field_widths}};

  my $id;
  my $fileid;
  my $pos;
  my $length;
    
  sysread($fh,$id,    $field_width[0]);
  sysread($fh,$fileid,$field_width[1]);
  sysread($fh,$pos,   $field_width[2]);
  sysread($fh,$length,$field_width[3]);

  return ($id,$fileid,$pos,$length);


}

sub get_Seqs_by_id_array {
    my ($self,@ids) = @_;
    
    my %entry;

    foreach my $id (@ids) {
	$entry{$id} = $self->get_Seq_by_id($id);
    }
    return %entry;
}



sub get_Seq_by_id {
    my ($self,$id) = @_;

    my $indexfh = $self->primary_index_filehandle;

    my $recsize = $self->record_size;

    sysseek ($indexfh,0,2);

    my $filesize = (tell $indexfh);

    my $end = ($filesize-$self->{_start_pos})/$recsize;


    print "File size is $filesize $end " . $self->{_start_pos} . "\n";

    my ($id,$fileid,$fpos,$length) = $self->find_entry($indexfh,$recsize,0,$end,$id);

    $id      =~ s/ //g;
    $fileid  =~ s/ //g;
    $fpos    =~ s/ //g;
    $length  =~ s/ //g;
    $length  =~ s/ +$//;

    if (!defined($id)) {
      return;
    }

    my $fh = $self->get_filehandle_by_fileid($fileid);

    return $self->get_entry($fh,$fpos,$length);
}

sub get_entry {
  my ($self,$fh,$pos,$length) = @_;

  
  my $entry;

  sysseek ($fh,$pos,0);
  sysread($fh,$entry,$length);

  return $entry;
}

sub find_entry {
    my ($self,$fh,$record_size,$start,$end,$id) = @_;

    my $mid = int(($end+1+$start)/2);

    my $pos = ($mid -1)* $record_size;

    sysseek($fh,$self->{_start_pos}+ $pos ,0);

    my @field_width =  @{$self->{_field_widths}};

    my $entryid;

    sysread($fh,$entryid,$field_width[0]);

    $entryid =~ s/ //g;

    my ($first,$second) = sort { $a cmp $b} ($id,$entryid);
#    print "$mid $pos $start $end $id $entryid\n";
    if ($id eq $entryid) {
      my $pos;
      my $length;
      my $fileid;

      sysread($fh,$fileid,$field_width[1]);
      sysread($fh,$pos,$field_width[2]);
      sysread($fh,$length,$field_width[3]);

      return ($id,$fileid,$pos,$length);

    } elsif ($first eq $id) {

      if ($end-$start <= 1) {
	return;
      }
	my $end = $mid;

      $self->find_entry($fh,$record_size,$start,$end,$id);

    } elsif ($second eq $id ) {

      if ($end-$start <= 1) {
	return;
      }

      $start = $mid;
      
      $self->find_entry($fh,$record_size,$start,$end,$id);
    } else {
	print "No match\n";
    }

 }   

sub make_index {

  my ($self,$dbname,$format,@files) = @_;;
    
  if (!defined(@files)) {
    $self->throw("Must enter an array of filenames to index");
  }

  if (!defined($dbname)) {
    $self->throw("Must enter an index name for your files");
  }

  foreach my $file (@files) {
    if (! -e $file) {
      $self->throw("Can't index file [$file] as it doesn't exist");
    }
  }

  print "Database is $dbname\n";

  $self->database($dbname);
  $self->_make_indexdir;
  

  # Check the available disk space


  $self->_make_BIOINDEX_file;
  $self->_make_config_file;
  $self->_make_fileid_file(@files);

    
  # Finally lets index
  foreach my $file (@files) {
    $self->_index_file($file);
  }
 
  $self->write_primary_index;

  # print some stats
    
}

sub _index_file {
    my ($self,$file) = @_;

    open(FILE,"<$file") || $self->throw("Can't open file [$file]");

    my $recstart = 0;
    my $fileid = $self->get_fileid_by_filename($file);
    my $found = 0;
    my $id;
    my $count;

    LINE: while (<FILE>) {
       if (/^>/) {
           $count++;

           my $newid = $_;
           $newid =~ s/^(\S+).*/$1/;
           chomp($newid);

           my $begin = tell(FILE) - length( $_ );
           my $length = ($begin - $recstart);
    
        if ($found) {
	if (!defined($id)) {
	    $self->throw("No id defined for sequence");
	}
	if (!defined($fileid)) {
	    $self->throw("No fileid defined for file $file");
	}
	if (!defined($recstart)) {
	    $self->throw("No position defined for " . $id . "\n");
	}
	if (!defined($length)) {
	    $self->throw("No length defined for " . $id . "\n");
	} 
        $self->_add_id_position($id,$recstart,$fileid,$length);
         } else {
             $found = 1;
         }
        $recstart = $begin;
        $id       = $newid;

	if ($count%1000 == 0) {
	    print "Indexed $count ids. Time " . time . "\n";
	}
	 }
     }
    close(FILE);
}

sub write_primary_index {
    my ($self) = @_;

    my @ids = sort keys %{$self->{_id}};

    print "Number of ids = " . scalar(@ids) . "\n";

    open (INDEX,">" . $self->primary_index_file);

    my $recordlength = $self->{_maxidlength} +
	               $self->{_maxfileidlength} + 
	               $self->{_maxposlength} +
   		       $self->{_maxlengthlength} + 6;
	
    foreach my $id (@ids) {
	if (!defined($self->{_id}{$id}{_fileid})) {
	    $self->throw("No fileid for $id\n");
	}
	if (!defined($self->{_id}{$id}{_pos})) {
	    $self->throw("No position for $id\n");
	}
	if (!defined($self->{_id}{$id}{_length})) {
	    $self->throw("No length for $id");
	}

	my $record =  $id . "\t" . 
	    $self->{_id}{$id}{_fileid} . "\t" .
	    $self->{_id}{$id}{_pos} .    "\t" .
	    $self->{_id}{$id}{_length};

	if (length($record) < $recordlength) {
	    my $pad = ' ' x ($recordlength - length($record));
	    $record = $record . $pad;
	}
	print INDEX $record . "\n";

    }
    close(INDEX);
}
	    
sub primary_index_file {
    my ($self) = @_;

    return $self->database . "/key_" . $self->namespace . ".key";
}

sub primary_index_filehandle {
    my ($self) = @_;

    if (!defined ($self->{_primary_index_handle})) {
	$self->{_primary_index_handle} = new FileHandle("<" . $self->primary_index_file);
    }
    return $self->{_primary_index_handle};
}

sub _add_id_position {
    my ($self,$id,$pos,$fileid,$length) = @_;

    if (!defined($id)) {
	$self->throw("No id defined. Can't add id position");
    }
    if (!defined($pos)) {
	$self->throw("No position defined. Can't add id position");
    }
    if (!defined($fileid)) {
	$self->throw("No fileid defined. Can't add id position");
    }
    if (!defined($length) || $length <= 0) {
	$self->throw("No length defined or <= 0 [$length]. Can't add id position");
    }

    $self->{_id}{$id}{_pos}    = $pos;
    $self->{_id}{$id}{_length} = $length;
    $self->{_id}{$id}{_fileid} = $fileid;

    if (length($id) >= $self->{_maxidlength}) {
	$self->{_maxidlength} = length($id);
    }

    if (length($pos) >= $self->{_maxposlength}) {
	$self->{_maxposlength} = length($pos);
    }

    if (length($length >= $self->{_maxlengthlength})) {
	$self->{_maxlengthlength} = length($length);
    }
}


sub _make_indexdir {
    my ($self) = @_;

    if (! -e $self->database) {
	mkdir $self->database,0755;
    } else {
	$self->throw("Index directory " . $self->database . " already exists. Exiting\n");
    }
}

sub _make_BIOINDEX_file {
    my ($self) = @_;
}

sub _make_config_file {
    my ($self) = @_;
}

sub _make_fileid_file {
    my ($self,@files) = @_;

    my $dir = $self->database;
    
    if (! -d $dir) {
	$self->throw("[$dir] is not a directory.  Can't write fileid file");
    }

    my $fileid_name = $dir . "/fileid";
    
    
    open(FILEID,">$fileid_name") || $self->throw("Can't write fileid file [$fileid_name]");

    my $count = 1;
    foreach my $file (@files) {

	print FILEID $count . "\t" . $file . "\n";
        $count++;
    }

    close(FILEID);

    $self->read_fileids;
}

sub get_fileid_by_filename {
    my ($self,$file) = @_;
    
    if (!defined($self->{_dbfile})) {
	$self->throw("No file to fileid mapping present.  Has the fileid file been read?");
    }

    
    return $self->{_dbfile}{$file};
}

sub get_filehandle_by_fileid {
    my ($self,$fileid) = @_;

    if (!defined($self->{_fileid}{$fileid})) {
	$self->throw("ERROR: undefined fileid in index [$fileid]");
    }

    return $self->{_fileid}{$fileid};
}


sub database {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_database} = $arg;
    }
    return $self->{_database};

}

sub record_size {
    my ($self,$arg) = @_;

    if (defined($arg)) {
      $self->{_record_size} = $arg;
    }
    return $self->{_record_size};
}

sub namespace {
  my ($self,$arg) =  @_;

  if (defined($arg)) {
    $self->{_namespace} =  $arg;
  }
  return $self->{_namespace};
}

1;

	

    



