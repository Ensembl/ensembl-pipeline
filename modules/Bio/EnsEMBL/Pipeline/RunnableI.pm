#
# Interface for running programs
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableI

=head1 SYNOPSIS

=head1 DESCRIPTION

Interface for running external programs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableI;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Programs;
use Bio::EnsEMBL::Root;
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::Config::General;

@ISA = qw(Bio::EnsEMBL::Root);


=head1 ABSTRACT METHODS

These methods need to be implemented in any
module which implements
C<Bio::EnsEMBL::Pipeline::RunnableI>.

=head2 run

    $self->run();

Actually runs the analysis programs.  If the
analysis has fails, it should throw an
exception.  It should also remove any temporary
files created (before throwing the exception!).

=head2 output

    @output = $self->output();

Return a list of objects created by the analysis
run (eg: C<Bio::EnsEMBL::FeaturePair> objects).

=cut

sub run {
    my ($self) = @_;

    $self->throw("run not implemented");
}

sub output {
    my ($self) = @_;

		if (!defined($self->{_output})) {
			$self->{_output} = [];
		}

		return @{$self->{_output}};
}

sub locate_executable {
    my ($self, $binary) = @_;

    if ($binary) {
        Bio::EnsEMBL::Analysis::Programs->import($binary);
        return $Bio::EnsEMBL::Analysis::Programs::Program_Paths { $binary };
    }
}

sub filename {
    my ($self, $filename) = @_;
    $self->{_filename} = $filename if ($filename);
    return $self->{_filename};
}

sub results {
    my ($self, $results) = @_;
    $self->{_results} = $results if ($results);
    return $self->{_results};
}

sub queryname {
    my ($self, $queryname) = @_;
    $self->{_queryname} = $queryname if ($queryname);
    return $self->{_queryname};
}


=head2 protect

    Title   :   protect
    Usage   :   $obj->protect('.masked', '.p');
    Function:   Protects files with suffix from deletion when execution ends
    Args    :   File suffixes

=cut

sub protect {
    my ($self, @filename) =@_;
    
    if (!defined($self->{_protected})) {
      $self->{_protected} = [];
    }

    push (@{$self->{_protected}}, @filename) if (@filename);

    return @{$self->{_protected}};
}

=head2 threshold

    Title   :   threshold
    Usage   :   $obj->threshold($value);
    Function:   Get/set method for threshold score required for outputting
                Feature/FeaturePair
    Args    :   Optional value (depends on type of Analysis)

=cut

sub threshold {
  my ($self, $value) = @_;
 
  
  if (defined ($value)) {
    $self->{_threshold} = $value;
  }
    
  return  $self->{'_threshold'};
}


sub threshold_type {
    my($self, $value) = @_;

    if (defined($value)) {
			$self->{'_threshold_type'} = $value;
    }

    return $self->{'_threshold_type'};
}


=head2 workdir

    Title   :   workdir
    Usage   :   $obj->workdir('~humpub/temp');
    Function:   Get/set method for the location of a directory to contain temp files
    Args    :   File path (optional)

=cut

sub workdir {
    my ($self, $directory) = @_;
		
    if (!defined($self->{_workdir})) {
        if ($directory) {
            mkdir ($directory, '777') unless (-d $directory);
            $self->throw ("$directory doesn't exist\n") unless (-d $directory);
            $self->{_workdir} = $directory;
}	elsif ($PIPELINE_WORK_DIR)	{
            $self->{_workdir}= $PIPELINE_WORK_DIR;
        } else {
            $self->{_workdir} = '/tmp';
        }
    }
    return $self->{_workdir};
}


sub get_tmp_file {
    my ($self,$dir,$stub,$ext) = @_;

   
    if ($dir !~ /\/$/) {
        $dir = $dir . "/";
    }

		# This is not good

    my $num = int(rand(100000));
    my $file = $dir . $stub . "." . $num . "." . $ext;
    #print STDERR "have filename ".$file."\n";
    while (-e $file) {
        $num = int(rand(100000));
        $file = $dir.$stub . "." . $num . "." . $ext;
    }
    #print STDERR "returning filename ".$file."\n";
    return $file;
}
   

sub write_sequence_to_file {
    my ($self, $seqobj) = @_;
  
    if (!defined($seqobj)) {
	$self->throw("Must enter a Bio::Seq or a Bio::PrimarySeq object to the write_sequence_to_file");
    }
    if (!$seqobj->isa("Bio::Seq") && !$seqobj->isa("Bio::PrimarySeq")) {
        $self->throw("Must enter a Bio::Seq or a Bio::PrimarySeq object to the write_sequence_to_file. Currently [$seqobj]");
    }

    my $file      = $self->get_tmp_file($self->workdir,"seq","fa");
    my $clone_out = Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
      
    $clone_out->write_seq($seqobj);

    $self->file($file);
    return $file;
}

sub file {
    my ($self,$arg) = @_;

    if (!defined($self->{_files})) {
	$self->{_files} = [];
    }
    if (defined($arg)) {
	push(@{$self->{_files}},$arg);
    }
    return @{$self->{_files}};
}

sub deletefiles {
    my ($self) = @_;

    foreach my $result ($self->file) {

        next unless -e $result;

        my $protected = undef; #flag for match found in $protected

        foreach my $suffix ($self->protect) {
            $protected = 'true' if ($result eq $self->filename.$suffix);
        }
        unless ($protected)
        {
            unlink ($result) or $self->throw ("Couldn't delete $result :$!");
        }
    }
}

sub checkdir {
    my ($self) = @_;

    my $spacelimit = 0.01;

    my $dir = $self->workdir;

    $self->throw("Not enough disk space ($spacelimit required):$!\n") 
                        unless ($self->diskspace($dir, $spacelimit));

    chdir ($dir) or $self->throw("Cannot change to directory $dir ($!)\n");

}

sub diskspace {
    my ($self, $dir, $limit) =@_;
    my $block_size; #could be used where block size != 512 ?
    my $Gb = 1024 ** 3;
    

    open DF, "df $dir |" or $self->throw ("Can't open 'df' pipe ($!)\n");

    while (<DF>) {
			if ($block_size)  {
				my @L = split;
				my $space_in_Gb = $L[3] * 512 / $Gb;
				return 0 if ($space_in_Gb < $limit);
				return 1;
			} else {
				($block_size) = /(\d+).+blocks/i
					|| $self->throw ("Can't determine block size from:\n$_");
			}
    }
    close DF || $self->throw("Error from 'df' : $!\n");
}


# FeaturePair methods:

sub create_FeaturePair {
    my ($self, $feat1, $feat2) = @_;
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => $feat2->{db},
                            -db_version      => $feat2->{db_version},
                            -program         => $feat2->{program},
                            -program_version => $feat2->{p_version},
                            -gff_source      => $feat2->{source},
                            -gff_feature     => $feat2->{primary},
                            -logic_name      => $feat2->{logic_name} );
    
    #create and fill Bio::EnsEMBL::Seqfeature objects
    my $seqfeature1 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname        => $feat1->{name},
                            -start          => $feat1->{start},
                            -end            => $feat1->{end},
                            -strand         => $feat1->{strand},
                            -score          => $feat1->{score},
                            -percent_id     => $feat1->{percent},
                            -p_value        => $feat1->{p},
                            -analysis       => $analysis_obj);
    
    my $seqfeature2 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname        => $feat2->{name},
                            -start          => $feat2->{start},
                            -end            => $feat2->{end},
                            -strand         => $feat2->{strand},
                            -score          => $feat2->{score},
                            -percent_id     => $feat2->{percent},
                            -p_value        => $feat2->{p},
                            -analysis       => $analysis_obj);
    #create featurepair
    my $fp = Bio::EnsEMBL::FeaturePair->new  (  -feature1 => $seqfeature1,
                                                -feature2 => $seqfeature2 ) ;

    #print "Feature pair " . $fp->gffstring . "\n";

    $self->growfplist($fp);                             

    return $fp;
}

sub create_Repeat {
    my ($self, $feat1, $feat2) = @_;
    
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => $feat2->{db},
                            -db_version      => $feat2->{db_version},
                            -program         => $feat2->{program},
                            -program_version => $feat2->{p_version},
                            -gff_source      => $feat2->{source},
                            -gff_feature     => $feat2->{primary},
                            -logic_name      => $feat2->{logic_name});
    
    my $rc = Bio::EnsEMBL::RepeatConsensus->new
                        (   -seqname        => $feat1->{name},
                            -start          => $feat1->{start},
                            -end            => $feat1->{end},
                            -strand         => $feat1->{strand},
                            -score          => $feat1->{score},
                            -source_tag     => $feat1->{source},
                            -primary_tag    => $feat1->{primary},
                            -percent_id     => $feat1->{percent},
                            -p_value        => $feat1->{p},
                            -analysis       => $analysis_obj);
    
    my $f = new Bio::EnsEMBL::RepeatFeature
                        (   -seqname        => $feat2->{name},
                            -start          => $feat2->{start},
                            -end            => $feat2->{end},
                            -strand         => $feat2->{strand},
                            -score          => $feat2->{score},
                            -source_tag     => $feat2->{source},
                            -primary_tag    => $feat2->{primary},
                            -percent_id     => $feat2->{percent},
                            -p_value        => $feat2->{p},
                            -analysis       => $analysis_obj);


    #$self->growfplist($fp);                             

}

sub growfplist {
    my ($self, $fp) =@_;    
    #load fp onto array using command _grow_fplist
    push(@{$self->{'_fplist'}}, $fp);
}

sub shrinkfplist {
    my ($self, $fp) =@_;    
    #load fp onto array using command _grow_fplist
    return pop(@{$self->{'_fplist'}});
}


sub find_executable {
  my ($self,$name) = @_;

  my $bindir = $BIN_DIR || undef;

  if (-x $name) {
    return $name;
  } elsif ($bindir && -x "$bindir/$name") {
    return "$bindir/$name";
  } else {
    eval {
      $name = $self->locate_executable($name);
    };
    if ($@) {
      $self->throw("Can't find executable [$name]");
    }
    return $name;
  }
}

sub writefile {
    my ($self, $seqobj, $seqfilename) = @_;
 
  if (defined($seqobj)) {
    $seqfilename = 'filename' unless ($seqfilename);
    #print "Writing sequence to ".$self->$seqfilename()."\n";
    #create Bio::SeqIO object and save to file
    my $clone_out = Bio::SeqIO->new(-file => ">".$self->$seqfilename(), '-format' => 'Fasta')

      or $self->throw("Can't create new Bio::SeqIO from ".$self->$seqfilename().":$!\n");

    $clone_out->write_seq($self->$seqobj())
      or $self->throw("Couldn't write to file ".$self->$seqfilename().":$!");
      $self->file($seqfilename);


  } else {
    #print "Writing sequence to ".$self->filename."\n";
    #create Bio::SeqIO object and save to file
    my $clone_out = Bio::SeqIO->new(-file => ">".$self->filename , '-format' => 'Fasta')
      or $self->throw("Can't create new Bio::SeqIO from ".$self->filename.":$!\n");

    # This is bad.  The subclass has the query method not this interface.
    $clone_out->write_seq($self->query)  or $self->throw("Couldn't write to file ".$self->filename.":$!");
    $self->file($self->filename);
  }
}

sub find_file {
  my ($self,$name) = @_;

  my $datadir = $DATA_DIR || undef;
  my $libdir  = $LIB_DIR  || undef;

  my $full_name;

  if (-e $name) {
    return $name;
    
  } elsif ($datadir && -e ($full_name = "$datadir/$name")) {
    return $full_name;
  } elsif ($libdir && -e ($full_name = "$libdir/$name")) {
    return $full_name;
  } else {
    $self->throw("Can't find file [$name]");
  }
}

=head2 options

    Title   :   options
    Usage   :   $obj->options('M=500');
    Function:   Get/set method for binary options
                e.g. M=? for e-PCR, B=? for blast
    Args    :   option string

=cut

sub options {
    my ($self, $args) = @_;
    if ($args)
    {
        $self->{'_options'} = $args ;
    }
    return $self->{'_options'};
}


=head2 create_SimpleFeature

=cut

sub create_SimpleFeature {
    my ($self, $feat) = @_;

    my $analysis_obj = Bio::EnsEMBL::Analysis->new(
        -db              => undef,
        -db_version      => undef,
        -program         => $feat->{'program'},
        -program_version => $feat->{'program_version'},
        -gff_source      => $feat->{'source'},
        -gff_feature     => $feat->{'primary'}
    );

    my $sf = Bio::EnsEMBL::SimpleFeature->new(
        -seqname     => $feat->{'name'},
        -start       => $feat->{'start'},
        -end         => $feat->{'end'},
        -strand      => $feat->{'strand'},
        -score       => $feat->{'score'},
        -analysis    => $analysis_obj
    );  

    # display_label must be a null string, and not undef
    # can't be set above as it is not known to SeqFeature
    # (SimpleFeature->new uses SeqFeature->new)
    $sf->display_label($feat->{'hit'});

    if ($sf) {
	$sf->validate();

	# add to _sflist
	push(@{$self->{'_sflist'}}, $sf);
    }
}


=head2 trunc_float_3

  Arg  1      float $value
              floating point number
  Function    Truncates a float to 3 places of decimal
  Returntype  float
  Exceptions  none
  Caller      Anyone!

=cut

sub trunc_float_3 {
    my ($self, $arg) = @_;

    # deal only with valid numbers
    # and only need cases of the form [+/-]xx.yyyyy
    return $arg unless $arg =~ /^[+-]?\d*\.\d+$/;

    return 0.001 * int (1000 * $arg);
}

1;
