
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::BlastDB

=head1 SYNOPSIS

=head2 Methods:

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BlastDB;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    
           
    my ($sequences,$dbfile, $type)  = $self->_rearrange([qw(SEQUENCES
																														DBFILE
																														TYPE
																													 )], @args);
		
		
		if (!defined($dbfile) && !defined($sequences)) {
			$self->throw("No dbfile or sequence array ref entered for indexing");
		}
		if (!defined($type))  { 
			$self->throw("Not database type entered");
		}

		$self->dbfile($dbfile)        if defined($dbfile);
		$self->sequences($sequences)  if defined($sequences);
		$self->type($type)            if defined($type);

		return $self;

}
sub sequences {
	my ($self,$sequences) = @_;

	if (!defined($self->{_sequences})) {
		$self->{_sequences} = [];
	}
	if (defined($sequences)) {
		if (ref($sequences) eq "ARRAY") {
			push(@{$self->{_sequences}},@$sequences);
		} else {
			$self->throw("Argument must be an array ref . Currently is [$sequences]");
		}
	}
	return @{$self->{_sequences}};
}

sub run {
	my ($self) = @_;

	if (!defined($self->dbfile)) {
		my $seqfile = $self->make_seqfile;
		$self->dbfile($seqfile);
	}

	my $seqfile = $self->dbfile;

	if (! -e $seqfile) {
		$self->throw("Database [$seqfile] doesn't exist. Can't index");
	}

	if ($self->type eq 'PROTEIN') {
		my $status = system("setdb $seqfile");
	} elsif ($self->type eq 'DNA') {
    my $status = system("pressdb $seqfile");
	}

	return 1;
}


sub make_seqfile {
	my ($self) = @_;
	
	my $tmpdir = '/tmp';

	my $blastfile = $self->get_tmp_file($tmpdir,'blast','fa');

	my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
                               -file    => ">$blastfile");

	foreach my $seq ($self->sequences) {
		$seqio->write_seq($seq);
	}

	close($seqio->_filehandle);

	return $blastfile;
}

sub dbname {
	my ($self) = @_;

	if (!defined($self->dbfile)) {
		$self->throw("No database file defined - can't removed index files");
	}

	my $dbname = $self->dbfile;

	$dbname =~ s/.*\/(.*?)/$1/;

	return $dbname;
}


sub remove_index_files {
	my ($self) = @_;

	if (!defined($self->dbfile)) {
		$self->throw("No database file defined - can't removed index files");
	}

	if ($self->type eq 'DNA') {
		unlink $self->dbfile . ".csq";
		unlink $self->dbfile . ".nhd";
		unlink $self->dbfile . ".ntb";
	} elsif ($self->type eq 'PROTEIN') {
		unlink $self->dbfile . ".bsq";
		unlink $self->dbfile . ".ahd";
		unlink $self->dbfile . ".atb";
	} else {
		$self->throw("Type [" . $self->type . "] not recognised");
	}
}
#################
# get/set methods 
#################

=head2 dbfile

 Title   : dbfile
 Usage   : $obj->dbfile($newval)
 Function: 
 Example : 
 Returns : value of dbfile
 Args    : newvalue (optional)


=cut

sub dbfile{
   my ($obj,$value) = @_;
   if( defined $value) {
		 if ($value !~ /^\//) {
			 my $pwd = `pwd`;
			 chomp($pwd);
			 $value = $pwd . "/" . $value;
		 }
      $obj->{'dbfile'} = $value;
    }
    return $obj->{'dbfile'};

}
=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: 
 Example : 
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub type{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'type'} = $value;
    }
    return $obj->{'type'};

}

1;
