package Bio::EnsEMBL::Pipeline::SubseqFetcher;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);


sub new {

    my ($class, $file) = @_;

    my $self = {};
    bless $self,$class;
    
    open(DATA_HANDLE, "<$file") or $self->throw("Cant open file $file");

    my $data_fh = \*DATA_HANDLE;

    $self->data_file($data_fh);

    my $index_file = "/tmp/subseqfetch.idx";  # Ooh, must make this nicer

    open(INDEX_HANDLE, "+>$index_file") or $self->throw("Cant open file to write index");

    my $index_fh = \*INDEX_HANDLE;

    $self->index_file($index_fh);

    $self->build_index;

    return $self;
}



sub build_index {
    my ($self) = @_;

    my $data_file = *{$self->data_file};
    my $index_file = *{$self->index_file};

    my $offset     = 0;

    while (<$data_file>) {
        print $index_file pack("N", $offset);
        $offset = tell($data_file);
    }
}

# usage: line_with_index($LINE_NUMBER)
# returns line or undef if LINE_NUMBER was out of range
sub line_with_index {
    my ($self, $line_number) = @_;

    my $data_file = $self->data_file;
    my $index_file = $self->index_file;

    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file

    $size = length(pack("N", 0));
    $i_offset = $size * ($line_number-1);
    seek($index_file, $i_offset, 0) or return;
    read($index_file, $entry, $size);
    $d_offset = unpack("N", $entry);
    seek($data_file, $d_offset, 0);
    return scalar(<$data_file>);
}

sub data_file {
    my ($self, $value) = @_;

    if (defined $value) {
	$self->{'_data_file'} = $value;
    } elsif (defined $self->{'_data_file'}) {
	return $self->{'_data_file'};
    }

    return 1;
}

sub index_file {
    my ($self, $value) = @_;

    if (defined $value) {
        $self->{'_index_file'} = $value;
    } elsif (defined $self->{'_index_file'}) {
        return $self->{'_index_file'};
    }

    return 1;
}


sub header_offset {

    return 1;

}

sub nucleotides_per_line {

    return 60;

}

sub subseq {
    my ($self, $start, $end) = @_;

#print "Start $start\tEnd $end\n";

    my $line_offset = $self->header_offset;
    my $line_length = $self->nucleotides_per_line;

    my $start_line = (($start / $line_length) - (($start % $line_length)/$line_length)) + 1 + $line_offset;

    my $end_line;

    if ($end % $line_length != 0){
	$end_line = (($end / $line_length) - (($end % $line_length)/$line_length)) + 1 + $line_offset;
    } else {
	$end_line = ($end / $line_length) + $line_offset;
    }
#print "Start line $start_line\tEnd line $end_line\n";

    my $raw_lines = '';

    for (my $line = $start_line; $line <= $end_line; $line++) {

	$raw_lines .= $self->line_with_index($line);
#print $self->line_with_index($line) . "\n";
    }

    my $raw_start = (($start_line - $line_offset) * $line_length) - ($line_length - 1);
    my $nucleotide_offset = $start - $raw_start;
    my $subseq_length = $end - $start + 1;
#print "Raw start $raw_start\tNucl offset $nucleotide_offset\tLength $subseq_length\n";
    $raw_lines =~ s/\n//;
    my @seq_array = split //, $raw_lines; 

    my $subseq = '';
    for (my $i = $nucleotide_offset + 1; $i <= ($nucleotide_offset + $subseq_length); $i++){
	$subseq .= $seq_array[$i-1];
    }
#print "Subseq $subseq\n";
    return $subseq;
}


sub DESTROY {
    my ($self) = @_;

    unlink $self->index_file;

}

return 1;


