package Bio::EnsEMBL::Pipeline::Tools::Embl;

### embl ###

use strict;
#no warnings;

sub new{
    my ($proto) = shift;
    my $class = ref($proto) || $proto;
    my $p;
    if(ref($_[0]) eq "HASH"){($p) = @_ ;}
    else{ $p = {@_} ;}
    my $self = {
        _debug => 0,
        _acc   => undef,
        _fh    => undef,
        map { "_".substr($_,1) => $p->{$_} } keys(%$p)
    };
    bless ($self, $class);
    if($self->{'_debug'} > 2){
        require 'Data/Dumper.pm';
        print Data::Dumper->Dump([$self],['self']);
    }
    return $self;
}
sub _fh{
    my $self = shift;
    if(my $fh = shift){
        map { delete $self->{$_} } keys(%{$self});
        die "Not a filehandle [_fh]" unless ref($fh) eq "GLOB" || ref($fh) eq 'IO::String';
        $self->{'_fh'} = $fh;
    }
    return $self->{'_fh'};
}

sub parse{
    my ($self,$string) = @_;
    if($string){$self->_parse_record($string)}
    elsif($self->_fh){
        local $/ = "//\n";
        my $FH = $self->_fh();
        my $single_check = 0;
        while(<$FH>){
            next unless $single_check < 1;
            $self->_parse_record($_);
            $single_check++;
        }
    }
}
sub _parse_record{
    my ($self,$string) = @_;
    $self->clean();
    $string =~ s/^SQ\s{3}(.+;)(.+)/$self->sequence($2,$1)/ems;
    $string =~ s/^([A-Z]{2})\s{3}(.+\n)/$self->{"_" . $1} .= $2/egm;
    if($self->{'_FH'}){
        chop $self->{'_FT'};
        my @keys = split(/\n\b/, $self->{'_FT'});
        delete $self->{'_FT'};
        foreach(@keys){
            my @qualifiers = split(/\s+\//, $_);
            my ($key,$location) = split(/\s+/,shift(@qualifiers));
#	        my $hash = { map { ($a, $b) = split('=') } @qualifiers }; # unfortunately not, [non unique $a!] :(
#	        $self->{'_FT'}->{$key}->{$location} = $hash || {};
            my $hash;
            foreach(@qualifiers){
                my ($qualifier, $value) = (split(/=/),'')[0,1];
                $value =~ s/^\"|\"$//g;
                $value =~ s/(\s?)\n\s+(\B)/$1/g;
                push(@{$self->{'_FT'}->{$key}->{$location}->{$qualifier}}, $value);
            }
        }
    }
    else{
        $self->_non_EMBL(1);
    }
    #   print "found accession: <" . join(" ", @{$self->accession()}) .">\n";
    #   print ", sequence version: <". $self->sequence_version() . ">";
    #   print ", description: <" . $self->description() . ">\n";
    #   print " * fasta:\n" . $self->fasta() . "\n";
}

sub feature_table{
    my $self = shift;
    return $self->{'_FT'};
}
sub comment{
    my $self = shift;
    return $self->{1}->{'_CC'};
}
sub description{
    my $self = shift;
    $self->{'_DE'} =~ s/\n/ /g;
    $self->{'_DE'} =~ s/\.\s$//;
    return $self->{'_DE'};
}
sub date{
    my $self = shift;
    unless( ref($self->{'_DT'}) eq "ARRAY"){
        $self->{'_DT'} = [ split("\n",$self->{'_DT'}) ];
    }
    return $self->{'_DT'};
}
sub sequence{
    my ($self, $sequence, $sequence_line) = @_;
    if(@_){
        $self->sequence_line($sequence_line);
        $self->{'_sequence'} = $sequence;
        $self->{'_sequence'} =~ s/[\s0-9\/]//g;
        return '';
    }
    return $self->{'_sequence'};
}
sub sequence_line{
    my $self = shift;
    $self->{'_SQ'} = shift if @_;
    return $self->{'_SQ'};
}
sub _non_EMBL{
    my $self = shift;
    $self->{'_non_EMBL'} = shift if @_;
    return $self->{'_non_EMBL'} ? 1 : 0;
}
sub accession{
    my $self = shift;
    $self->{'_AC'} =~ s/[;]//g;
    unless(ref($self->{'_AC'}) eq "ARRAY"){
        $self->{'_AC'} = [ split(/\s+/, $self->{'_AC'}) ];
    }
    return $self->{'_AC'};
}
sub sequence_version{
    my $self = shift;
    $self->{'_SV'} =~ s/\n//g if $self->{'_SV'};
    return $self->{'_SV'};
}
sub ox{
    # what's OX? a chop'ed OXO??
    my $self = shift;
    $self->{'_OX'} =~ s/[;\n]//g;
    return $self->{'_OX'};
}
sub os{
    my $self = shift;
    $self->{'_OS'} =~ s/[\.$|\n]//g;
    return $self->{'_OS'};
}
sub keywords{
    my $self = shift;
    unless(ref($self->{'_KW'}) eq "ARRAY"){
        $self->{'_KW'} =~ s/\n|\.$//g;
        $self->{'_KW'} = [ split(/;\s?/, $self->{'_KW'}) ];
    }
    return $self->{'_KW'};
}
sub fasta{
    my $self = shift;
    my ($length, $fasta) = (60,undef);
    my $seq = $self->sequence();
    unless($fasta = [ $seq =~ /\G.{1,$length}/g ]){
        $fasta = [ $seq ];
    }
    return ">" . $self->description() . "\n" . join("\n",@{$fasta});
}
sub clean{
    my $self = shift;
    map { delete $self->{$_} } keys(%{$self});
}
sub DESTROY{## DESTRUCTOR if needed
    my ($self) = @_;
    ## Clear/Correct Class Variables
}
# +-------------------------------------------------+
# | useful methods                                  |
# |                                                 |
# +-------------------------------------------------+
sub taxon{
    my ($self, $obj) = @_;
    my $ft = $self->feature_table();
    my $ret;
    if(ref($ft) eq "HASH"){
        foreach my $location(keys %{$ft->{'source'}}){
            foreach(@{$ft->{'source'}->{$location}->{'db_xref'}}){
                return $1 if /taxon:(\d+)/;
            }
        }
        }else{
        return $1 if $self->ox =~ /NCBI_TaxID=(\d+)/;
    }
    return undef;
}
sub seq_length{
    my ($self, $obj) = @_;
    my $arr = [ split(";", $self->sequence_line) ];
    my $first = shift @{$arr};
    my $length = $1 if $first =~ /(\d+)/;
    return $length;
}
sub which_database{
    my ($self, $obj, $ret) = @_;
    if($self->_non_EMBL){
        $ret = 'SWISSPROT';
        my $date = shift( @{$self->date} );
        $ret = 'TrEMBL' if $date =~ /TrEMBL/;
        }else{
        $ret = 'EMBL'; # default assumption
    }
    return $ret;
}
1;

=head1 

=head2 DESCRIPTION

Designed to only read one embl entry @ a time, as they can be big.

=cut
