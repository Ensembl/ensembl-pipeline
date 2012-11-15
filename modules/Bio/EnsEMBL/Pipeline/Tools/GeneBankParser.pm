package Bio::EnsEMBL::Pipeline::Tools::GeneBankParser;

=pod

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::GeneBankParser

=head1 SYNOPSIS

  	
=head1 DESCRIPTION

A collection of subroutines to parse GeneBank files

=head1 METHODS

See subroutines.

=head1 MAINTAINER

$Author: th3 $

=head1 VERSION

$Revision: 1.3 $

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::IO qw/:iterate/;

my $HASH;

sub new {
    my $class = shift;
    my $self = {
        filehandle => shift,
    };
    bless $self, $class;
    if (!defined($self->{'filehandle'})) {
        throw("GeneBankParser requires a valid filehandle to a GeneBank formatted file"); 
    }
    return $self;

}

sub next_entry {
    my $self = shift;
    $HASH = { __feature_index => -1,
              __header         => 1,
            };
    iterate_lines($self->{filehandle}, \&_parse_entry);
    to_ensEMBL();
    return $HASH;
}

sub next_seq {
    my $self = shift;
    $HASH = { __feature_index => -1,
              __header        => 1,
            };
    iterate_lines($self->{filehandle}, \&_parse_seq);
    to_ensEMBL();
    return $HASH;
}

sub _parse_seq {
    my $line = shift;
    chomp($line);
    _get_header($line) if ($HASH->{__header});
    _get_sequence($line) unless ($HASH->{__header});
}

sub _parse_entry {
    my $line = shift;
    chomp($line);
    _get_header($line) if ($HASH->{__header});
    _get_features($line) unless ($HASH->{__header});
    _get_sequence($line) if ($HASH->{__seq});
}

sub _get_header {
    my $line = shift;
    if ($line =~ /^LOCUS\s+(\S+)\s+(\d+)\s+bp\s+(\w+)\s+(\w+)\s+(\w+)\s+(\S+)/i) {
        $HASH->{_locus_id} = $1;
        $HASH->{_length}   = $2;
        $HASH->{_molecule} = $3;
        $HASH->{_tax}      = $5;
        $HASH->{_date}     = $6;
        if ($4 eq 'circular') {
            $HASH->{-is_circular} = 1;
        }
        else {
            $HASH->{-is_circular} = 0;
        }
    }
    elsif ($line =~ /^DEFINITION\s+(\S+.*\S)\s*$/i) {
        $HASH->{-definition} = $1;
    }
    elsif ($line =~ /^ACCESSION\s+(\S+)/i) {
       $HASH->{-accession} = $1;
    }
    elsif ($line =~ /^VERSION\s+\S+\.(\d+)\s+GI:(\d+)/i) {
        $HASH->{-version} = $1;
        $HASH->{-secondary_id} = $2;
    }
    elsif ($line =~ /^DBLINK\s+(\S+.*\S)\s*$/i) {
        $HASH->{_dblink} = $1;
    }
    elsif ($line =~ /^KEYWORDS\s+([^.]*)(\.)?\s*$/i) {
        push(@{$HASH->{_keywords}}, split(/;/, $1)) if ($1 ne '');
        if (!defined $2) {
            $HASH->{__hopen} = 1;
            $HASH->{__current} = '_keywords';
        }
    }
    elsif ($line =~ /^SOURCE\s+(\S+.*\S)\s*$/i) {
        $HASH->{_source} = $1;
    }
    elsif ($line =~ /ORGANISM\s+(\S+.*\S)\s*$/i) {
        $HASH->{__hopen} = 1;
        $HASH->{__current} = '_taxonomy';
    }
    elsif ($line =~ /^COMMENT\s+(\S+.+)/i) {
        $HASH->{__hopen} = 1;
        $HASH->{_comment} = $1;
        $HASH->{__current} = '_comment';
    }
    elsif ($line =~ /^REFERENCE/) {
        delete $HASH->{__hopen} if ($HASH->{__hopen});
    }
    elsif (exists $HASH->{__hopen}) {
        if ($HASH->{__current} eq '_comment') {
            if ($line =~ /^FEATURES/) {
                delete $HASH->{__hopen};
                $HASH->{__header} = 0;
            }
            else {
                $line =~ s/^\s+//;
                $HASH->{$HASH->{__current}} .= ' '.$line;
            }
        }
        else {
            $line =~ /\s+([^.]*)(\.)?\s*$/;
            delete $HASH->{__hopen} if (defined $2);
            push(@{$HASH->{$HASH->{__current}}}, split(/; /, $1)) if (defined $1);
        }
    }
}

sub _get_features {
    my $line = shift;
    if ($line =~ /FEATURES/) {
        $HASH->{__features} = 1;
    }
    if ($line =~ /ORIGIN/) {
        $HASH->{__features} = 0;
        $HASH->{__seq} = 1;
    }
    return unless $HASH->{__features};
    if ($line =~ /(\w+)\s+(<)?(\d+)\.\.(\d+)/) {
        $HASH->{__current} = $1;
        $HASH->{-features}->[++$HASH->{__feature_index}]->{$HASH->{__current}}->{exons}->[0]->{start} = $3;
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{fragment} = 1 if (defined $2);
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{exons}->[0]->{end} = $4;
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{exons}->[0]->{strand} = 1;
    }
    elsif ($line =~ /(\w+)\s+complement\((\d+)\.\.(\d+)/) {
        $HASH->{__current} = $1;
        $HASH->{-features}->[++$HASH->{__feature_index}]->{$HASH->{__current}}->{exons}->[0]->{start} = $2;
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{fragment} = 1 if ($line =~ '<');
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{exons}->[0]->{end} = $3;
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{exons}->[0]->{strand} = -1;
    }
    elsif ($line =~ /(\w+)\s+(complement\()?join\((.+)/) {
        $HASH->{__current} = $1;
        ++$HASH->{__feature_index};
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{fragment} = 1 if ($line =~ '<');
        my $strand = 1;
        $strand = -1 if ($2);
        my $string = $3;
        while ($string =~ /(\d+)\.\.(\d+)/g) {
            push(@{$HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{exons}}, {start => $1, end => $2, strand => $strand});
        }
    }
    elsif ($line =~ /\/(\w+)="?([^"]+)(")?/) {
        if ($1 eq 'db_xref') {
            $2 =~ /([^:]+):(\S+)/;
            $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{db_xref}->{$1} = $2;
        }
        else {
            $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{$1} = $2;
            $HASH->{__fopen} = $1 unless $3;
        }
    }
    elsif (exists $HASH->{__fopen}) {
        $line =~ /(\S.*)(\S)\s*$/;
        $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{$HASH->{__fopen}} .= ' '.$1;
        if ($2 eq '"') {
            delete $HASH->{__fopen};
        }
        else {
            $HASH->{-features}->[$HASH->{__feature_index}]->{$HASH->{__current}}->{$HASH->{__fopen}} .= $2;
        }
    }
}

sub _get_sequence {
    my $line = shift;
    if ($line =~ /ORIGIN/) {
        $HASH->{__seq} = 1;
        return;
    }
    if ($line =~ /\/\//) {
        $HASH->{__seq} = 0;
        $HASH->{__header} = 1;
    }
    return unless $HASH->{__seq};
    $line =~ s/[\d\s]+//g;
    $HASH->{-seq} .= $line;
}

sub to_ensEMBL {
}


1;
