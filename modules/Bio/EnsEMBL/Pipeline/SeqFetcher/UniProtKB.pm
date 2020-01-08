=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::SeqFetcher::UniProtKB -

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::UniProtKB->new(
							     );
    my $seq = $obj->get_Entry_fields($acc);

=head1 DESCRIPTION

  Object to retrieve protein classification

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
package Bio::EnsEMBL::Pipeline::SeqFetcher::UniProtKB;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use LWP::UserAgent;
use Bio::EnsEMBL::IO::Parser::EMBL;
use vars qw(@ISA);
@ISA = qw();

use constant BASE_URL => 'http://www.uniprot.org/uniprot/';

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($options ) = rearrange(['OPTIONS'], @args);

  if (defined $options) {
    $options=~s/^\s+//g;
    $self->options($options);
  }

  return $self; # success - we hope!
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

    throw("No accession input") if (!defined($acc_to_fetch));
    print "Fields to get : " . join ( " " , @$fields )."\n" if $self->{verbose};

    my @queries;
    my %acc_hash;
    if ( ref($acc_to_fetch)=~m/ARRAY/ ) {
        throw("No accession input") unless (scalar(@$acc_to_fetch));
        print "BatchFetch fields to get : " . join ( " " , @$fields )."\n" if $self->{verbose};
        my @accessions;
        my $index = 0;
        my $batch_size = 20;
        while (my $acc = pop(@$acc_to_fetch)) {
            my ($tmp_acc) = $acc =~ /^(\w+)/;
            $acc_hash{$tmp_acc} = $acc;
            push(@accessions, $tmp_acc);
            $index++;
            if ($batch_size == $index) {
                push(@queries, '?query='.join('+OR+', map {'accession:'.$_} @accessions).'&format=txt');
                @accessions = ();
                $index = 0;
            }
        }
        push(@queries, '?query='.join('+OR+', map {'accession:'.$_} @accessions).'&format=txt') if(@accessions);
    }
    else {
        print " try to fetch $acc_to_fetch\n" if ($self->{verbose});
        my ($tmp_acc) = $acc_to_fetch =~ /^(\w+)/;
        $acc_hash{$tmp_acc} = $acc_to_fetch;
        push(@queries, $tmp_acc.'.txt');
    }

    my %all_entries;
    my @entries_not_found;
    foreach my $query (@queries) {
        my $parser = Bio::EnsEMBL::IO::Parser::EMBL->open('wget -qq -O - "'.BASE_URL.$query.'" |');
        while ($parser->next) {
            my $accession;
            foreach my $acc (@{$parser->get_accessions}) {
                if (exists $acc_hash{$acc}) {
                    $accession = $acc_hash{$acc};
                    delete $acc_hash{$acc};
                }
            }
            foreach my $field (@$fields) {
                my ($key, $value) = $self->get_field($parser, $field) =~ /(\w{2})\s+(.*)/;
                $all_entries{$accession}{$key} = $value;
            }
        }
    }
    foreach my $entry (values %acc_hash) {
        push(@entries_not_found, $entry);
    }
    return [\%all_entries , \@entries_not_found ] ;
}

=head2 get_field

 Arg [1]    : Bio::EnsEMBL::IO::Parser::EMBL object
 Arg [2]    : String, the field you want to retrieve in the mfetch format
 Example    : $self->get_field($parser, $field);
 Description: We need to retrieve the field the same way Mfetch does. Here is the list of field
              that mfetch uses. You can verify this list with: mfetch -w line
            OC  Taxon
            AC  acc
            CC  cc
            CC  cct
            SQ  cks
            DT  crd
            DR  dbn
            DE  des
            ID  div
            DR  dr
            FT  ftd
            FT  ftk
            FT  ftl
            GN  gen
            ID  id
            KW  key
            DT  lau
            DT  lsu
            ID  mol
            OS  org
            AC  pac
            PE  pe
            DR  prd
            ID  sl
            DT  sv
            OC  tax
            OX  txi
 Returntype : String, ([two letter EMBL ID]  [free text])
 Exceptions : Throw if the field you're asking has no handler


=cut

sub get_field {
    my ($self, $parser, $field) = @_;

    if ($field eq 'acc') {
        return 'AC  '.join(' ', @{$parser->get_raw_accessions});
    }
    elsif ($field eq 'org') {
        return 'OS  '.$parser->get_species;
    }
    elsif ($field eq 'pe') {
        return 'PE  '.$parser->get_pelevel;
    }
    elsif ($field eq 'crd') {
        return 'DT  '.join(' ', @{$parser->get_date});
    }
    elsif ($field eq 'des') {
        return 'DE  '.$parser->get_description;
    }
    elsif ($field =~ /^tax/i) {
        return 'OC  '.join(' ', @{$parser->get_classification});
    }
    else {
        throw("You need to add a handler for your field: $field!\n");
    }
}

sub uniprot_file {
    my($self, $uniprot_file) = @_;

    if (defined $uniprot_file) {
        my @files;
        foreach my $file (@$uniprot_file) {
            if (-e $uniprot_file) {
                push(@files, $uniprot_file);
            }
        }
        if (scalar(@files)) {
            $self->{uniprot_file} = \@files;
        }
        else {
            throw("Cannot hammer the UniProt website or use -website");
        }
    }
    return $self->{uniprot_file};
}

sub use_website {
    my($self, $use_web) = @_;

    $self->{use_web} = $use_web if (defined $use_web);
    return $self->{use_web};
}

sub verbose {
    my($self, $verbosity) = @_;

    $self->{verbose} = $verbosity if (defined $verbosity);
    return $self->{verbose};
}

1;
