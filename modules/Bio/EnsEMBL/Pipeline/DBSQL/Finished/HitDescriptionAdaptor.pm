=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::Finished::HitDescriptionAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


### Bio::EnsEMBL::Pipeline::DBSQL::Finished::HitDescriptionAdaptor

package Bio::EnsEMBL::Pipeline::DBSQL::Finished::HitDescriptionAdaptor;

use strict;
use Bio::EnsEMBL::Pipeline::Finished::HitDescription;
use base 'Bio::EnsEMBL::DBSQL::BaseAdaptor';

sub fetch_HitDescriptions_into_hash {
    my( $self, $hash ) = @_;
    
    return unless %$hash;
    
    my $sql = qq{
        SELECT hit_name
          , hit_length
          , hit_description
          , hit_taxon
          , hit_db
        FROM hit_description
        WHERE hit_name IN (
        };
    $sql .= join(',', map "'$_'", keys %$hash);
    $sql .= qq{\n)};
    #warn $sql;
    
    my $sth = $self->prepare($sql);
    $sth->execute;

    my( $name, $length, $desc, $taxon_id, $db_name );
    $sth->bind_columns(\$name, \$length, \$desc, \$taxon_id, \$db_name);

    while ($sth->fetch) {
        $hash->{$name} = bless
            {
                _hit_length     => $length,
                _description    => $desc,
                _taxon_id       => $taxon_id,
                _db_name        => $db_name,
            }, 'Bio::EnsEMBL::Pipeline::Finished::HitDescription';
    }
}

1;

__END__

=head1 NAME - Bio::EnsEMBL::Pipeline::DBSQL::Finished::HitDescriptionAdaptor

=head2 AUTHOR
James Gilbert B<email> jgrg@sanger.ac.uk    - original implementation

=head2 AUTHOR
Mustapha Larbaoui B<email> ml6@sanger.ac.uk - new pipeline

