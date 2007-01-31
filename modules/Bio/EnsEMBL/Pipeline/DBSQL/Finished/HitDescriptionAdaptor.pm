
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

