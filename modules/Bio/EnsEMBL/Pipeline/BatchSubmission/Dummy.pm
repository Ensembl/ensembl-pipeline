

package Bio::EnsEMBL::Pipeline::BatchSubmission::Dummy;


use Bio::EnsEMBL::Pipeline::BatchSubmission;
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Pipeline::BatchSubmission);

#this module is currently just present for running the RuleManager locally
#could potentially implement some checking code but currently not necessary


sub temp_filename{
  return '';
}

1;
