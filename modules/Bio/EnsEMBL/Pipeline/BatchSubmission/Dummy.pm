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

=head1 NAME

Bio::EnsEMBL::Pipeline::BatchSubmission::Dummy - 

=head1 SYNOPSIS


=head1 DESCRIPTION

this module is currently just present for running the RuleManager locally
could potentially implement some checking code but currently not necessary

=head1 METHODS

=cut


# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/modules/Bio/EnsEMBL/Pipeline/BatchSubmission/Dummy.pm,v $
# $Version: $
package Bio::EnsEMBL::Pipeline::BatchSubmission::Dummy;


use warnings ;
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
