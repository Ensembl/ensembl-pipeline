# Script for operating the analysis pipeline
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 05.09.2000
# Last modified : 05.09.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this code under the same terms as perl itself




# outline of the algorithm

use strict;

use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::;


__END__

fetch_all_rules ..

fetch_object_attributes



for hotJobs (jobs which have been made and are running)
  check successful state,
  add the analysis to the inputIds attributes
  add the inputId to the hotIds list

for hotIds
  for each of the ids get the full attribute covering
  check all rule which could apply
  create all jobs and add them to hotJob

for all failed hot jobs
  reissue if appropriate
  mail if mailcount is not exceeded
  die if mailcount is exceeded

