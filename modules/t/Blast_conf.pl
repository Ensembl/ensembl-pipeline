# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# Author: Simon Potter
# Creation: 11.03.2002

# configuration information for Blast header parsing
#
# each blast database needs a regular expression to
# enable the correct accession to be parsed from the
# header
# the code will throw if a (legal) regex is not
# defined for a database when it is used
#
# Example RE's are below. If you don't understand
# regexes, e-mail http://lists.ensembl.org/mailman/listinfo/dev for help!


BEGIN {
package main;

# Regex's for Blast fasta headers

%fasta_header_re = (
		    'embl_vertrna'   => '\w+\s+(\w+)',
		    'embl_vertrna-1'   => '\w+\s+(\w+)',
		    'mini_mrna.fa'   => '\w+\s+(\w+)',
		    'swall'          => '^\w+\s+(\w+)',
		    'swall-1'          => '^\w+\s+(\w+)',
		    'mini_protein.fa'          => '^\w+\s+(\w+)',
		    'dbHUMAN-EST'    => '^\w+\|\w+\|\w+\|([\w\.]+)',
		    'dbMOUSE-EST'    => '^\w+\|\w+\|\w+\|([\w\.]+)',
		    'dbOTHERS-EST'   => '^\w+\|\w+\|\w+\|([\w\.]+)',
		    'unigene.seq'    => '\/ug=([\w\.]+)\s/',
		    'C.elegans_nematode_mRNAs'  => '^(\w+)\s+',
		    'AI053588.fa' =>	'(\w+)',
		    'AP000074.pep' => '(\w+)\|',

);

}

1;


__END__


^\w+\s+(\w+)'
>143B_HUMAN P31946 CAA40621.1 CAA15497.1 AAH01359.1 Desc: 14-3-3
PROTEIN BETA/ALPHA (PROTEIN KINASE C INHIBITOR PROTEIN-1) (KCIP-1)
(PROTEIN 1054).

'\/ug=(\w+)'
>gnl|UG|At#S146917 Arabidopsis thaliana AT3g12580/T2E22_110 mRNA,
complete cds /cds=(102,2054) /gb=AY054183 /gi=15809831 /ug=At.1
/len=2284

>AP000074|GENSCAN_predicted_peptide_1|956_aa
