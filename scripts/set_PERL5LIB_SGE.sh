
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute
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


# MAIN CONFIGURATION TO CHECKOUT-ROOT
##########################################

  # symlink to workdir

  setenv MYHOME '/home/ensembl' 
 
 
  setenv refdbhost  127.0.0.1 
  setenv db  'jhv_homo_sapiens_core_59_37d_cloud_ref' 
  setenv dbcon '-dbhost 127.0.0.1  -dbuser ens-training  -dbport 3306  -dbpass workshop  -dbname ' 
  setenv gb ref


  set prompt = "$db  %c2>"


  #
  # set PERL5LIB to ensembl-checkouts ( diff. checkouts below ) 
  # 

  setenv ANA_CVS  $MYHOME/cvs_checkout 
  
  setenv ENSEMBL           $ANA_CVS/ensembl_59 
  setenv ENSEMBL_ANALYSIS  $ANA_CVS/ensembl-analysis
  setenv ENSEMBL_PIPELINE  $ANA_CVS/ensembl-pipeline 
  setenv ORGANISM_CONFIG   $ANA_CVS/ensembl-config/aws_cloud/test-body-map


  
echo " "
echo "checking pathes...."
echo " " 



set path_to_check = (  ENSEMBL   ORGANISM_CONFIG  ENSEMBL_PIPELINE  ENSEMBL_ANALYSIS   ) 


set vars_to_check = ( $ENSEMBL  $ORGANISM_CONFIG  $ENSEMBL_PIPELINE $ENSEMBL_ANALYSIS  ) 








set i = 0 
foreach var ( $vars_to_check )
@ i++
  if (! -e $var ) then
     echo "HEY ! One of your pathes does not exist: "
     echo "ERROR: $path_to_check[$i] $var "
     exit
  else
     echo "checking "'$'"$path_to_check[$i]  $var  ======> OK "
  endif
end 



# set the pathes supplied 
setenv PERL5LIB ''


# ENSEMBL MODULES 

setenv PERL5LIB $ENSEMBL/modules
setenv PERL5LIB ${PERL5LIB}:$ENSEMBL_PIPELINE/modules
setenv PERL5LIB ${PERL5LIB}:$ENSEMBL_ANALYSIS/modules
setenv PERL5LIB ${PERL5LIB}:$ENSEMBL_PIPELINE/scripts

setenv OC $ORGANISM_CONFIG


# bioperl

setenv PERL5LIB ${PERL5LIB}:${ANA_CVS}/cvs_checkout/bioperl-live-12
setenv PERL5LIB ${PERL5LIB}:${ANA_CVS}/cvs_checkout/bioperl-run-1.5.1


# own modules

setenv PERL5LIB $ORGANISM_CONFIG/modules:${PERL5LIB}

setenv ES   $ENSEMBL/scripts
setenv PS   $ENSEMBL_PIPELINE/scripts
setenv AS   $ENSEMBL_ANALYSIS/scripts
setenv TS   $ENSEMBL_PIPELINE/misc-scripts/

setenv BQ   $OC/modules/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm 
setenv DB   $OC/modules/Bio/EnsEMBL/Analysis/Config/Databases.pm 

setenv ES   $ENSEMBL/scripts
setenv PS   $ENSEMBL_PIPELINE/scripts
setenv AS   $ENSEMBL_ANALYSIS/scripts
setenv TS   $ENSEMBL_PIPELINE/misc-scripts/


echo " " 
echo "All pathes are existing and ok, all previous vars are overwritten" 
echo " " 
echo " " 

set scripts = ( ES PS AS TS  )

echo " " 
echo " Shortcuts to Scripts: "
echo "-------------------------"
echo " " 

foreach SC ( $scripts )
  set variable_content  = (` env | grep -w  $SC= ` ) 
  echo '$'"${SC} ==>   $variable_content " 
end 

unset scripts  


echo " " 
echo " " 
echo " " 

set shortcuts = ( PI PE EC )

echo " " 
echo " Shortcuts to configs "
echo "-------------------------"
echo " " 

foreach SC ( $shortcuts )
  set variable_content  = (` env | grep -w  $SC= ` ) 
  echo '$'"${SC} ==>   $variable_content " 
end 

 
echo " " 
echo " " 
echo " " 
echo " " 
echo " " 
echo " setenv PATH /usr/local/ensembl/bin/" 


setenv PATH /usr/local/ensembl/bin/:${PATH} 

