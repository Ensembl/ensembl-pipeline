#!/bin/bash
export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-2-3:$PWD/ensembl-test/modules:$PWD/ensembl/modules:$PWD/ensembl-external/modules:$PWD/ensembl-analysis/scripts:$PWD/ensembl-analysis/modules:$PWD/ensembl-compara/modules:$PWD/modules:$PWD/ensembl-killlist/modules:$PWD/bioperl-run/lib

echo "Running test suite"
echo "Using $PERL5LIB"
rt=0
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t $SKIP_TESTS
else
# Preparing configs
  mkdir modules/Bio/EnsEMBL/Pipeline/Config/GeneBuild
  echo "1;" > modules/Bio/EnsEMBL/Pipeline/Config/GeneBuild/Databases.pm
  echo "1;" > modules/Bio/EnsEMBL/Pipeline/Tools/TranscriptUtils.pm
  # just test the basic syntax for all the scripts and modules - rename .example files first
  find $PWD/ensembl-analysis/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  find $PWD/modules -type f -name '*.example' | while read f; do mv "$f" "${f%.example}"; done
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  printf "\e[31mWe will not test:\e[0m\n - \e[33m%s\e[0m\n" "Annacode scripts"
  RES="! -path \"*Finished*\""
  P=( "Hmms" \
    "Pseudogenes" \
    "EST")
  for D in `seq 0 $((${#P[@]}-1))`; do
      printf " - \e[33m%s\n\e[0m" "${P[$D]}"
      RES=${RES}" ! -path \"*`basename ${P[$D]}`*\""
  done
  M=( "scripts/post_GeneBuild/post_GeneBuild_checks_denormalised.pl" )
  for S in `seq 0 $((${#M[@]}-1))`; do
      printf " - \e[33m%s\n\e[0m" "${M[$S]}"
      RES=${RES}" ! -name `basename ${M[$S]}`"
  done
  printf "  RES = \e[33m%s\e[0m\n" "$RES"
  find $PWD/scripts -type f -name "*.pl" `echo "$RES"` | xargs -i perl -c {}
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
  # Directory we don't want to check
  printf "\e[31mWe will not test:\e[0m\n - \e[33m%s\e[0m\n" "Annacode modules"
  RES="! -path \"*Finished*\""
  P=( "GeneDuplication" )
  for D in `seq 0 $((${#P[@]}-1))`; do
      printf " - \e[33m%s\n\e[0m" "${P[$D]}"
      RES=${RES}" ! -path \"*`basename ${P[$D]}`*\""
  done
# As long as EMBL parser has not been merge in ensembl-io master, we will avoid modules/Bio/EnsEMBL/Pipeline/SeqFetcher/UniProtKB.pm \ "scripts/post_GeneBuild/post_GeneBuild_checks_denormalised.pl"
  M=( "Bio/EnsEMBL/Pipeline/SeqFetcher/UniProtKB.pm" )
  for S in `seq 0 $((${#M[@]}-1))`; do
      printf " - \e[33m%s\n\e[0m" "${M[$S]}"
      RES=${RES}" ! -name `basename ${M[$S]}`"
  done
  printf "  RES = \e[33m%s\e[0m\n" "$RES"
  find $PWD/modules -type f -name "*.pm" `echo "$RES"` | xargs -i perl -c {}
  EXIT_CODE=$?
  if [ "$EXIT_CODE" -ne 0 ]; then
      rt=$EXIT_CODE
  fi
fi
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
