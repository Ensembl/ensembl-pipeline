#!/bin/sh

if [ $# -ne 1 ]
then
   echo "Usage: run_job.sh <jobid>"
   exit 1
fi

JOBID=$1

CONFIG_FILE=`dirname $0`/config.txt
. $CONFIG_FILE

$DUSTEXE ${SEQFILE} 

