#!/bin/sh

if [ $# -ne 2 ]
then
   echo "Usage: check_node.sh <fam_dir> <jobid>"
   exit 1
fi

FAM_DIR=$1
JOBID=$2

#FAM_DIR=`dirname $0`
CONFIG_FILE=${FAM_DIR}/config.txt
. $CONFIG_FILE

if [ -f $DUSTEXE ]
then
    echo Found $DUSTEXE on `hostname`
else
    echo Cannot find $DUSTEXE on `hostname`
    exit 1
fi

if [ -f $SEQFILE ]
then
    echo Found $SEQFILE on `hostname`
else
    echo Cannot find $SEQFILE on `hostname`
    exit 1
fi

echo Host `hostname` is OK
exit 0

