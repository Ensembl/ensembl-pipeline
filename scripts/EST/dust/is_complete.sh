#!/bin/sh

if [ $# -ne 1 ]
then
   echo "Usage: is_complete.sh <jobid>"
   exit 1
fi

JOBID=$1

FAM_DIR=`dirname $0`
CONFIG_FILE=${FAM_DIR}/config.txt
. $CONFIG_FILE

OUTPUT_FILE="${FAM_DIR}/result/stdout/${JOBID}.gz"

if [ -f $OUTPUT_FILE ]
then
    echo Output present : $OUTPUT_FILE
else
    echo Output absent : $OUTPUT_FILE
    exit 1
fi

INPUTSEQS=`grep -c '^>' $SEQFILE`
OUTPUTSEQS=`zgrep -c '^>' $OUTPUT_FILE`

if [ $INPUTSEQS -eq $OUTPUTSEQS ]
then
    echo Job successful: $JOBID
    exit 0
else
    echo Job incomplete: $OUTPUT_FILE
    echo input seqs $INPUTSEQS output seqs $OUTPUTSEQS
    exit 1
fi

exit 1

