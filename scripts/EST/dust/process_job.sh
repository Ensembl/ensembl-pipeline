#!/bin/sh

if [ $# -ne 2 ]
then
    echo "Usage: process_job.sh <fam_dir> <job_id>"
    exit 1
fi

FAM_DIR=$1
JOB_ID=$2

FAM_CONFIG_FILE=${FAM_DIR}/fam_config.sh
. $FAM_CONFIG_FILE

TMP_WORK_DIR=/tmp/fam_${JOB_ID}_$$_`hostname`

TMP_STDOUT_FILE=${TMP_WORK_DIR}.stdout.gz
TMP_STDERR_FILE=${TMP_WORK_DIR}.stderr

${FAM_DIR}/run_job.sh ${JOB_ID} 2> $TMP_STDERR_FILE \
                       | gzip -9 > $TMP_STDOUT_FILE

mv $TMP_STDOUT_FILE ${FAM_DIR}/result/stdout/${JOB_ID}.gz
mv $TMP_STDERR_FILE ${FAM_DIR}/result/stderr/$JOB_ID
