#!/bin/sh
#

if [ $# -ne 2 ]
then
    echo "usage: fam.sh <fam_dir> <mode>"
    echo
    echo "Modes available: [ init | submit | clean | collate | time ]"
    exit 1
fi

FAM_DIR=$1
FAM_MODE=$2
FAM_JOB_IDS=${FAM_DIR}/job_ids.txt
FAM_CONFIG=${FAM_DIR}/config.txt
FAM_RUN_JOB=${FAM_DIR}/run_job.sh
FAM_CHECK_NODE=${FAM_DIR}/check_node.sh
FAM_IS_COMPLETE=${FAM_DIR}/is_complete.sh
#FAM_PROCESS_JOB=`dirname $0`/process_job.sh
FAM_PROCESS_JOB=${FAM_DIR}/process_job.sh
FAM_BSUB_TO_SUBMIT=${FAM_DIR}/to_submit.sh

fam_check_file(){
    if [ -f $1 ]
    then
        return 0
    else
        echo "Missing file:" $1
        exit 1
    fi
    return
    }

fam_dir_is_valid(){
    fam_check_file $FAM_JOB_IDS
    fam_check_file $FAM_CONFIG
    fam_check_file $FAM_RUN_JOB
    fam_check_file $FAM_IS_COMPLETE
    fam_check_file $FAM_CHECK_NODE
#    echo Number of jobs: `wc -l ${FAM_JOB_IDS}`
    return
    }

fam_init(){
    if [ -f ${FAM_DIR}/init ]
    then
        echo Already initialised
    else
        echo "Initialising:" $FAM_DIR
        touch ${FAM_DIR}/init
        mkdir ${FAM_DIR}/result
        mkdir ${FAM_DIR}/result/stderr
        mkdir ${FAM_DIR}/result/stdout
        mkdir ${FAM_DIR}/bsub
        mkdir ${FAM_DIR}/bsub/stderr
        mkdir ${FAM_DIR}/bsub/stdout
        mkdir ${FAM_DIR}/completed
    fi
    return
    }

fam_clean(){
    echo Removing all results
    rm -f ${FAM_DIR}/init
    rm -fr ${FAM_DIR}/result
    rm -fr ${FAM_DIR}/bsub
    rm -fr ${FAM_DIR}/completed
    return
    }

fam_clean_job(){
    rm -f ${FAM_DIR}/bsub/stdout/$JOBID
    rm -f ${FAM_DIR}/bsub/stderr/$JOBID
    rm -f ${FAM_DIR}/result/stdout/${JOBID}.gz
    rm -f ${FAM_DIR}/result/stderr/$JOBID
    rm -f ${FAM_DIR}/completed/$JOBID
    return
    }

fam_submit(){
    rm $FAM_BSUB_TO_SUBMIT
    cat $FAM_JOB_IDS | while read JOBID
    do
        JOB_MARKER=${FAM_DIR}/completed/$JOBID
        if [ -f $JOB_MARKER ]
        then
            echo Already completed: $JOBID
        else
            $FAM_IS_COMPLETE $JOBID
            STATUS=$?
            if [ $STATUS -eq 0 ]
            then
                touch $JOB_MARKER
            else
                fam_clean_job $JOBID
                echo Submitting: $JOBID
                echo bsub $BSUB_OPTIONS \
                 -E \"${FAM_DIR}/check_node.sh ${FAM_DIR} ${JOBID} \" \
                 -o ${FAM_DIR}/bsub/stdout/$JOBID \
                 -e ${FAM_DIR}/bsub/stderr/$JOBID \
                 ${FAM_PROCESS_JOB} ${FAM_DIR} $JOBID \
                 >> $FAM_BSUB_TO_SUBMIT
            fi
        fi
    done
    echo Jobs to submit: `wc -l $FAM_BSUB_TO_SUBMIT`
    chmod +x $FAM_BSUB_TO_SUBMIT
    return
    }

fam_collate(){
    cat $FAM_JOB_IDS | while read JOBID
    do
        gzip -dc ${FAM_DIR}/result/stdout/${JOBID}.gz
    done
    return
    }

fam_time(){
    cat $FAM_JOB_IDS | while read JOBID
    do
        grep "CPU time" ${FAM_DIR}/bsub/stdout/$JOBID \
        | awk '{print $(NF-1)}'
    done
    return
    }

# -----------------
# START OF COMMANDS

fam_dir_is_valid
. $FAM_CONFIG

case $FAM_MODE in
    "init")
    fam_init $FAM_DIR
    ;;
    "submit")
    fam_submit $FAM_DIR
    ;;
    "clean")
    fam_clean $FAM_DIR
    ;;
    "collate")
    fam_collate $FAM_DIR
    ;;
    "time")
    fam_time $FAM_DIR
    ;;
    *)
    echo Unknown fam mode: $FAM_MODE
    ;;
esac


