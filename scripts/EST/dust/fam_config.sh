
# ASSUMES $FAM_DIR IS SET

if [ -z $FAM_DIR ]
then
    echo FAM_DIR is not set
    exit 1
fi

# REQUIRED FILES AND SCRIPTS

FAM_JOB_IDS=${FAM_DIR}/job_ids.txt
FAM_CONFIG=${FAM_DIR}/config.txt

FAM_CHECK_NODE=${FAM_DIR}/check_node.sh
FAM_RUN_JOB=${FAM_DIR}/run_job.sh
FAM_IS_COMPLETE=${FAM_DIR}/is_complete.sh

FAM_PROCESS_JOB=${FAM_DIR}/process_job.sh

# GENERATED FILES

FAM_INIT_FILE=${FAM_DIR}/init

FAM_BSUB_STDERR=${FAM_DIR}/bsub_stderr.txt
FAM_BSUB_STDOUT=${FAM_DIR}/bsub_stdout.txt

FAM_INITIAL_TO_SUBMIT=${FAM_DIR}/initial_to_submit.sh
FAM_TO_SUBMIT=${FAM_DIR}/to_submit.sh

# -//-
