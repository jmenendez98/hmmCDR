#!/bin/bash

rm -r cdr_detect_jobstore
rm -r wdl_temp

mkdir -p cdr_detect_output
mkdir -p wdl_temp

toil-wdl-runner \
	--jobStore cdr_detect_jobstore \
	--stats \
	--batchSystem single_machine \
	--maxCores 16 \
	--runLocalJobsOnWorkers \
	--disableProgress \
	--workDir wdl_temp \
	--caching false \
	--container singularity \
	--clean=never \
	--stats \
	cdr_detection.wdl inputs.json \
	-o outputs \
	-m outputs.json \
	2>&1 | tee "toil_runner_log.txt"
