# Makefile for running the NGS607 Nextflow DSL1 pipeline on UltraViolet
SHELL := /bin/bash

NXF_VER := 19.07.0
NXF_ANSI_LOG := false
export NXF_VER NXF_ANSI_LOG

TIMESTAMP ?= $(shell date +"%Y%m%d-%H%M%S")
DIRNAME ?= $(notdir $(CURDIR))

NXF_PROFILE := ultraviolet
NXF_SCRIPT ?= main.nf
NF_NGS_ONLY := main-sans-sophia.nf

EXTRA_PARAMS ?=

LOGFILE := logs/log.$(TIMESTAMP).out

SUBMITTED := .nextflow.submitted
SUBMITTED_NGS_ONLY := .nextflow.ngs-only.submitted
SUBMITTED_BAF_ONLY := .baf-only.submitted

BAF_SCRIPT := /gpfs/data/molecpathlab/development/post_pact_qc_workflow/baf_nf_sbatch.sh

.PHONY: backup run submit submit-ngs-only submit-baf-only run-hpc install remove-framework clean fix-permissions

backup:
	@if [ -d "output" ]; then \
		echo ">>> Backing up output to output-backup/$(TIMESTAMP)" ; \
		mkdir -p "output-backup/$(TIMESTAMP)" ; \
		mv "output" "output-backup/$(TIMESTAMP)" ; \
	fi

run: backup submit

remove-framework:
	@if [ -e "$(HOME)/.nextflow/framework/$(NXF_VER)" ]; then \
		moved="$(HOME)/.nextflow/framework/$(NXF_VER).$$(date +%s)" ; \
		echo ">>> Moving stale framework $(HOME)/.nextflow/framework/$(NXF_VER) -> $${moved}" ; \
		mv "$(HOME)/.nextflow/framework/$(NXF_VER)" "$${moved}" ; \
	fi

./nextflow:
	@[ -d "$(HOME)/.nextflow/framework/$(NXF_VER)" ] && $(MAKE) remove-framework || :
	curl -fsSL get.nextflow.io | bash
	@echo ">>> Warming Nextflow $(NXF_VER) framework cache" ; \
	./nextflow -version >/dev/null 2>&1 || true

install: ./nextflow

submit: install
	@mkdir -p "logs"
	@if ! ( set -o noclobber ; printf 'PENDING\t%s\n' "$(TIMESTAMP)" > "$(SUBMITTED)" ) 2>/dev/null; then \
		echo ">>> Refusing to submit: $(SUBMITTED) already exists." ; \
		echo ">>> remove $(SUBMITTED) before submitting again" ; \
		exit 1 ; \
	fi
	@job_id="$$(sbatch --parsable \
		-D "$(CURDIR)" \
		-o "$(CURDIR)/logs/slurm-%j.log.$(TIMESTAMP).out" \
		-J "NGS607-$(DIRNAME)" \
		-p intellispace \
		--time=5-00:00:00 \
		--nodes=1 \
		--ntasks=1 \
		-c 8 \
		--mem=16G \
		--export=ALL \
		submit.nextflow.sbatch.sh "$(TIMESTAMP)" "$(LOGFILE)" "output/logs" "$(EXTRA_PARAMS)" "main.nf")" && \
	baf_job_id="$$(sbatch --parsable \
		--dependency=afterok:$${job_id} \
		"$(BAF_SCRIPT)" "$(DIRNAME)")" && \
	printf '%s\t%s\t%s\n' "$${job_id}" "$(TIMESTAMP)" "$${baf_job_id}" > "$(SUBMITTED)" && \
	printf '>>> Submitted driver job:            %s\n' "$${job_id}" && \
	printf '>>> Submitted BAF/CNV dependent job: %s\n' "$${baf_job_id}" && \
	printf '>>> BAF/CNV dependency:              afterok:%s\n' "$${job_id}" && \
	printf '>>> Slurm log:                       %s\n' "$(CURDIR)/logs/slurm-%j.log.$(TIMESTAMP).out" && \
	printf '>>> Nextflow stdout log:             %s\n' "$(LOGFILE)" || \
	{ status="$$?" ; rm -f "$(SUBMITTED)" ; exit "$${status}" ; }

submit-ngs-only: install
	@mkdir -p "logs"
	@if ! ( set -o noclobber ; printf 'PENDING\t%s\n' "$(TIMESTAMP)" > "$(SUBMITTED_NGS_ONLY)" ) 2>/dev/null; then \
		echo ">>> Refusing to submit: $(SUBMITTED_NGS_ONLY) already exists." ; \
		echo ">>> remove $(SUBMITTED_NGS_ONLY) before submitting again" ; \
		exit 1 ; \
	fi
	@job_id="$$(sbatch --parsable \
		-D "$(CURDIR)" \
		-o "$(CURDIR)/logs/slurm-%j.log.$(TIMESTAMP).out" \
		-J "NGS607-$(DIRNAME)-ngs-only" \
		-p intellispace \
		--time=5-00:00:00 \
		--nodes=1 \
		--ntasks=1 \
		-c 8 \
		--mem=16G \
		--export=ALL \
		submit.nextflow.sbatch.sh "$(TIMESTAMP)" "$(LOGFILE)" "output/logs" "$(EXTRA_PARAMS)" "$(NF_NGS_ONLY)")" && \
	printf '%s\t%s\t%s\n' "$${job_id}" "$(TIMESTAMP)" "$(NF_NGS_ONLY)" > "$(SUBMITTED_NGS_ONLY)" && \
	printf '>>> Submitted NGS-only driver job: %s\n' "$${job_id}" && \
	printf '>>> Nextflow script:              %s\n' "$(NF_NGS_ONLY)" && \
	printf '>>> Nextflow stdout log:          %s\n' "$(LOGFILE)" || \
	{ status="$$?" ; rm -f "$(SUBMITTED_NGS_ONLY)" ; exit "$${status}" ; }

submit-baf-only:
	@mkdir -p "logs"
	@if ! ( set -o noclobber ; printf 'PENDING\t%s\n' "$(TIMESTAMP)" > "$(SUBMITTED_BAF_ONLY)" ) 2>/dev/null; then \
		echo ">>> Refusing to submit: $(SUBMITTED_BAF_ONLY) already exists." ; \
		echo ">>> remove $(SUBMITTED_BAF_ONLY) before submitting again" ; \
		exit 1 ; \
	fi
	@baf_job_id="$$(sbatch --parsable "$(BAF_SCRIPT)" "$(DIRNAME)")" && \
	printf '%s\t%s\n' "$${baf_job_id}" "$(TIMESTAMP)" > "$(SUBMITTED_BAF_ONLY)" && \
	printf '>>> Submitted BAF/CNV-only job: %s\n' "$${baf_job_id}" && \
	printf '>>> BAF script:                 %s\n' "$(BAF_SCRIPT)" && \
	printf '>>> Run directory name:         %s\n' "$(DIRNAME)" || \
	{ status="$$?" ; rm -f "$(SUBMITTED_BAF_ONLY)" ; exit "$${status}" ; }

fix-permissions:
	@find . -type f -executable -exec chmod ug+X {} \;
	@find . -type d -exec chmod ug+rwxs {} \;
	@find . -type f -exec chmod ug+rw {} \;
	@find . ! -group "molecpathlab" -exec chgrp "molecpathlab" {} \;

run-hpc: install
	@mkdir -p "logs"
	@set -o pipefail ; \
	./nextflow run "$(NXF_SCRIPT)" \
		-profile "$(NXF_PROFILE)" \
		-resume \
		-with-dag flowchart.dot \
		$(EXTRA_PARAMS) 2>&1 | tee -a "$(LOGFILE)" ; \
	status="$${PIPESTATUS[0]}" ; \
	if [ "$${status}" -ne 0 ]; then exit "$${status}" ; fi ; \
	$(MAKE) fix-permissions

clean:
	@echo ">>> Removing run intermediates from $(CURDIR)"
	rm -rf \
		.nextflow.log .nextflow.log.* \
		flowchart.dot flowchart.dot.* flowchart.html flowchart*.dot flowchart*.png \
		nextflow.html nextflow.html.* report.html report.html.* \
		timeline.html timeline.html.* \
		trace.txt trace.txt.* trace*.txt trace*.txt.* \
		pipeline_info work.rm.txt .trace.hash.txt .nextflow \
		logs \
		.nextflow.jobid .nextflow.ngs-only.jobid .baf.jobid .baf-only.jobid \
		.nextflow.pid "$(SUBMITTED)" "$(SUBMITTED_NGS_ONLY)" "$(SUBMITTED_BAF_ONLY)"