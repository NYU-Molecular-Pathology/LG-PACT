# Makefile for running the NGS607 Nextflow DSL1 pipeline on UltraViolet
SHELL := /bin/bash

NXF_VER := 19.07.0
NXF_ANSI_LOG := false
export NXF_VER NXF_ANSI_LOG

TIMESTAMP ?= $(shell date +"%Y%m%d-%H%M%S")
DIRNAME ?= $(notdir $(CURDIR))

REMOTE_http := https://github.com/NYU-Molecular-Pathology/LG-PACT.git

NXF_PROFILE := ultraviolet
NXF_SCRIPT ?= main.nf
NF_NGS_ONLY := main-sans-sophia.nf

EXTRA_PARAMS ?=

LOGFILE := logs/log.$(TIMESTAMP).out

SUBMITTED := .nextflow.submitted
SUBMITTED_NGS_ONLY := .nextflow.ngs-only.submitted
SUBMITTED_BAF_ONLY := .baf-only.submitted

BAF_SCRIPT := /gpfs/data/molecpathlab/development/post_pact_qc_workflow/baf_nf_sbatch.sh

.PHONY: backup run submit submit-ngs-only submit-baf-only run-hpc install remove-framework clean fix-permissions update remote deploy check-runid check-fastqdir config config-add-fastqdirs samplesheet pairs

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

# ~~~~~ DEPLOY ANALYSIS DIRECTORY ~~~~~ #

FASTQDIR :=
FASTQDIRS :=
DEMUX_SAMPLESHEET :=
DEMUX_SAMPLESHEET_output := demux-samplesheet.csv
RUNID :=
PRODDIR := /gpfs/data/molecpathlab/production/NGS607

check-runid:
	@[ -z "$(RUNID)" ] && printf "invalid RUNID specified: $(RUNID)\n" && exit 1 || :

check-fastqdir:
	@[ -z "$(FASTQDIR)" ] && printf "invalid FASTQDIR specified: $(FASTQDIR)\n" && exit 1 || :
	@[ ! -d "$(FASTQDIR)" ] && printf "FASTQDIR is not a valid directory: $(FASTQDIR)\n" && exit 1 || :

deploy:
	@$(MAKE) check-runid
	@$(MAKE) check-fastqdir
	@repo_dir="$${PWD}" && \
	output_dir="$(PRODDIR)/$(RUNID)" && \
	echo ">>> Setting up new repo in location: $${output_dir}" && \
	git clone --recursive "$${repo_dir}" "$${output_dir}" && \
	cd "$${output_dir}" && \
	echo ">>> Linking input directory: $(FASTQDIR)" && \
	ln -s "$(FASTQDIR)" input && \
	if [ -e "$(DEMUX_SAMPLESHEET)" ]; then \
		/bin/cp "$(DEMUX_SAMPLESHEET)" "$(DEMUX_SAMPLESHEET_output)" ; \
	fi && \
	echo ">>> Creating input fastq sheet" && \
	python generate-samplesheets.py --name-mode "$(NAMEMODE)" "$(FASTQDIR)" && \
	echo ">>> Creating config file..." && \
	$(MAKE) config CONFIG_OUTPUT="$${output_dir}/config.json" && \
	if [ -e "$(DEMUX_SAMPLESHEET_output)" ]; then \
		echo ">>> Adding demux samplesheet to config" && \
		$(MAKE) config DEMUX_SAMPLESHEET="$(DEMUX_SAMPLESHEET_output)" CONFIG_OUTPUT="$${output_dir}/config.json" ; \
	fi && \
	if [ -e "$(DEMUX_SAMPLESHEET_output)" ]; then \
		$(MAKE) pairs PAIRS_SHEET="$(DEMUX_SAMPLESHEET_output)" PAIRS_MODE=demux ; \
	fi && \
	printf ">>> NGS607 analysis directory prepared: $${output_dir}\n"

CONFIG_INPUT := .config.json
CONFIG_OUTPUT := config.json

$(CONFIG_OUTPUT):
	@echo ">>> Creating $(CONFIG_OUTPUT)"
	@cp "$(CONFIG_INPUT)" "$(CONFIG_OUTPUT)"

SAMPLESHEET :=

config: $(CONFIG_OUTPUT)
	@[ -n "$(RUNID)" ] && echo ">>> Updating runID config" && python config.py --update "$(CONFIG_OUTPUT)" --runID "$(RUNID)" || :
	@[ -n "$(SAMPLESHEET)" ] && echo ">>> Updating samplesheet config" && python config.py --update "$(CONFIG_OUTPUT)" --samplesheet "$(SAMPLESHEET)" || :
	@[ -n "$(DEMUX_SAMPLESHEET)" ] && echo ">>> Updating demultiplexing samplesheet config" && python config.py --update "$(CONFIG_OUTPUT)" --demux-samplesheet "$(DEMUX_SAMPLESHEET)" || :
	@[ -n "$(FASTQDIR)" ] && echo ">>> Updating fastqDirs config" && python config.py --update "$(CONFIG_OUTPUT)" --fastqDirs "$(FASTQDIR)" || :
	@[ -n "$(FASTQDIRS)" ] && echo ">>> Adding fastq dirs to config" && python config.py --update "$(CONFIG_OUTPUT)" --fastqDirs $(FASTQDIRS) || :

config-add-fastqdirs:
	@if [ ! -z "$(FASTQDIRS)" ]; then \
		for fastqdir in $(FASTQDIRS); do \
			echo ">>> Adding $${fastqdir} to $(CONFIG_OUTPUT)" ; \
			$(MAKE) config FASTQDIR="$${fastqdir}" ; \
		done ; \
	else \
		echo ">>> ERROR: no FASTQDIRS passed" ; \
	fi

NAMEMODE := noLaneSplit
SAMPLESHEET_OUTPUT := samples.analysis.tsv

samplesheet:
	@echo ">>> Getting fastqdirs from config file: $(CONFIG_OUTPUT)" && \
	fastqdirs="$$(python -c 'import json; fastq_dirs = json.load(open("$(CONFIG_OUTPUT)")).get("fastqDirs", None); print(" ".join(fastq_dirs) if fastq_dirs else "")')" && \
	echo ">>> loaded fastqdirs: $${fastqdirs}" && \
	echo ">>> Generating samplesheet '$(SAMPLESHEET_OUTPUT)' for fastqdirs" && \
	python generate-samplesheets.py $(EP) --samples-analysis-tsv "$(SAMPLESHEET_OUTPUT)" --name-mode "$(NAMEMODE)" $${fastqdirs}

PAIRS_SHEET := samples.pairs.csv
PAIRS_MODE := sns

pairs:
	@if [ ! -e "$(SAMPLESHEET_OUTPUT)" ]; then $(MAKE) samplesheet; fi && \
	if [ ! -e "$(PAIRS_SHEET)" ]; then \
		echo ">>> ERROR: PAIRS_SHEET does not exist: $(PAIRS_SHEET)" ; \
		exit 1 ; \
	fi && \
	if [ "$(PAIRS_MODE)" == "sns" ]; then \
		echo ">>> Updating samplesheet with sample pairs from sns style sheet" ; \
		python update-samplesheets.py --tumor-normal-sheet "$(PAIRS_SHEET)" --pairs-tumor-colname '#SAMPLE-T' --pairs-normal-colname '#SAMPLE-N' ; \
	elif [ "$(PAIRS_MODE)" == "demux" ]; then \
		echo ">>> Updating samplesheet with sample pairs from demultiplexing style sheet" ; \
		python bin/demux2tumor_normal_sheet.py "$(PAIRS_SHEET)" samples.tumor.normal.csv && \
		python update-samplesheets.py --tumor-normal-sheet samples.tumor.normal.csv ; \
	else \
		echo ">>> ERROR: PAIRS_MODE not recognized: $(PAIRS_MODE)" ; \
		exit 1 ; \
	fi

remote:
	@echo ">>> Setting git remote origin to $(REMOTE_http)"
	@git remote set-url origin "$(REMOTE_http)"

update: remote
	@echo ">>> Updating repo"
	@git pull
	@echo ">>> Updating git submodules"
	@git submodule update --recursive --remote --init
	@if [ -f nextflow ]; then \
		echo ">>> Removing old Nextflow" && \
		rm -f nextflow && \
		echo ">>> Reinstalling Nextflow" && \
		$(MAKE) install ; \
	else \
		$(MAKE) install ; \
	fi

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
		--ntasks-per-node=1 \
		-c 8 \
		--mem=48G \
		--exclusive \
		--export=ALL \
		submit.nextflow.sbatch.sh "$(TIMESTAMP)" "$(LOGFILE)" "output/logs" "$(EXTRA_PARAMS)" "main.nf")" && \
	baf_job_id="$$(sbatch --parsable \
		--dependency=afterok:$${job_id} \
		--kill-on-invalid-dep=yes \
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
		--ntasks-per-node=1 \
		-c 8 \
		--mem=48G \
		--exclusive \
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