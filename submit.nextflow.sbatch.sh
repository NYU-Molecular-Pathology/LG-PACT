#!/usr/bin/env bash
set -euo pipefail

timestamp="${1:?missing timestamp}"
run_log="${2:?missing Nextflow stdout log}"
publish_logs="${3:-output/logs}"
extra_params="${4:-}"
nxf_script="${5:-main.nf}"
nxf_pid_file=".nextflow.pid"

terminate_nextflow() {
    status="${1:-143}"
    trap - TERM INT HUP EXIT

    if [[ -f "${nxf_pid_file}" ]]; then
        nxf_pid="$(head -n 1 "${nxf_pid_file}" 2>/dev/null || true)"

        if [[ "${nxf_pid}" =~ ^[0-9]+$ ]] && \
        kill -0 "${nxf_pid}" 2>/dev/null; then
            echo ">>> Forwarding term signal to Nextflow PID ${nxf_pid}"
            kill "${nxf_pid}" 2>/dev/null || true

            for _ in $(seq 1 60); do
                kill -0 "${nxf_pid}" 2>/dev/null || break
                sleep 5
            done

            kill -0 "${nxf_pid}" 2>/dev/null && \
            kill -9 "${nxf_pid}" 2>/dev/null || true
        fi
    fi

    exit "${status}"
}

cleanup_on_exit() {
    status="$?"
    [[ "${status}" -ne 0 ]] && terminate_nextflow "${status}"
    exit 0
}

trap 'terminate_nextflow 143' TERM
trap 'terminate_nextflow 130' INT
trap 'terminate_nextflow 129' HUP
trap cleanup_on_exit EXIT

mkdir -p "${publish_logs}"

make run-hpc \
    TIMESTAMP="${timestamp}" \
    LOGFILE="${run_log}" \
    EXTRA_PARAMS="${extra_params}" \
    NXF_SCRIPT="${nxf_script}"