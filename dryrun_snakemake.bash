#!/bin/bash
set -xeuo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: bash $0 PATH_TO_CONFIG"
    exit 1
fi

CONFIG_FP=$1

snakemake \
    --jobs 300 \
    --configfile ${CONFIG_FP} \
    --cluster-config cluster.json \
    --latency-wait 180 \
    --notemp \
    --printshellcmds \
    --cluster \
    "qsub -cwd -r n -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -pe smp {threads}" \
    --dryrun
