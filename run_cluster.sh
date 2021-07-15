#!/bin/bash
set -eou pipefail

num_jobs=300
snakemake --cluster-config cluster.yaml --cluster-sync "srun --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --partition {cluster.partition} --output {cluster.output} --error {cluster.error} --job-name {cluster.name}" -j $num_jobs --keep-going

