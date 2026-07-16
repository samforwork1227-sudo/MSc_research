#!/bin/bash
#SBATCH --job-name=fgfr2_nf
#SBATCH --output=fgfr2_nf_%j.out
#SBATCH --error=fgfr2_nf_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

export BASE_DIR="/scratch/users/k25053715/MSc_research/0_project_env"
export JAVA_HOME="${BASE_DIR}/nextflow_env/jdk"
export PATH="${BASE_DIR}/nextflow_env/bin:${JAVA_HOME}/bin:${PATH}"
export NXF_JAVA_HOME="${JAVA_HOME}"

cd /scratch/users/k25053715/MSc_research/nextflow/2_pipeline/013_nextflow_pipeline

nextflow run main.nf -c nextflow.config
