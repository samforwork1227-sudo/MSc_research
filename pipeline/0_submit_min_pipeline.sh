#!/bin/bash
set -euo pipefail

BASE_DIR="/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1"
OUTDIR="${BASE_DIR}/results/05_openmm"

INPUT_PDB="${BASE_DIR}/results/04_biopython/FGFR1_clean_linkers.pdb"

STEP1_OUT="${OUTDIR}/3_fgfr1_clean_minimized.pdb"
STEP2_PDB="${OUTDIR}/4_with_membrane.pdb"
STEP2_XML="${OUTDIR}/4_with_membrane_system.xml"
STEP3_OUT="${OUTDIR}/5_fgfr1_membrane_min_coarse_full.pdb"

mkdir -p "${OUTDIR}" "${OUTDIR}/python"

if [ ! -f "${INPUT_PDB}" ]; then
    echo "ERROR: input PDB not found: ${INPUT_PDB}"
    exit 1
fi

echo "Submitting OpenMM minimization pipeline..."
echo "Input PDB: ${INPUT_PDB}"
echo "Output dir: ${OUTDIR}"

# 1. protein minimisation
J1=$(sbatch --parsable "${BASE_DIR}/1_protein_min.sh" "${INPUT_PDB}")

# 2. membrane construction
J2=$(sbatch --parsable --dependency=afterok:${J1} \
    "${BASE_DIR}/2_protein_membrane.sh" \
    "${STEP1_OUT}")

# 3. system coarse minimisation
J3=$(sbatch --parsable --dependency=afterok:${J2} \
    "${BASE_DIR}/3_system_min_coarse.sh" \
    "${STEP2_PDB}" \
    "${STEP2_XML}")

# 4. system fine minimisation
J4=$(sbatch --parsable --dependency=afterok:${J3} \
    "${BASE_DIR}/4_system_min_fine.sh" \
    "${STEP3_OUT}" \
    "${STEP2_XML}")

echo "Submitted:"
echo " 1_protein_min.sh        : ${J1}"
echo " 2_protein_membrane.sh   : ${J2} (afterok:${J1})"
echo " 3_system_min_coarse.sh  : ${J3} (afterok:${J2})"
echo " 4_system_min_fine.sh    : ${J4} (afterok:${J3})"

echo
echo "Expected outputs:"
echo " step1 -> ${STEP1_OUT}"
echo " step2 -> ${STEP2_PDB}"
echo " step2 -> ${STEP2_XML}"
echo " step3 -> ${STEP3_OUT}"
echo " step4 -> ${OUTDIR}/6_fgfr1_membrane_minimized_full.pdb"
echo " step4 -> ${OUTDIR}/6_fgfr1_membrane_minimized_no_water.pdb"
