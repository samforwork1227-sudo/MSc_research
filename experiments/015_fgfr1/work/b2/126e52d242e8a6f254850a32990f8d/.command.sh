#!/bin/bash -euo pipefail
set -euo pipefail

echo "Using clean PDB in Nextflow task dir: FGFR1_clean_linkers.pdb"
test -s "FGFR1_clean_linkers.pdb"

CLEAN_PDB_PATH=$(realpath "FGFR1_clean_linkers.pdb")
echo "Resolved clean PDB path: $CLEAN_PDB_PATH"

J1=$(sbatch --parsable /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/1_protein_min.sh "/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/04_biopython/FGFR1_clean_linkers.pdb")
J2=$(sbatch --parsable --dependency=afterok:${J1} /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/2_protein_membrane.sh /scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/3_fgfr1_clean_minimized.pdb)
J3=$(sbatch --parsable --dependency=afterok:${J2} /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/3_system_min_coarse.sh /scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/4_with_membrane.pdb /scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/4_with_membrane_system.xml)
J4=$(sbatch --parsable --dependency=afterok:${J3} /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/4_system_min_fine.sh /scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/5_fgfr1_membrane_min_coarse_full.pdb /scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/4_with_membrane_system.xml)

echo "OpenMM job chain for fgfr1_20lipids: ${J1} -> ${J2} -> ${J3} -> ${J4}" > openmm_jobs_submitted.txt
