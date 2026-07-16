#!/bin/bash -euo pipefail
echo "Using clean PDB: FGFR1_clean_linkers.pdb"
J1=$(sbatch --parsable /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/1_protein_min.sh FGFR1_clean_linkers.pdb)
J2=$(sbatch --parsable --dependency=afterok:${J1} /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/2_protein_membrane.sh)
J3=$(sbatch --parsable --dependency=afterok:${J2} /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/3_system_min_coarse.sh)
J4=$(sbatch --parsable --dependency=afterok:${J3} /cephfs/volumes/hpc_data_usr/k25053715/3a800abf-2370-46d6-96d9-b8284b670110/MSc_research/nextflow/2_pipeline/015_fgfr1/4_system_min_fine.sh)
echo "OpenMM job chain for fgfr1_20lipids: ${J1} -> ${J2} -> ${J3} -> ${J4}" > openmm_jobs_submitted.txt
