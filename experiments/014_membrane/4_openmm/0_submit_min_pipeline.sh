#!/bin/bash
# submit_min_pipeline.sh

# 1. protein minimisation
#J1=$(sbatch --parsable 1_protein_min.sh)

# 2. membrane construction
#J2=$(sbatch --parsable --dependency=afterok:${J1} 2_protein_membrane.sh)
J2=$(sbatch --parsable 2_protein_membrane.sh)
# 3. system coarse minimisation
J3=$(sbatch --parsable --dependency=afterok:${J2} 3_system_min_coarse.sh)

# 4. system fine minimisation
J4=$(sbatch --parsable --dependency=afterok:${J3} 4_system_min_fine.sh)

echo "Submitted:"
#echo "  1_protein_min.sh        : $J1"
echo "  2_protein_membrane.sh   : $J2 (afterok:$J1)"
echo "  3_system_min_coarse.sh  : $J3 (afterok:$J2)"
echo "  4_system_min_fine.sh    : $J4 (afterok:$J3)"
