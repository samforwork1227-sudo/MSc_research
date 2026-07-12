#!/bin/bash -euo pipefail
start_ts=$(date +%s)
    echo "[RUN_BOLTZ2] Started at $(date)"

    module purge
    module load python/3.11.6-gcc-13.2.0 cuda/12.2.1-gcc-13.2.0
    source /scratch/users/k25053715/MSc_research/0_project_env/boltz2_venv/bin/activate

    export PATH="/usr/local/cuda/12.2/bin:$PATH"
    export LD_LIBRARY_PATH="/usr/local/cuda/12.2/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    export TORCH_HOME="/scratch/users/k25053715/MSc_research/0_project_env/.cache/torch"
    export HF_HOME="/scratch/users/k25053715/MSc_research/0_project_env/.cache/huggingface"

    mkdir -p boltz_out
    boltz predict fgfr1_20lipids.yaml       --out_dir boltz_out       --use_msa_server       --accelerator gpu       --devices 1       --seed 42

    find boltz_out -maxdepth 3 -type f

    cif_file=$(find boltz_out -name "*.cif" | head -n 1)
    if [ -z "$cif_file" ]; then
      echo "ERROR: no CIF file produced by boltz" >&2
      exit 1
    fi

    cp "$cif_file" fgfr1_20lipids_model_0.cif
    test -s fgfr1_20lipids_model_0.cif
    
    end_ts=$(date +%s)
    elapsed=$(( end_ts - start_ts ))
    printf "[RUN_BOLTZ2] Finished in %d min %d s
"       $((elapsed/60)) $((elapsed%60))
