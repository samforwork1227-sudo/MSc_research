#!/bin/bash
#SBATCH --job-name=fgfr2_x20_gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --time=02:00:00

echo "=== Load modules ==="
module purge
module load python/3.11.6-gcc-13.2.0 cuda/12.2.1-gcc-13.2.0

# ====== Path Configuration ======
WORKDIR="/scratch/users/k25053715/MSc_research/0_project_env"
RUN_DIR="/scratch/users/k25053715/MSc_research/1_boltz2_project/2_run/010_fgfr2_20lipids"
YAML_PATH="$RUN_DIR/fgfr2_20lipids.yaml"

source "$WORKDIR/boltz2_venv/bin/activate"

export PATH="/usr/local/cuda/12.2/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/cuda/12.2/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "=== GPU check (run time) ==="
nvidia-smi 2>/dev/null || echo "nvidia-smi: $(which nvidia-smi)"
python -c "
import torch
print('Torch:', torch.__version__)
print('CUDA:', torch.cuda.is_available())
print('CUDA version:', torch.version.cuda if torch.cuda.is_available() else 'N/A')
"

mkdir -p "$RUN_DIR"
cat > "$YAML_PATH" << EOF
version: 1
sequences:
  - protein:
      id: A
      msa: empty
      sequence: "MVSWGRFICLVVVTMATLSLARPSFSLVEDTTLEPEEPPTKYQISQPEVYVAAPGESLEVRCLLKDAAVISWTKDGVHLGPNNRTVLIGEYLQIKGATPRDSGLYACTASRTVDSETWYFMVNVTDAISSGDDEDDTDGAEDFVSENSNNKRAPYWTNTEKMEKRLHAVPAANTVKFRCPAGGNPMPTMRWLKNGKEFKQEHRIGGYKVRNQHWSLIMESVVPSDKGNYTCVVENEYGSINHTYHLDVVERSPHRPILQAGLPANASTVVGGDVEFVCKVYSDAQPHIQWIKHVEKNGSKYGPDGLPYLKVLKAAGVNTTDKEIEVLYIRNVTFEDAGEYTCLAGNSIGISFHSAWLTVLPAPGREKEITASPDYLEIAIYCIGVFLIACMVVTVILCRMKNTTKKPDFSSQPAVHKLTKRIPLRRQVTVSAESSSSMNSNTPLVRITTRLSSTADTPMLAGVSEYELPEDPKWEFPRDKLTLGKPLGEGCFGQVVMAEAVGIDKDKPKEAVTVAVKMLKDDATEKDLSDLVSEMEMMKMIGKHKNIINLLGACTQDGPLYVIVEYASKGNLREYLRARRPPGMEYSYDINRVPEEQMTFKDLVSCTYQLARGMEYLASQKCIHRDLAARNVLVTENNVMKIADFGLARDINNIDYYKKTTNGRLPVKWMAPEALFDRVYTHQSDVWSFGVLMWEIFTLGGSPYPGIPVEELFKLLKEGHRMDKPANCTNELYMMMRDCWHAVPSQRPTFKQLVEDLDRILTLTTNEEYLDLSQPLEQYSPSYPDTRSSCSSGDDSVFSPDPMPYEPCLPQYPHINGSVKT"
EOF

LIGAND_SMILES="CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC"
for i in {B..U}; do
  cat >> "$YAML_PATH" << EOF
  - ligand:
      id: $i
      smiles: "$LIGAND_SMILES"
EOF
done

cat >> "$YAML_PATH" << 'EOF'
constraints:
  - pocket:
      binder: B
      contacts: &id001
      - - A
        - 378
      - - A
        - 398
      max_distance: 6.0
EOF

for i in {C..U}; do
  cat >> "$YAML_PATH" << EOF
  - pocket:
      binder: $i
      contacts: *id001
      max_distance: 6.0
EOF
done

echo "Fixed YAML (B-U letters): $YAML_PATH"
head -30 "$YAML_PATH" | tail -10

echo "=== FGFR2 + Lipid analysis ==="
export TORCH_HOME="$WORKDIR/.cache/torch"
export HF_HOME="$WORKDIR/.cache/huggingface"

srun boltz predict "$RUN_DIR/fgfr2_20lipids.yaml" \
  --out_dir "$RUN_DIR" \
  --use_msa_server \
  --accelerator gpu --devices 1

echo "DONE! Check $RUN_DIR/predictions/*.cif"
