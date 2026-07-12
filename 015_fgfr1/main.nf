params.project_root   = '/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1'
params.outdir         = "${params.project_root}/results"

params.boltz_venv     = '/scratch/users/k25053715/MSc_research/0_project_env/boltz2_venv'
params.boltz_seed = 42
params.biopython_env  = '/scratch/users/k25053715/MSc_research/0_project_env/biopython_env'
params.openmm_env     = '/scratch/users/k25053715/MSc_research/0_project_env/openmm'
params.torch_cache    = '/scratch/users/k25053715/MSc_research/0_project_env/.cache/torch'
params.hf_cache       = '/scratch/users/k25053715/MSc_research/0_project_env/.cache/huggingface'
params.sequence       = 'MWSWKCLLFWAVLVTATLCTARPSPTLPEQAQPWGAPVEVESFLVHPGDLLQLRCRLRDDVQSINWLRDGVQLAESNRTRITGEEVEVQDSVPADSGLYACVTSSPSGSDTTYFSVNVSDALPSSEDDDDDDDSSSEEKETDNTKPNRMPVAPYWTSPEKMEKKLHAVPAAKTVKFKCPSSGTPNPTLRWLKNGKEFKPDHRIGGYKVRYATWSIIMDSVVPSDKGNYTCIVENEYGSINHTYQLDVVERSPHRPILQAGLPANKTVALGSNVEFMCKVYSDPQPHIQWLKHIEVNGSKIGPDNLPYVQILKTAGVNTTDKEMEVLHLRNVSFEDAGEYTCLAGNSIGLSHHSAWLTVLEALEERPAVMTSPLYLEIIIYCTGAFLISCMVGSVIVYKMKSGTKKSDFHSQMAVHKLAKSIPLRRQVTVSADSSASMNSGVLLVRPSRLSSSGTPMLAGVSEYELPEDPRWELPRDRLVLGKPLGEGCFGQVVLAEAIGLDKDKPNRVTKVAVKMLKSDATEKDLSDLISEMEMMKMIGKHKNIINLLGACTQDGPLYVIVEYASKGNLREYLQARRPPGLEYCYNPSHNPEEQLSSKDLVSCAYQVARGMEYLASKKCIHRDLAARNVLVTEDNVMKIADFGLARDIHHIDYYKKTTNGRLPVKWMAPEALFDRIYTHQSDVWSFGVLLWEIFTLGGSPYPGVPVEELFKLLKEGHRMDKPSNCTNELYMMMRDCWHAVPSQRPTFKQLVEDLDRIVALTSNQEYLDLSMPLDQYSPSFPDTRSSTCSSGEDSVFSHEPLPEEPCLPRHPAQLANGGLKRR'
params.ligand_smiles  = 'CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC'
params.sample_id      = 'fgfr1_20lipids'
params.openmm_input_name = 'FGFR1_clean_linkers.pdb'
params.fixed_boltz_cif = '/scratch/users/k25053715/MSc_research/1_boltz2_project/2_run/010_fgfr1_20lipids/boltz_results_fgfr1_20lipids/predictions/fgfr1_20lipids/fgfr1_20lipids_model_0.cif'
params.tm_start = 377
params.tm_end   = 397

process MAKE_BOLTZ_INPUT {
  tag "${params.sample_id}"
  publishDir "${params.outdir}/01_boltz_input", mode: 'copy'

  output:
  tuple val(params.sample_id), path("${params.sample_id}.yaml")

  script:
  def ligands = ['B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U']
  def seq = params.sequence
  def smiles = params.ligand_smiles
  def tmStart = params.tm_start
  def tmEnd   = params.tm_end

  def seq_block = """\
version: 1
sequences:
  - protein:
      id: A
      msa: empty
      sequence: "${seq}"
"""

  def ligand_blocks = ligands.collect { x -> """\
  - ligand:
      id: ${x}
      smiles: "${smiles}"
""" }.join("")

  def constraint_blocks = """\
constraints:
""" + ligands.collect { x -> """\
  - pocket:
      binder: ${x}
      contacts:
        - - A
          - ${tmStart}
        - - A
          - ${tmEnd}
      max_distance: 6.0
""" }.join("")

  """
  cat > ${params.sample_id}.yaml <<'EOF'
${seq_block}${ligand_blocks}${constraint_blocks}
EOF

  echo "===== GENERATED YAML ====="
  cat ${params.sample_id}.yaml
  """
}

process RUN_BOLTZ2 {
    tag "$sample_id"
    publishDir "${params.outdir}/02_boltz2", mode: 'copy'

    input:
    tuple val(sample_id), path(yaml_file)

    output:
    tuple val(sample_id), path("${sample_id}_model_0.cif"), emit: cif
    tuple val(sample_id), path('boltz_out'), emit: boltz_dir

    script:
    """
    set -euo pipefail

    start_ts=\$(date +%s)
    echo "[RUN_BOLTZ2] Started at \$(date)"
    echo "[RUN_BOLTZ2] Input YAML: ${yaml_file}"

    module purge
    module load python/3.11.6-gcc-13.2.0 cuda/12.2.1-gcc-13.2.0
    source ${params.boltz_venv}/bin/activate

    export PATH="/usr/local/cuda/12.2/bin:\$PATH"
    export LD_LIBRARY_PATH="/usr/local/cuda/12.2/lib64\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
    export TORCH_HOME="${params.torch_cache}"
    export HF_HOME="${params.hf_cache}"

    mkdir -p boltz_out
    boltz predict "${yaml_file}" \\
        --out_dir boltz_out \\
        --use_msa_server \\
        --accelerator gpu \\
        --devices 1 \\
        --seed ${params.boltz_seed}

    find boltz_out -maxdepth 5 -type f

    cif_file=\$(find boltz_out -name "*.cif" | head -n 1)
    if [ -z "\$cif_file" ]; then
        echo "ERROR: no CIF file produced by boltz" >&2
        exit 1
    fi

    cp "\$cif_file" "${sample_id}_model_0.cif"
    test -s "${sample_id}_model_0.cif"

    end_ts=\$(date +%s)
    elapsed=\$(( end_ts - start_ts ))
    printf "[RUN_BOLTZ2] Finished in %d min %d s\n" \$((elapsed/60)) \$((elapsed%60))
    """
}

process RUNCHIMERAX {
    tag "$sample_id"
    publishDir params.outdir + "/03_chimerax", mode: 'copy'

    input:
    tuple val(sample_id), path(ciffile)

    output:
    tuple val(sample_id), path("protein_preMD.pdb")

    script:
    """
    start_ts=\$(date +%s)
    echo "[RUNCHIMERAX] Started at \$(date)"

    export HOME=\$PWD
    export XDG_CONFIG_HOME=\$PWD/.config
    export XDG_DATA_HOME=\$PWD/.local/share
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$HOME" "\$XDG_CONFIG_HOME" "\$XDG_DATA_HOME" "\$XDG_CACHE_HOME"

    echo "PWD=\$PWD"
    echo "ciffile=${ciffile}"
    ls -l
    ls -l "${ciffile}"

    cat > run.cxc <<EOF
close session
open ${ciffile}
cartoon
color /A:1-376 blue
color /A:377-397 red
color /A:398-822 green

show #1/B-U
style #1/B-U stick
color #1/B-U gold

measure inertia #1/A:377-397@CA

shape cylinder fromPoint 8.6225,-11.792,-4.3179 axis -0.422,0.112,-0.900 height 40 radius 0.3 color red coordinateSystem #1
shape cylinder fromPoint 8.6225,-11.792,-4.3179 axis 0.422,-0.112,0.900 height 40 radius 0.3 color red coordinateSystem #1

# ICD
torsion #1:A:408@C #1:A:409@N #1:A:409@CA #1:A:409@C 150 move large
torsion #1:A:419@C #1:A:420@N #1:A:420@CA #1:A:420@C -8 move small
torsion #1:A:399@C #1:A:400@N #1:A:400@CA #1:A:400@C 10 move large

# ECD hinges
torsion #1:A:253@C #1:A:254@N #1:A:254@CA #1:A:254@C -10 move small
torsion #1:A:376@C #1:A:377@N #1:A:377@CA #1:A:377@C 30 move small
torsion #1:A:154@C #1:A:155@N #1:A:155@CA #1:A:155@C 15 move small
torsion #1:A:376@C #1:A:377@N #1:A:377@CA #1:A:377@C 10 move small
torsion #1:A:376@C #1:A:377@N #1:A:377@CA #1:A:377@C 55 move small
torsion #1:A:250@C #1:A:251@N #1:A:251@CA #1:A:251@C 12 move small
torsion #1:A:255@C #1:A:256@N #1:A:256@CA #1:A:256@C 25 move small
torsion #1:A:360@C #1:A:361@N #1:A:361@CA #1:A:361@C 20 move small
torsion #1:A:370@C #1:A:371@N #1:A:371@CA #1:A:371@C 30 move small

save "\$PWD/protein_preMD.pdb" format pdb models #1
exit
EOF

    echo "===== run.cxc ====="
    cat run.cxc

    chimerax --nogui --script "\$PWD/run.cxc"
    echo "ChimeraX exit code: \$?"

    ls -l
    test -s protein_preMD.pdb

    end_ts=\$(date +%s)
    elapsed=\$(( end_ts - start_ts ))
    printf "[RUNCHIMERAX] Finished in %d min %d s\\n" \
      \$((elapsed/60)) \$((elapsed%60))
    """
}

process STRAIGHTEN_LINKERS {
    tag "$sample_id"
    publishDir "${params.outdir}/04_biopython", mode: 'copy'

    input:
    tuple val(sample_id), path(pre_md_pdb)

    output:
    tuple val(sample_id), path('FGFR1_clean_linkers.pdb')

    script:
    """
    start_ts=\$(date +%s)
    echo "[STRAIGHTEN_LINKERS] Started at \$(date)"
    
    source "\$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ${params.biopython_env}


    cat > straighten_linker.py << 'PYEOF'
from Bio.PDB import PDBParser, PDBIO
import numpy as np
import os

def stretch_linker_selective(chain, linker_range, target_length=25.0, bend_strength=0.25):
    start, end = linker_range
    linker_res = []
    for r in chain.get_residues():
        if start <= r.id[1] <= end and 'CA' in r:
            linker_res.append(r)
    if len(linker_res) < 2:
        print(f"Warning: linker {start}-{end} only {len(linker_res)} CA atoms")
        return
    n_res = len(linker_res)
    print(f"Stretching linker {start}-{end}: {n_res} residues -> {target_length:.1f}A")
    first_ca = linker_res[0]['CA'].get_coord()
    last_ca = linker_res[-1]['CA'].get_coord()
    axis = last_ca - first_ca
    orig_length = np.linalg.norm(axis)
    if orig_length == 0:
        print(f"Zero length linker {start}-{end}")
        return
    axis_unit = axis / orig_length
    for i, res in enumerate(linker_res):
        frac = i / (n_res - 1)
        bend = bend_strength * (2*frac - 1)**3 * (1 - abs(2*frac - 1))
        target_frac = np.clip(frac + bend, 0, 1)
        new_ca_pos = first_ca + target_frac * target_length * axis_unit
        old_ca = res['CA'].get_coord()
        delta = new_ca_pos - old_ca
        for atom in res:
            atom.set_coord(atom.get_coord() + delta)
    new_first = linker_res[0]['CA'].get_coord()
    new_last = linker_res[-1]['CA'].get_coord()
    new_length = np.linalg.norm(new_last - new_first)
    print(f"Result: {new_length:.1f}A (bend={bend_strength})")

parser = PDBParser(QUIET=True)
s = parser.get_structure('FGFR1', 'protein_preMD.pdb')
model = s[0]
chainA = model['A']
linkers = [
    ('ECD_linker', (359, 377), 18.0, 0.2),
    ('ICD_linker', (399, 418), 25.0, 0.3)
]
for name, range_, target_len, bend in linkers:
    stretch_linker_selective(chainA, range_, target_len, bend)
    print(f"{name} completed")

output = 'FGFR1_clean_linkers.pdb'
io = PDBIO()
io.set_structure(s)
io.save(output)
print(f"Saved {output} ({os.path.getsize(output)/1024:.0f}KB)")
PYEOF

    python straighten_linker.py
    test -s FGFR1_clean_linkers.pdb

    end_ts=\$(date +%s)
    elapsed=\$(( end_ts - start_ts ))
    printf "[RUNCHIMERAX] Finished in %d min %d s\n" \
      \$((elapsed/60)) \$((elapsed%60))

    """
}

process SUBMIT_OPENMM_JOBS {

  tag "${sample_id}"

  input:
  tuple val(sample_id), path(clean_pdb)

  output:
  tuple val(sample_id), path('openmm_jobs_submitted.txt')

  script:
  """
  set -euo pipefail

  echo "Using clean PDB in Nextflow task dir: ${clean_pdb}"
  test -s "${clean_pdb}"

  CLEAN_PDB_PATH=\$(realpath "${clean_pdb}")
  echo "Resolved clean PDB path: \$CLEAN_PDB_PATH"

  J1=\$(sbatch --parsable ${projectDir}/1_protein_min.sh "/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/04_biopython/FGFR1_clean_linkers.pdb")
  J2=\$(sbatch --parsable --dependency=afterok:\${J1} ${projectDir}/2_protein_membrane.sh /scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/3_fgfr1_clean_minimized.pdb)
  J3=\$(sbatch --parsable --dependency=afterok:\${J2} ${projectDir}/3_system_min_coarse.sh \
/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/4_with_membrane.pdb \
/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/4_with_membrane_system.xml)
  J4=\$(sbatch --parsable --dependency=afterok:\${J3} ${projectDir}/4_system_min_fine.sh \
/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/5_fgfr1_membrane_min_coarse_full.pdb \
/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/015_fgfr1/results/05_openmm/4_with_membrane_system.xml)

  echo "OpenMM job chain for ${sample_id}: \${J1} -> \${J2} -> \${J3} -> \${J4}" > openmm_jobs_submitted.txt
  """
}

workflow {
  boltz_input = MAKE_BOLTZ_INPUT()
  boltz_ch = RUN_BOLTZ2(boltz_input)
  chim_ch = RUNCHIMERAX(boltz_ch.cif)
  clean_ch = STRAIGHTEN_LINKERS(chim_ch)
  submit_openmm = SUBMIT_OPENMM_JOBS(clean_ch)
}
