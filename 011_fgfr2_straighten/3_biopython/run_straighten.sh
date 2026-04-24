#!/bin/bash
#SBATCH --job-name=FGFR2_perfect
#SBATCH --output=fgfr2_perfect_%j.out
#SBATCH --error=fgfr2_perfect_%j.err
#SBATCH --partition=cpu
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

echo "=== FGFR2 Linker Straightener (PERFECT) ==="
echo "Start: $(date)"

ENV_PATH="/scratch/users/k25053715/MSc_research/0_project_env/biopython_env"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_PATH"

cat > straighten_linker.py << 'EOF'
#!/usr/bin/env python3
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
    print(f"Stretching linker {start}-{end}: {n_res} residues → {target_length:.1f}A")

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
    print(f"  Result: {new_length:.1f}A (bend={bend_strength})")

parser = PDBParser(QUIET=True)
s = parser.get_structure('FGFR2', './010_protein_preMD.pdb')
model = s[0]
chainA = model['A']

print(f"Total residues in chain A: {len(list(chainA.get_residues()))}")

linkers = [
    ('ECD_linker', (359, 377), 18.0, 0.2),
    ('ICD_linker', (399, 418), 25.0, 0.3)
]

for name, range_, target_len, bend in linkers:
    stretch_linker_selective(chainA, range_, target_len, bend)
    print(f"✓ {name} completed")

# Save
output = 'FGFR2_clean_linkers.pdb'
io = PDBIO()
io.set_structure(s)
io.save(output)

print(f"\nSAVED: {output} ({os.path.getsize(output)/1024:.0f}KB)")
print("\nChimeraX commands:")
print(f"open {output}")
print("select :359-377.A,:399-418.A; style cartoon")
print("color /A:359-377#1 yellow; color /A:399-418#1 cyan")
print("select clear; transparency 50")
print("measure distance :378.A@CA :398.A@CA")
print("lighting soft")
EOF

python straighten_linker.py
rm straighten_linker.py

if [ -f "FGFR2_straightened_linkers.pdb" ]; then
    ls -lh FGFR2_straightened_linkers.pdb
    echo "PERFECT! Linkers straightened. Ready for MD/ChimeraX."
    echo "Next: sbatch openmm_md.slurm"
else
    echo "Failed!"
    exit 1
fi

echo "Done: $(date)"
