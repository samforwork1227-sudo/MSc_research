params.project_root   = '/scratch/users/k25053715/MSc_research/nextflow/2_pipeline/013_nextflow_pipeline'
params.outdir         = "${params.project_root}/results"

params.boltz_venv     = '/scratch/users/k25053715/MSc_research/0_project_env/boltz2_venv'
params.biopython_env  = '/scratch/users/k25053715/MSc_research/0_project_env/biopython_env'
params.openmm_env     = '/scratch/users/k25053715/MSc_research/0_project_env/openmm'
params.torch_cache    = '/scratch/users/k25053715/MSc_research/0_project_env/.cache/torch'
params.hf_cache       = '/scratch/users/k25053715/MSc_research/0_project_env/.cache/huggingface'
params.sequence       = 'MVSWGRFICLVVVTMATLSLARPSFSLVEDTTLEPEEPPTKYQISQPEVYVAAPGESLEVRCLLKDAAVISWTKDGVHLGPNNRTVLIGEYLQIKGATPRDSGLYACTASRTVDSETWYFMVNVTDAISSGDDEDDTDGAEDFVSENSNNKRAPYWTNTEKMEKRLHAVPAANTVKFRCPAGGNPMPTMRWLKNGKEFKQEHRIGGYKVRNQHWSLIMESVVPSDKGNYTCVVENEYGSINHTYHLDVVERSPHRPILQAGLPANASTVVGGDVEFVCKVYSDAQPHIQWIKHVEKNGSKYGPDGLPYLKVLKAAGVNTTDKEIEVLYIRNVTFEDAGEYTCLAGNSIGISFHSAWLTVLPAPGREKEITASPDYLEIAIYCIGVFLIACMVVTVILCRMKNTTKKPDFSSQPAVHKLTKRIPLRRQVTVSAESSSSMNSNTPLVRITTRLSSTADTPMLAGVSEYELPEDPKWEFPRDKLTLGKPLGEGCFGQVVMAEAVGIDKDKPKEAVTVAVKMLKDDATEKDLSDLVSEMEMMKMIGKHKNIINLLGACTQDGPLYVIVEYASKGNLREYLRARRPPGMEYSYDINRVPEEQMTFKDLVSCTYQLARGMEYLASQKCIHRDLAARNVLVTENNVMKIADFGLARDINNIDYYKKTTNGRLPVKWMAPEALFDRVYTHQSDVWSFGVLMWEIFTLGGSPYPGIPVEELFKLLKEGHRMDKPANCTNELYMMMRDCWHAVPSQRPTFKQLVEDLDRILTLTTNEEYLDLSQPLEQYSPSYPDTRSSCSSGDDSVFSPDPMPYEPCLPQYPHINGSVKT'
params.ligand_smiles  = 'CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC'
params.sample_id      = 'fgfr2_20lipids'
params.openmm_input_name = '012_FGFR2_clean_linkers.pdb'

process MAKE_BOLTZ_INPUT {
    tag "${params.sample_id}"
    publishDir "${params.outdir}/01_boltz_input", mode: 'copy'

    output:
    tuple val(params.sample_id), path('fgfr2_20lipids.yaml')

    script:
    def ligands = (['B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U'])
    def seq = params.sequence
    def smiles = params.ligand_smiles

    def seq_block = """
version: 1
sequences:
  - protein:
      id: A
      msa: empty
      sequence: "${seq}"
"""

    def ligand_blocks = ligands.collect { x ->
        """  - ligand:
      id: ${x}
      smiles: "${smiles}"
"""
    }.join("")

    def constraint_blocks = """constraints:
  - pocket:
      binder: B
      contacts: &id001
      - - A
        - 378
      - - A
        - 398
      max_distance: 6.0
""" + ligands.findAll { it != 'B' }.collect { x ->
        """  - pocket:
      binder: ${x}
      contacts: *id001
      max_distance: 6.0
"""
    }.join("")

    """
    cat > fgfr2_20lipids.yaml <<'EOF'
${seq_block}${ligand_blocks}${constraint_blocks}
EOF
    """
}

process RUN_BOLTZ2 {
    tag "$sample_id"
    publishDir "${params.outdir}/02_boltz2", mode: 'copy'

    input:
    tuple val(sample_id), path(yaml_file)

    output:
    tuple val(sample_id), path('fgfr2_20lipids_model_0.cif'), emit: cif
    tuple val(sample_id), path('boltz_out'), emit: boltz_dir

    script:
    """
    module purge
    module load python/3.11.6-gcc-13.2.0 cuda/12.2.1-gcc-13.2.0
    source ${params.boltz_venv}/bin/activate

    export PATH="/usr/local/cuda/12.2/bin:\$PATH"
    export LD_LIBRARY_PATH="/usr/local/cuda/12.2/lib64\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
    export TORCH_HOME="${params.torch_cache}"
    export HF_HOME="${params.hf_cache}"

    mkdir -p boltz_out
    boltz predict ${yaml_file} \
      --out_dir boltz_out \
      --use_msa_server \
      --accelerator gpu --devices 1

    find boltz_out -maxdepth 3 -type f

    cif_file=\$(find boltz_out -name "*.cif" | head -n 1)
    if [ -z "\$cif_file" ]; then
      echo "ERROR: no CIF file produced by boltz" >&2
      exit 1
    fi

    cp "\$cif_file" fgfr2_20lipids_model_0.cif
    test -s fgfr2_20lipids_model_0.cif
    """
}

process RUNCHIMERAX {
    tag "$sample_id"
    publishDir params.outdir + "/03chimerax", mode: 'copy'

    input:
    tuple val(sample_id), path(ciffile)

    output:
    tuple val(sample_id), path("013_protein_preMD.pdb")

    script:
    """
    export HOME=\$PWD
    export XDG_CONFIG_HOME=\$PWD/.config
    export XDG_DATA_HOME=\$PWD/.local/share
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$HOME" "\$XDG_CONFIG_HOME" "\$XDG_DATA_HOME" "\$XDG_CACHE_HOME"

    cat > run.cxc <<EOF
    # subdomains & linkers
    #1:A:25-125      # sub1
    #1:A:126-153     # linker12
    #1:A:154-247     # sub2
    #1:A:248-255     # linker23
    #1:A:256-358     # sub3
    #1:A:481-670     # ICD
    #1:A:378-398     # TM

    close session
    open ${ciffile}
    cartoon
    color /A:1-377 blue
    color /A:378-398 red  
    color /A:399-821 green

    show #1/B-U
    style #1/B-U stick
    color #1/B-U gold

    measure inertia #1/A:378-398@CA

    shape cylinder fromPoint 7.4935,11.63,9.6412 axis -0.582,0.739,-0.338 height 40 radius 0.3 color red coordinateSystem #1
    shape cylinder fromPoint 7.4935,11.63,9.6412 axis 0.582,-0.739,0.338 height 40 radius 0.3 color red coordinateSystem #1

    # ICD
    torsion #1:A:408@C #1:A:409@N #1:A:409@CA #1:A:409@C 90 move large
    torsion #1:A:399@C #1:A:400@N #1:A:400@CA #1:A:400@C 150 move large
    torsion #1:A:470@C #1:A:471@N #1:A:471@CA #1:A:471@C 10 move small
    torsion #1:A:419@C #1:A:420@N #1:A:420@CA #1:A:420@C -8 move small
    torsion #1:A:420@C #1:A:421@N #1:A:421@CA #1:A:421@C -10 move small
    torsion #1:A:418@C #1:A:419@N #1:A:419@CA #1:A:419@C 10 move small

    getcrd #1:A:401@CA
    turn 0,1,0 -60 atoms #1:A:401-412 center 18.061,-0.777,8.716 coordinateSystem #1

    # ECD
    turn 0.582,-0.739,0.338 8 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:361-377
    turn 0.582,-0.739,0.338 10 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:256-360
    turn 0.582,-0.739,0.338 10 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:256-360
    turn -0.582,0.739,-0.338 5 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:256-360

    turn 0.582,-0.739,0.338 12 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:154-255
    turn 0.582,-0.739,0.338 8 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:154-255
    turn 0.582,-0.739,0.338 10 center 7.4935,11.63,9.6412 coordinateSystem #1 atoms #1/A:25-153

    torsion #1/A:250@C #1/A:251@N #1/A:251@CA #1/A:251@C 10 move small
    torsion #1/A:253@C #1/A:254@N #1/A:254@CA #1/A:254@C -10 move small

    torsion #1/A:376@C #1/A:377@N #1/A:377@CA #1/A:377@C 30 move small

    torsion #1/A:152@C #1/A:153@N #1/A:153@CA #1/A:153@C 20 move small
    torsion #1/A:154@C #1/A:155@N #1/A:155@CA #1/A:155@C 15 move small

    save "\$PWD/013_protein_preMD.pdb" models #1
    exit
    EOF

    chimerax --nogui --script "\$PWD/run.cxc"
    ls -l 013_protein_preMD.pdb
    test -s 013_protein_preMD.pdb
    """
}

process STRAIGHTEN_LINKERS {
    tag "$sample_id"
    publishDir "${params.outdir}/04_biopython", mode: 'copy'

    input:
    tuple val(sample_id), path(pre_md_pdb)

    output:
    tuple val(sample_id), path('012_FGFR2_clean_linkers.pdb')

    script:
    """
    source "\$(conda info --base)/etc/profile.d/conda.sh"
    conda activate ${params.biopython_env}

    cp ${pre_md_pdb} 012_protein_preMD.pdb

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
s = parser.get_structure('FGFR2', '012_protein_preMD.pdb')
model = s[0]
chainA = model['A']
linkers = [
    ('ECD_linker', (359, 377), 18.0, 0.2),
    ('ICD_linker', (399, 418), 25.0, 0.3)
]
for name, range_, target_len, bend in linkers:
    stretch_linker_selective(chainA, range_, target_len, bend)
    print(f"{name} completed")

output = '012_FGFR2_clean_linkers.pdb'
io = PDBIO()
io.set_structure(s)
io.save(output)
print(f"Saved {output} ({os.path.getsize(output)/1024:.0f}KB)")
PYEOF

    python straighten_linker.py
    test -s 012_FGFR2_clean_linkers.pdb
    """
}

process OPENMM_MINIMIZE {
    tag "$sample_id"
    publishDir "${params.outdir}/05_openmm", mode: 'copy'

    input:
    tuple val(sample_id), path(clean_pdb)

    output:
    tuple val(sample_id), path('012_3_fgfr2_clean_minimized.pdb')
    tuple val(sample_id), path('012_1_clean_fixed.pdb')
    tuple val(sample_id), path('012_2_clean_withH.pdb')

    script:
    """
    module purge
    module load anaconda3/2022.10-gcc-13.2.0
    eval "\$(conda shell.bash hook)"
    conda activate ${params.openmm_env}

    cat > fgfr2_minimize_clean.py << 'PYEOF'
import openmm
import openmm.app
import openmm.unit as unit
import pdbfixer

print('=== FGFR2 Clean Linkers Minimization ===')
fixer = pdbfixer.PDBFixer(filename='012_FGFR2_clean_linkers.pdb')
fixer.removeHeterogens(False)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
print(f'Repaired: {fixer.topology.getNumResidues()} residues')
openmm.app.PDBFile.writeFile(fixer.topology, fixer.positions, open('012_1_clean_fixed.pdb', 'w'))

pdb = openmm.app.PDBFile('012_1_clean_fixed.pdb')
modeller = openmm.app.Modeller(pdb.topology, pdb.positions)
forcefield = openmm.app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

disulfide_count = 0
for bond in list(modeller.topology.bonds()):
    a1, a2 = bond
    if a1.name == 'SG' and a2.name == 'SG':
        modeller.topology._bonds.remove(bond)
        disulfide_count += 1
print(f'Removed {disulfide_count} bad disulfides')

modeller.addHydrogens(forcefield, pH=7.4)
openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, open('012_2_clean_withH.pdb', 'w'))

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=openmm.app.CutoffNonPeriodic,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=openmm.app.HBonds,
    removeCMMotion=True
)

integrator = openmm.LangevinMiddleIntegrator(
    300*unit.kelvin,
    1/unit.picosecond,
    0.002*unit.picoseconds
)

simulation = openmm.app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

state = simulation.context.getState(getEnergy=True)
print(f"Initial energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.1f} kJ/mol")

print('Coarse min...')
simulation.minimizeEnergy(tolerance=100*unit.kilojoule_per_mole/unit.nanometer, maxIterations=500)

state = simulation.context.getState(getEnergy=True)
print(f"After coarse: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.1f} kJ/mol")

print('Fine min...')
simulation.minimizeEnergy(tolerance=1*unit.kilojoule_per_mole/unit.nanometer, maxIterations=2000)

state = simulation.context.getState(getPositions=True, getEnergy=True)
openmm.app.PDBFile.writeFile(simulation.topology, state.getPositions(), open('012_3_fgfr2_clean_minimized.pdb', 'w'))
print(f"Saved minimized structure, final energy = {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.1f} kJ/mol")
PYEOF

    python fgfr2_minimize_clean.py
    test -s 012_3_fgfr2_clean_minimized.pdb
    """
}

workflow {
    boltz_input = MAKE_BOLTZ_INPUT()

    boltz_ch = RUN_BOLTZ2(boltz_input)
    chim_ch  = RUNCHIMERAX(boltz_ch.cif)
    clean_ch = STRAIGHTEN_LINKERS(chim_ch)
    OPENMM_MINIMIZE(clean_ch)
}
