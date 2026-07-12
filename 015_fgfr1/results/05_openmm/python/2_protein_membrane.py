import os
import openmm
import openmm.app as app
import openmm.unit as unit
from openmm import XmlSerializer

print("=== FGFR1 membrane build + minimization ===")

INPUT_PDB = os.environ["INPUT_PDB"]
OUTDIR = os.environ["OUTDIR"]

TM_START = 378
TM_END = 398
TARGET_CHAIN = None
PH = 7.4

print(f"Input file: {INPUT_PDB}")
print(f"Output dir: {OUTDIR}")

############################################
# Step 1: Build modeller
############################################

pdb = app.PDBFile(INPUT_PDB)
modeller = app.Modeller(pdb.topology, pdb.positions)

############################################
# Step 2: Put residues 378-398 near membrane center (z=0)
############################################

forcefield = app.ForceField(
    'charmm36.xml',
    'charmm36/water.xml'
)

tm_atom_indices = []
tm_residue_ids = []

for chain in modeller.topology.chains():
    if TARGET_CHAIN is not None and chain.id != TARGET_CHAIN:
        continue
    for residue in chain.residues():
        try:
            resid = int(residue.id)
        except ValueError:
            continue
        if TM_START <= resid <= TM_END:
            tm_residue_ids.append((chain.id, residue.id, residue.name))
            for atom in residue.atoms():
                if atom.name == 'CA':
                    tm_atom_indices.append(atom.index)

if len(tm_atom_indices) == 0:
    raise ValueError(f"No CA atoms found in residues {TM_START}-{TM_END}. Check residue numbering/chain ID.")

positions_nm = modeller.positions.value_in_unit(unit.nanometer)

z_vals = [positions_nm[i][2] for i in tm_atom_indices]
tm_center_z = sum(z_vals) / len(z_vals)

print(f"TM residues found: {len(tm_residue_ids)}")
print(f"TM center z before shift: {tm_center_z:.3f} nm")

shifted_positions = [
    openmm.Vec3(pos[0], pos[1], pos[2] - tm_center_z)
    for pos in positions_nm
]

modeller = app.Modeller(modeller.topology, shifted_positions)

print("Shifted protein so residues 378-398 are centered near z=0")

############################################
# Step 3: Add hydrogens and membrane
############################################

print("Adding membrane...")
modeller.addMembrane(
    forcefield,
    lipidType='POPC',
    membraneCenterZ=0.0 * unit.nanometer,
    minimumPadding=0.8 * unit.nanometer,
    positiveIon='Na+',
    negativeIon='Cl-',
    ionicStrength=0.15 * unit.molar
)

lipid_names = ('POPC', 'POPE', 'DPPC')
solvent_names = ('HOH', 'WAT', 'TIP3', 'TIP3P', 'Na+', 'Cl-')

positions_nm = modeller.positions.value_in_unit(unit.nanometer)

protein_atom_indices = [
    atom.index for atom in modeller.topology.atoms()
    if atom.residue.name not in lipid_names + solvent_names
]

x_vals = [positions_nm[i][0] for i in protein_atom_indices]
y_vals = [positions_nm[i][1] for i in protein_atom_indices]
x_min, x_max = min(x_vals), max(x_vals)
y_min, y_max = min(y_vals), max(y_vals)

margin = 2.5

x_min_box = x_min - margin
x_max_box = x_max + margin
y_min_box = y_min - margin
y_max_box = y_max + margin

to_delete = []
for res in modeller.topology.residues():
    if res.name in lipid_names + solvent_names:
        in_box = False
        for atom in res.atoms():
            pos = positions_nm[atom.index]
            if (x_min_box <= pos[0] <= x_max_box) and (y_min_box <= pos[1] <= y_max_box):
                in_box = True
                break
        if not in_box:
            to_delete.append(res)

print(f"Deleting {len(to_delete)} lipid/water residues far from the protein")
modeller.delete(to_delete)

with_membrane_pdb = os.path.join(OUTDIR, '4_with_membrane.pdb')
app.PDBFile.writeFile(
    modeller.topology, modeller.positions, open(with_membrane_pdb, 'w')
)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=5e-4
)

system_xml = os.path.join(OUTDIR, '4_with_membrane_system.xml')
with open(system_xml, 'w') as f:
    f.write(XmlSerializer.serialize(system))

print(f"Wrote {system_xml}")

slim = app.Modeller(modeller.topology, modeller.positions)

waters_ions = [res for res in slim.topology.residues()
               if res.name in ('HOH', 'WAT', 'TIP3', 'TIP3P', 'Na+', 'Cl-')]
slim.delete(waters_ions)

membrane_no_water_pdb = os.path.join(OUTDIR, '4_membrane_no_water.pdb')
app.PDBFile.writeFile(
    slim.topology, slim.positions, open(membrane_no_water_pdb, 'w')
)

print(f"Saved {with_membrane_pdb}")
print(f"Saved {membrane_no_water_pdb}")
