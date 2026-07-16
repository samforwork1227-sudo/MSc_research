import os
import openmm
import openmm.app as app
import openmm.unit as unit
from openmm import XmlSerializer

print("=== FGFR1 membrane minimization: FINE stage ===")

INPUT_PDB = os.environ["INPUT_PDB"]
SYSTEM_XML = os.environ["SYSTEM_XML"]
OUTDIR = os.environ["OUTDIR"]

print(f"Input PDB: {INPUT_PDB}")
print(f"System XML: {SYSTEM_XML}")
print(f"Output dir: {OUTDIR}")

# Step 1: load coarse-minimized PDB + same System XML
pdb = app.PDBFile(INPUT_PDB)
modeller = app.Modeller(pdb.topology, pdb.positions)
print(f"Loaded coarse-minimized system: {modeller.topology.getNumAtoms()} atoms")

with open(SYSTEM_XML) as f:
    system = XmlSerializer.deserialize(f.read())
print(f"Loaded System from {SYSTEM_XML}")

# Step 2: weaker restraints on protein heavy atoms
restraint = openmm.CustomExternalForce(
    "0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
)
restraint.addGlobalParameter("k", 100.0)
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

positions = modeller.positions
for atom in modeller.topology.atoms():
    if atom.residue.name not in ('POPC', 'POPE', 'DPPC', 'HOH', 'WAT', 'TIP3', 'TIP3P', 'Na+', 'Cl-'):
        if atom.element.symbol != app.element.hydrogen.symbol:
            idx = atom.index
            pos = positions[idx]
            restraint.addParticle(idx, pos)

system.addForce(restraint)

integrator = openmm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1 / unit.picosecond,
    0.002 * unit.picoseconds
)

simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Step 3: fine minimisation
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Energy before fine: {energy:.1f} kJ/mol")

print("Fine min with weaker restraints (k=100, maxIterations=1000)...")
simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,
    maxIterations=1000
)

state = simulation.context.getState(getEnergy=True, getPositions=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
positions = state.getPositions()
print(f"Final energy after fine: {energy:.1f} kJ/mol")

# Step 4: save final minimized structures
full_pdb = os.path.join(OUTDIR, '6_fgfr1_membrane_minimized_full.pdb')
app.PDBFile.writeFile(
    simulation.topology,
    positions,
    open(full_pdb, 'w')
)
print(f"Saved minimized full system to {full_pdb}")

slim_min = app.Modeller(simulation.topology, positions)
waters_ions_min = [res for res in slim_min.topology.residues()
                   if res.name in ('HOH', 'WAT', 'TIP3', 'TIP3P', 'Na+', 'Cl-')]
slim_min.delete(waters_ions_min)

no_water_pdb = os.path.join(OUTDIR, '6_fgfr1_membrane_minimized_no_water.pdb')
app.PDBFile.writeFile(
    slim_min.topology, slim_min.positions, open(no_water_pdb, 'w')
)
print(f"Saved minimized protein+membrane to {no_water_pdb}")
