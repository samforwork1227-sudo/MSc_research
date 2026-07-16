
import openmm
import openmm.app
import openmm.unit
import pdbfixer

print("=== FGFR2 structure preparation ===")

############################################
# Step 1  PDBFixer repair
############################################

fixer = pdbfixer.PDBFixer(filename='009_protein_preMD.pdb')

fixer.removeHeterogens(False)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# remove problematic disulfide bonds
for bond in list(fixer.topology.bonds()):
    a1, a2 = bond
    if a1.element.symbol == 'S' and a2.element.symbol == 'S':
        fixer.topology._bonds.remove(bond)

print("Protein repaired")

############################################
# Save repaired structure
############################################

openmm.app.PDBFile.writeFile(
    fixer.topology,
    fixer.positions,
    open('fgfr2_fixed.pdb','w')
)

print("Fixed structure saved")

############################################
# Step 2  Load structure
############################################

pdb = openmm.app.PDBFile('fgfr2_fixed.pdb')
modeller = openmm.app.Modeller(pdb.topology, pdb.positions)

############################################
# Step 3  Forcefield
############################################

forcefield = openmm.app.ForceField(
    'amber99sb.xml',
    'implicit/obc2.xml'
)

############################################
# Step 4  Add hydrogens
############################################

modeller.addHydrogens(forcefield, pH=7.0)

openmm.app.PDBFile.writeFile(
    modeller.topology,
    modeller.positions,
    open('fgfr2_withH.pdb','w')
)

print("Hydrogens added")

############################################
# Step 5  Build system
############################################

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=openmm.app.NoCutoff,
    constraints=openmm.app.HBonds,
    ignoreExternalBonds=True
)

############################################
# Step 6  Integrator
############################################

integrator = openmm.LangevinMiddleIntegrator(
    300*openmm.unit.kelvin,
    1/openmm.unit.picosecond,
    0.004*openmm.unit.picoseconds
)

############################################
# Step 7  Simulation
############################################

simulation = openmm.app.Simulation(
    modeller.topology,
    system,
    integrator
)

simulation.context.setPositions(modeller.positions)

############################################
# Step 8  Energy minimization
############################################

print("Starting minimization")

simulation.minimizeEnergy(maxIterations=500)

############################################
# Step 9  Save structure
############################################

positions = simulation.context.getState(getPositions=True).getPositions()

openmm.app.PDBFile.writeFile(
    simulation.topology,
    positions,
    open('fgfr2_minimized.pdb','w')
)

print("Minimization finished")

