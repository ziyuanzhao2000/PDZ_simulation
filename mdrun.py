from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
#from pdbfixer.pdbfixer import PDBFixer
#fixer = PDBFixer(pdbfile='1BFE')
#fixer.findMissingAtoms()
#fixer.addMissingAtoms()
#PDBFile.writeFile(fixer.topology, fixer.positions, open('output.pdb', 'w'))

pdb = PDBFile('1BFE_cleaned.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.deleteWater()
modeller.addSolvent(forcefield, padding=1.0*nanometers)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

platform = Platform.getPlatformByName('CUDA')

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output_long.pdb', 10000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(100000)
