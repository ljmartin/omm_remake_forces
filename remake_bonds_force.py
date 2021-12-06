from openmm.app import *
from openmm import *
from openmm.unit import *
import boilerplate


pdb = PDBFile('./data/input.pdb')

#make system with regular forces,
#and print potential energy:
system = boilerplate.makeSystem(pdb)
simulation = boilerplate.makeSimulation(system, pdb)
print('Total force using OpenMM-derived forces:\n\t', simulation.context.getState(getEnergy=True,).getPotentialEnergy())



#make system with CustomBondForce,
#and print potential energy:
def replaceBonds(system):
    forces = system.getForces()
    for c, f in enumerate(forces):
        if isinstance(f, HarmonicBondForce):  
            energy_expression = "0.5*k*(r-r0)^2;"      
            new_bond_force = CustomBondForce(energy_expression)          
            new_bond_force.addPerBondParameter("k")                                                                                          
            new_bond_force.addPerBondParameter("r0")   
            for bond in range(f.getNumBonds()):
                
                pars = f.getBondParameters(bond)
                new_bond_force.addBond(pars[0], pars[1], [pars[3], pars[2]])
            system.addForce(new_bond_force)
        
    for c, f in enumerate(forces):
        if isinstance(f, HarmonicBondForce):
            system.removeForce(c)

system = boilerplate.makeSystem(pdb)
replaceBonds(system)
simulation = boilerplate.makeSimulation(system, pdb)
print('Total force using custom forces:\n\t', simulation.context.getState(getEnergy=True,).getPotentialEnergy())

