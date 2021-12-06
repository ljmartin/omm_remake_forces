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



#make system with CustomRemakeForce,
#and print potential energy:
def replaceTorsions(system):
    forces = system.getForces()
    for c, f in enumerate(forces):
        if isinstance(f, PeriodicTorsionForce):  
            energy_expression = "k*(1+cos(periodicity*theta-theta0))"
            new_torsion_force = CustomTorsionForce(energy_expression)         
            new_torsion_force.addPerTorsionParameter("k");
            new_torsion_force.addPerTorsionParameter("periodicity")
            new_torsion_force.addPerTorsionParameter("theta0"); 
    
            for torsion_index in range(f.getNumTorsions()):
                a0, a1, a2, a3, periodicity, phase, k = f.getTorsionParameters(torsion_index)
                new_torsion_force.addTorsion(a0, a1, a2, a3, [k, periodicity,phase])
            system.addForce(new_torsion_force)
            
    for c, f in enumerate(forces):
        if isinstance(f, PeriodicTorsionForce):
            system.removeForce(c)

system = boilerplate.makeSystem(pdb)
replaceTorsions(system)
simulation = boilerplate.makeSimulation(system, pdb)
print('Total force using custom forces:\n\t', simulation.context.getState(getEnergy=True,).getPotentialEnergy())

