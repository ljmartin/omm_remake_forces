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
def replaceAngles(system):
    forces = system.getForces()
    for c, f in enumerate(forces):
        if isinstance(f, HarmonicAngleForce):  
            energy_expression = "0.5*k*(theta-theta0)^2;"      
            new_angle_force = CustomAngleForce(energy_expression)          
            new_angle_force.addPerAngleParameter("k")                                                                                          
            new_angle_force.addPerAngleParameter("theta0")   
            for angle in range(f.getNumAngles()):
                
                a1,a2,a3, theta0, k = f.getAngleParameters(angle)
                new_angle_force.addAngle(a1, a2, a3, [k,theta0])
            system.addForce(new_angle_force)
        
    for c, f in enumerate(forces):
        if isinstance(f, HarmonicAngleForce):
            system.removeForce(c)
            

system = boilerplate.makeSystem(pdb)
replaceAngles(system)
simulation = boilerplate.makeSimulation(system, pdb)
print('Total force using custom forces:\n\t', simulation.context.getState(getEnergy=True,).getPotentialEnergy())

