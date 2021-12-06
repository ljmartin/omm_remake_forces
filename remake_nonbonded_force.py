from openmm.app import *
from openmm import *
from openmm.unit import *
import boilerplate


pdb = PDBFile('./data/input.pdb')

#make system with regular forces,
#and print potential energy:
system = boilerplate.makeSystem(pdb)


for c, f in enumerate(system.getForces()):
    if isinstance(f, NonbondedForce):
        f.setUseDispersionCorrection(False)

        simulation = boilerplate.makeSimulation(system, pdb)
print('Total force using OpenMM-derived forces:\n\t', simulation.context.getState(getEnergy=True,).getPotentialEnergy())



#make system with CustomNonbondedForce + CustomBondedForce (for exceptions)
#and print potential energy:
def replaceNonbonded(system):
    forces = system.getForces()
    for c, f in enumerate(forces):
        if isinstance(f, NonbondedForce):
            original_nbforce = f
            ONE_4PI_EPS0 = 138.935456
            epsilon_solvent = original_nbforce.getReactionFieldDielectric()
            r_cutoff = original_nbforce.getCutoffDistance()
            k_rf = r_cutoff**(-3) * ((epsilon_solvent - 1) / (2*epsilon_solvent + 1))
            c_rf = r_cutoff**(-1) * ((3*epsilon_solvent) / (2*epsilon_solvent + 1))

            energy_expression = "electrostatics+sterics;"
            
            energy_expression += "electrostatics=ONE_4PI_EPS0*chargeprod*(r^(-1) + k_rf*r^2-c_rf);"
            energy_expression += "chargeprod = charge1*charge2;"
            energy_expression += "k_rf = %f;" % (k_rf.value_in_unit_system(md_unit_system))
            energy_expression += "c_rf = %f;" % (c_rf.value_in_unit_system(md_unit_system))
            energy_expression += "ONE_4PI_EPS0 = %f;" % ONE_4PI_EPS0 
            
            energy_expression += "sterics=4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
            energy_expression += "epsilon = sqrt(epsilon1*epsilon2);"
            energy_expression += "sigma = 0.5*(sigma1+sigma2);"
            
            new_custom_nonbonded_force = openmm.CustomNonbondedForce(energy_expression)
            new_custom_nonbonded_force.addPerParticleParameter('charge')
            new_custom_nonbonded_force.addPerParticleParameter('sigma')
            new_custom_nonbonded_force.addPerParticleParameter('epsilon')
            
            new_custom_nonbonded_force.setNonbondedMethod(original_nbforce.getNonbondedMethod())
            
            new_custom_nonbonded_force.setCutoffDistance(original_nbforce.getCutoffDistance())
            new_custom_nonbonded_force.setUseLongRangeCorrection(False) 
            
            energy_expression = "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod/r;"
            energy_expression += "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)  # already in OpenMM units
            new_custom_bond_force = openmm.CustomBondForce(energy_expression)
            new_custom_bond_force.addPerBondParameter('chargeprod')
            new_custom_bond_force.addPerBondParameter('sigma')
            new_custom_bond_force.addPerBondParameter('epsilon')
            
            for index in range(system.getNumParticles()):
                [charge, sigma, epsilon] = original_nbforce.getParticleParameters(index)
                new_custom_nonbonded_force.addParticle([charge, sigma, epsilon])
                
            for index in range(original_nbforce.getNumExceptions()):
                idx, jdx, c, s, eps = original_nbforce.getExceptionParameters(index)
                new_custom_nonbonded_force.addExclusion(idx, jdx)
                c_value = c/elementary_charge**2
                eps_value = eps/(kilojoule/mole)
                if c_value != 0 or eps_value!=0:
                    new_custom_bond_force.addBond(idx, jdx, [c, s, eps])

            system.addForce(new_custom_nonbonded_force)
            system.addForce(new_custom_bond_force)

            
    for c, f in enumerate(forces):
        if isinstance(f, NonbondedForce):
            system.removeForce(c)
        


system = boilerplate.makeSystem(pdb)
replaceNonbonded(system)
simulation = boilerplate.makeSimulation(system, pdb)
print('Total force using custom forces:\n\t', simulation.context.getState(getEnergy=True,).getPotentialEnergy())

