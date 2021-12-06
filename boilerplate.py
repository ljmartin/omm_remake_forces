from openmm.app import *
from openmm import *
from openmm.unit import *


def makeSystem(pdb):
    """
    Generate a `System` object
    Note:
        - I use CutoffPeriodic here - this uses Reaction Field for 
          long range forces, not PME. Reason for this is that it's easier
          to write a CustomNonbondedForce for reaction field. 
    """
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=0.9*nanometer
    )
    return system

def makeSimulation(system, pdb):
    """
    Generate an OpenMM `Simulation` object
    Note:
       - I use the `Reference` implementation here for simplicity.
         OpenCL or CUDA are preferred, particularly when replacing
         the nonbonded forces - CPU impl. is slow in this case. 
    """
    integrator = openmm.LangevinIntegrator(
        298*kelvin,
        1/picosecond,
        1*femtosecond
    )
    platform = openmm.Platform.getPlatformByName('Reference')
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    return simulation
