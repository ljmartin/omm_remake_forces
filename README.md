# omm_remake_forces
remaking OpenMM forces with Custom\*Forces

This repo demonstrates how to replace the 'native' forces in OpenMM - `HarmonicBondForce`, `HarmonicAngleForce`, `PeriodicTorsionForce`, and `NonbondedForce` - with Custom\*Force versions. To run one of these, simply run `python remake_bonds_force.py` (or any of the other force types). This generates `System` objects with or without a Custom\*Force and prints the energies to demonstrate they are equivalent. The input file is a PDB of protein in solvent+ions from the `share` installed with OpenMM (typically run using `simulatePdb.py`).

Some of the code is hidden away in `boilerplate.py`. This is just for readability. The code is just to instantiate the `System` and `Simulation` objects and so it's relatively uninteresting:
```python
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
```




Notes:
- :rotating_light: `remake_nonbonded_force.py` currently turns OFF the long-range dispersion correction. For the reason why, see: https://github.com/openmm/openmm/issues/3281#issuecomment-940593435 . 
- All `System` objects are generted with `app.CutoffPeriodic`, which applies reaction field, rather than PME, for long-range electrostatics. For the reason why, see this comment: https://github.com/openmm/openmm/issues/3281#issuecomment-940551386
- Use of CustomNonbondedForce is slow on the CPU platform. For the reason why, see this comment: https://github.com/openmm/openmm/issues/3162#issuecomment-871605371




