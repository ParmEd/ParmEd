Seeding OpenMM Simulations with PyRosetta
==========================================

Generate pose from primary structure
---

    from rosetta import init, pose_from_sequence

    init()

    seq = 12*'A'
    pose = pose_from_sequence(seq)


Generate pose an exisiting PDB
---

    from rosetta import init, pose_from_pdb

    init()

    pdb = 'dodecaalanine.pdb'
    pose = pose_from_sequence(seq)


Load Pose into ParmEd
---

    from chemistry.rosetta import load_rosetta

    struct = load_rosetta(pose)


Start an OpenMM simulation
---
    from simtk.openmm import *
    from simtk.openmm.app import *
    from simtk.unit import *

    ff = ForceField('amber99sbildn.xml', 'tip3p.xml')

    #solvate system
    mod = Modeller(struct.topology, struct.positions)
    mod.addSolvent(ff, model='tip3p', boxSize=Vec3(4, 4, 4)*nanometers,
                   positiveIon='Na+', negativeIon='Cl-',
                   ionicStrength=smolar*molar)

    system = ff.createSystem(mod.topology, nonbondedMethod=PME,
                             nonbondedCutoff=10*nanometers,
                             constraints=HBonds)
    integrator = LangevinIntegrator(300*kelvin, 1/picoseconds,
                                    2*femtoseconds)

    integrator.setConstraintTolerance(1e-5)

    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
    simulation = Simulation(mod.topology, system, integrator, platform, props)
    simulation.context.setPositions(positions)

    simulation.minimizeEnergy(maxIterations=1000)

    simulation.reporters.append(
        PDBReporter('dodecaalanine.solv.pdb', 50000))

    simulation.context.setVelocitiesToTemperature(300)

    simulation.step(50000)
