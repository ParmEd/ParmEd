#!/usr/bin/env python

"""
Test ParmEd's ability to process a Gromacs position/topology file
by comparing Gromacs energy/force to OpenMM-via-ParmEd energy/force.

This script contains bits of ForceBalance to obtain the Gromacs energy/force
and also reads parts of the Gromacs .mdp file to set up the system.

There are also some OpenMM imports for calculating the OpenMM energy/force..

To run this script, provide a gro, top and mdp file.  The difference from
test.py is that this script also runs AMBER.

Author: Lee-Ping Wang
"""

# General import
from collections import OrderedDict
import numpy as np
import os, sys, re, copy
import argparse

# ForceBalance convenience functions
from nifty import printcool, printcool_dictionary, _exec, which, wopen, isint, isfloat, logger
# Only needed for writing constrained .gro files
# from molecule import Molecule

# ParmEd import
from parmed import gromacs, amber
from parmed.amber.mdin import Mdin

# OpenMM import
import simtk.unit as u
import simtk.openmm as mm
import simtk.openmm.app as app

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--amber', action='store_true', help='Pass this flag to run AMBER tests along with Gromacs / OpenMM.')
args, sys.argv = parser.parse_known_args(sys.argv)

# Gromacs settings
gmxsuffix="_d"
if which('gmx'+gmxsuffix) != '':
    logger.info("Using double precision GROMACS version 5\n")
    gmxpath = which('gmx'+gmxsuffix)
    GMXVERSION = 5
elif which('mdrun'+gmxsuffix) != '':
    logger.info("Using double precision GROMACS version 4\n")
    gmxpath = which('mdrun'+gmxsuffix)
    GMXVERSION = 4
else:
    gmxsuffix=""
    if which('gmx'+gmxsuffix) != '':
        logger.info("Using single precision GROMACS version 5\n")
        gmxpath = which('gmx'+gmxsuffix)
        GMXVERSION = 5
    elif which('mdrun'+gmxsuffix) != '':
        logger.info("Using single precision GROMACS version 4\n")
        gmxpath = which('mdrun'+gmxsuffix)
        GMXVERSION = 4
    else:
        logger.error("Cannot find the GROMACS executables!\n")
        raise RuntimeError
os.environ["GMX_MAXBACKUP"] = "-1"
os.environ["GMX_NO_SOLV_OPT"] = "TRUE"
os.environ["GMX_NO_ALLVSALL"] = "TRUE"

gmxprogs = ["anadock", "anaeig", "analyze", "angle", "bar", "bond", "bundle", "chi", 
            "cluster", "clustsize", "confrms", "covar", "current", "enemat", "energy", 
            "filter", "gyrate", "h2order", "hbond", "helix", "helixorient", "hydorder", 
            "kinetics", "lie", "luck", "mdmat", "membed", "mindist", "morph", "msd", 
            "nmeig", "nmens", "nmtraj", "options", "order", "pme_error", "polystat", 
            "potential", "principal", "protonate", "rama", "rdf", "rms", "rmsdist", 
            "rmsf", "rotacf", "rotmat", "saltbr", "sans", "sas", "select", "sgangle", 
            "sham", "sigeps", "sorient", "spatial", "spol", "tcaf", "traj", "tune_pme", 
            "vanhove", "velacc", "wham", "wheel", "x2top"]

def edit_mdp(fin=None, fout=None, options={}, defaults={}, verbose=False):
    """
    Read, create or edit a Gromacs MDP file.  The MDP file contains GROMACS run parameters.
    If the input file exists, it is parsed and options are replaced where "options" overrides them.
    If the "options" dictionary contains more options, they are added at the end.
    If the "defaults" dictionary contains more options, they are added at the end.
    Keys are standardized to lower-case strings where all dashes are replaced by underscores.
    The output file contains the same comments and "dressing" as the input.
    Also returns a dictionary with the final key/value pairs.

    Parameters
    ----------
    fin : str, optional
        Input .mdp file name containing options that are more important than "defaults", but less important than "options"
    fout : str, optional
        Output .mdp file name.
    options : dict, optional
        Dictionary containing mdp options. Existing options are replaced, new options are added at the end, None values are deleted from output mdp.
    defaults : dict, optional
        defaults Dictionary containing "default" mdp options, added only if they don't already exist.
    verbose : bool, optional
        Print out additional information        

    Returns
    -------
    OrderedDict
        Key-value pairs combined from the input .mdp and the supplied options/defaults and equivalent to what's printed in the output mdp.
    """
    clashes = ["pbc"]
    # Make sure that the keys are lowercase, and the values are all strings.
    options = OrderedDict([(key.lower().replace('-','_'), str(val) if val is not None else None) for key, val in options.items()])
    # List of lines in the output file.
    out = []
    # List of options in the output file.
    haveopts = []
    # List of all options in dictionary form, to be returned.
    all_options = OrderedDict()
    if fin is not None and os.path.isfile(fin):
        for line in open(fin).readlines():
            line    = line.strip().expandtabs()
            # The line structure should look something like this:
            # key   = value  ; comments
            # First split off the comments.
            if len(line) == 0:
                out.append('')
                continue
            s = line.split(';',1)
            data = s[0]
            comms = s[1] if len(s) > 1 else None
            # Pure comment lines or empty lines get appended to the output.
            if set(data).issubset([' ']):
                out.append(line)
                continue
            # Now split off the key and value fields at the equals sign.
            keyf, valf = data.split('=',1)
            key = keyf.strip().lower().replace('-','_')
            haveopts.append(key)
            if key in options:
                val = options[key]
                val0 = valf.strip()
                if key in clashes and val != val0:
                    logger.error("edit_mdp tried to set %s = %s but its original value was %s = %s\n" % (key, val, key, val0))
                    raise RuntimeError
                # Passing None as the value causes the option to be deleted
                if val is None: continue
                if len(val) < len(valf):
                    valf = ' ' + val + ' '*(len(valf) - len(val)-1)
                else:
                    valf = ' ' + val + ' '
                lout = [keyf, '=', valf]
                if comms is not None:
                    lout += [';',comms]
                out.append(''.join(lout))
            else:
                out.append(line)
                val = valf.strip()
            all_options[key] = val
    for key, val in options.items():
        key = key.lower().replace('-','_')
        if key not in haveopts:
            haveopts.append(key)
            out.append("%-20s = %s" % (key, val))
            all_options[key] = val
    # Fill in some default options.
    for key, val in defaults.items():
        key = key.lower().replace('-','_')
        options[key] = val
        if key not in haveopts:
            out.append("%-20s = %s" % (key, val))
            all_options[key] = val
    if fout != None:
       file_out = wopen(fout) 
       for line in out:
           print >> file_out, line
       file_out.close()
    if verbose:
        printcool_dictionary(options, title="%s -> %s with options:" % (fin, fout))
    return all_options

def rm_gmx_baks(dir):
    # Delete the #-prepended files that GROMACS likes to make
    for root, dirs, files in os.walk(dir):
        for file in files:
            if re.match('^#',file):
                os.remove(file)

def callgmx(command, stdin=None, print_to_screen=False, print_command=False, **kwargs):
    # Remove backup files.
    rm_gmx_baks(os.getcwd())
    # Call a GROMACS program as you would from the command line.
    if GMXVERSION == 5:
        csplit = ('gmx ' + command.replace('gmx', '')).split()
    else:
        if command.split()[0] in gmxprogs:
            csplit = ('g_%s' % command).split()
        else:
            csplit = command.split()
    prog = os.path.join(gmxpath, csplit[0])
    csplit[0] = prog + gmxsuffix
    return _exec(' '.join(csplit), stdin=stdin, print_to_screen=print_to_screen, print_command=print_command, **kwargs)

def energy_termnames(edrfile):
    """ Get a list of energy term names from the .edr file by parsing a system call to g_energy. """
    if not os.path.exists(edrfile):
        logger.error('Cannot determine energy term names without an .edr file\n')
        raise RuntimeError
    ## Figure out which energy terms need to be printed.
    o = callgmx("energy -f %s -xvg no" % (edrfile), stdin="Total-Energy\n", copy_stdout=False, copy_stderr=True)
    parsemode = 0
    energyterms = OrderedDict()
    for line in o:
        s = line.split()
        if "Select the terms you want from the following list" in line:
            parsemode = 1
        if parsemode == 1:
            if len(s) > 0 and all([isint(i) for i in s[::2]]):
                parsemode = 2
        if parsemode == 2:
            if len(s) > 0:
                try:
                    if all([isint(i) for i in s[::2]]):
                        for j in range(len(s))[::2]:
                            num = int(s[j])
                            name = s[j+1]
                            energyterms[name] = num
                except: pass
    return energyterms

def energy_components(Sim, verbose=False):
    # Before using EnergyComponents, make sure each Force is set to a different group.
    EnergyTerms = OrderedDict()
    if type(Sim.integrator) in [mm.LangevinIntegrator, mm.VerletIntegrator]:
        for i in range(Sim.system.getNumForces()):
            EnergyTerms[Sim.system.getForce(i).__class__.__name__] = Sim.context.getState(getEnergy=True,groups=2**i).getPotentialEnergy() / u.kilojoules_per_mole
    EnergyTerms['Potential'] = Sim.context.getState(getEnergy=True).getPotentialEnergy() / u.kilojoules_per_mole
    return EnergyTerms

def interpret_mdp(mdp_file):
    # Keyword args to pass to createSystem()
    sysargs = {}
    # Read stuff from the Gromacs .mdp file
    # to inform how we build the OpenMM System
    mdp_opts = edit_mdp(mdp_file)
    if 'define' in mdp_opts:
        defines = dict([(k.replace("-D",''),1) for k in mdp_opts['define'].split()])
    else:
        defines = {}
    print "Defines:", defines
    sysargs['rigidWater'] = 'FLEXIBLE' not in defines
    # Constraints
    constraint_map = {'none':None,'h-bonds':app.HBonds,'all-bonds':app.AllBonds,'h-angles':app.HAngles}
    if 'constraints' in mdp_opts:
        omm_constraints = constraint_map[mdp_opts['constraints'].replace('_','-').lower()]
    else:
        omm_constraints = None
    print "Constraints", omm_constraints
    sysargs['constraints'] = omm_constraints
    # Periodic boundary conditions
    if mdp_opts['pbc'].lower() in ['none', 'no']:
        pbc = False
    elif mdp_opts['pbc'].lower() == 'xyz':
        pbc = True
    else:
        raise RuntimeError('Unsupported PBC')
    # Cut-off radii and nonbonded method
    if float(mdp_opts['rcoulomb']) != float(mdp_opts['rvdw']):
        raise RuntimeError('Please set rcoulomb to equal rvdw')
    if 'rvdw_switch' in mdp_opts:
        sysargs['switchDistance'] = mdp_opts['rvdw_switch'] * u.nanometer
    
    if mdp_opts['coulombtype'].lower() == 'cut-off':
        if float(mdp_opts['rcoulomb']) == 0.0:
            sysargs['nonbondedMethod'] = app.NoCutoff
        elif pbc:
            sysargs['nonbondedMethod'] = app.CutoffPeriodic
            sysargs['nonbondedCutoff'] = float(mdp_opts['rcoulomb'])*u.nanometer
        else:
            sysargs['nonbondedMethod'] = app.CutoffNonPeriodic
            sysargs['nonbondedCutoff'] = float(mdp_opts['rcoulomb'])*u.nanometer
    elif mdp_opts['coulombtype'].lower() == 'pme':
        sysargs['nonbondedMethod'] = app.PME
        sysargs['ewaldErrorTolerance'] = 1e-5
        sysargs['nonbondedCutoff'] = float(mdp_opts['rcoulomb'])*u.nanometer
    return defines, sysargs, mdp_opts

def Calculate_GMX(gro_file, top_file, mdp_file):
    #===============================#
    #| GROMACS energies and forces |#
    #===============================#
    # Create .mdp file for single-point energies and forces.
    shot_opts = OrderedDict([("nsteps", 0), ("nstxout", 0), ("nstxtcout", 0), ("nstenergy", 1), ("nstfout", 1)])
    edit_mdp(fin=mdp_file, fout="enerfrc.mdp", options=shot_opts)
    # Call grompp to set up calculation.
    callgmx("grompp -f enerfrc.mdp -c %s -p %s -maxwarn 1" % (gro_file, top_file))
    # Run gmxdump to determine which atoms are real.
    o = callgmx("gmxdump -s topol.tpr -sys", copy_stderr=True)
    AtomMask = []
    for line in o:
        line = line.replace("=", "= ")
        if "ptype=" in line:
            s = line.split()
            ptype = s[s.index("ptype=")+1].replace(',','').lower()
            AtomMask.append(ptype=='atom')
    # Get the energy and the forces.
    callgmx("mdrun -nt 1 -rerunvsite -rerun %s" % gro_file)
    callgmx("energy -xvg no -f ener.edr -o energy.xvg", stdin='Potential')
    Efile = open("energy.xvg").readlines()
    GMX_Energy = np.array([float(Eline.split()[1]) for Eline in Efile])
    callgmx("traj -xvg no -s topol.tpr -f traj.trr -of force.xvg -fp", stdin='System')
    GMX_Force = np.array([[float(j) for i, j in enumerate(line.split()[1:]) if AtomMask[i/3]] \
                              for line in open("force.xvg").readlines()])
    # Perform energy component analysis and return properties.
    energyterms = energy_termnames("ener.edr")
    ekeep = [k for k,v in energyterms.items() if v <= energyterms['Total-Energy']]
    callgmx("energy -f ener.edr -o energy.xvg -xvg no", stdin="\n".join(ekeep))
    ecomp = OrderedDict()
    for line in open("energy.xvg"):
        s = [float(i) for i in line.split()]
        for i in range(len(ekeep) - 2):
            val = s[i+1]
            if ekeep[i] in ecomp:
                ecomp[ekeep[i]].append(val)
            else:
                ecomp[ekeep[i]] = [val]
    Ecomps_GMX = OrderedDict([(key, val[0]) for key, val in ecomp.items()])
    return GMX_Energy[0], GMX_Force, Ecomps_GMX

def Calculate_ParmEd_OpenMM(gro_file, top_file, sysargs, defines):
    #===============================#
    #|   ParmEd object creation    |#
    #===============================#
    # Make sure the proper defines from the .mdp file are passed into the GromacsTopologyFile() :)
    ParmEd_GmxTop = gromacs.GromacsTopologyFile(top_file, defines=defines)
    ParmEd_GmxGro = gromacs.GromacsGroFile.parse(gro_file)
    ParmEd_GmxTop.box = ParmEd_GmxGro.box
    ParmEd_GmxTop.positions = ParmEd_GmxGro.positions
    #===============================#
    #|   OpenMM simulation setup   |#
    #===============================#
    # ParmEd creates System object
    system = ParmEd_GmxTop.createSystem(**sysargs)
    # Keep a record of which atoms are real (not virtual sites)
    isAtom = []
    for i in range(system.getNumParticles()):
        isAtom.append(system.getParticleMass(i).value_in_unit(u.dalton) > 0.0)
    # Setting force groups enables energy components analysis
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(True)
        elif isinstance(f, mm.CustomNonbondedForce):
            f.setUseLongRangeCorrection(True)
    integ = mm.VerletIntegrator(1.0*u.femtosecond)
    plat = mm.Platform.getPlatformByName('Reference')
    # Create Simulation object
    simul = app.Simulation(ParmEd_GmxTop.topology, system, integ, plat)
    simul.context.setPositions(ParmEd_GmxGro.positions)
    simul.context.applyConstraints(1e-12)
    # Obtain OpenMM potential energy
    state = simul.context.getState(getPositions=True,getEnergy=True,getForces=True)
    parmed_energy = state.getPotentialEnergy()
    parmed_forces = state.getForces()
    pos = np.array(state.getPositions().value_in_unit(u.angstrom)).reshape(-1,3)
    # Obtain and save constrained positions
    # M = Molecule(gro_file)
    # M.xyzs[0] = pos
    # M.write('constrained.gro')
    # Print OpenMM-via-ParmEd energy components
    Ecomps_OMM = energy_components(simul)
    printcool_dictionary(Ecomps_OMM, title="OpenMM energy components via ParmEd")
    parmed_forces = np.array([f for i, f in enumerate(parmed_forces.value_in_unit(u.kilojoule_per_mole/u.nanometer)) if isAtom[i]])
    return ParmEd_GmxTop, parmed_energy, parmed_forces, Ecomps_OMM

def Calculate_AMBER(Structure, mdp_opts):
    pbc = mdp_opts["pbc"].lower() == "xyz"
    # Create AMBER inpcrd file
    inpcrd = amber.AmberAsciiRestart("inpcrd", mode="w")
    inpcrd.coordinates = np.array(Structure.positions.value_in_unit(u.angstrom)).reshape(-1,3)
    inpcrd.box = Structure.box
    inpcrd.close()
    # sander insists on providing a trajectory to iterate over, 
    # so we feed it the same coordinates again. But we don't use it
    # because the positions are imprecise.
    mdcrd = amber.AmberMdcrd("mdcrd", natom=len(Structure.atoms), hasbox=pbc, mode="w")
    mdcrd.add_coordinates(np.array(Structure.positions.value_in_unit(u.angstrom)).reshape(-1,3))
    if pbc:
        mdcrd.add_box(Structure.box[:3])
    mdcrd.close()
    # Create AMBER prmtop object from ParmEd Structure :)
    prmtop = amber.AmberParm.from_structure(Structure)
    prmtop.write_parm("prmtop")
    # Create AMBER mdin file and append some stuff
    mdin = Mdin()
    # Single point energies?
    mdin.change('cntrl','imin','5')
    # Periodic boundary conditions?
    if pbc:
        mdin.change('cntrl','ntb','1')
    else:
        mdin.change('cntrl','ntb','0')
    # Cutoff zero is really infinite
    if float(mdp_opts['rlist']) == 0.0:
        mdin.change('cntrl','cut','9999')
    else:
        mdin.change('cntrl','cut',str(int(float(mdp_opts['rlist'])*10)))
    # Take zero MD steps
    mdin.change('cntrl','nstlim','0')
    # Don't update nonbond parameters
    mdin.change('cntrl','nsnb','0')
    # if mdp_opts['coulombtype'].lower() == 'pme':
    #     mdin.change('ewald','order',5)
    #     mdin.change('ewald','skinnb',0)
    mdin.write("mdin")
    # Nonbonded method
    if mdp_opts['coulombtype'].lower() == 'pme':
        with open("mdin",'a') as f:
            print >> f, """&ewald
 order=5, skinnb=0
/"""
    with open("mdin",'a') as f:
        print >> f, """&debugf
do_debugf=1, dumpfrc=1
/"""
    # Call sander for energy and force
    os.system('rm -f forcedump.dat')
    _exec("sander -O -y mdcrd", print_command=False)
    # Parse energy and force
    ParseMode = 0
    Energies = []
    Forces = []
    Force = []
    iatom = 0
    isAtom = [atom.atomic_number > 0 for atom in Structure.atoms]
    for line in open('forcedump.dat'):
        line = line.strip()
        sline = line.split()
        if ParseMode == 1:
            if len(sline) == 1 and isfloat(sline[0]):
                Energies.append(float(sline[0]) * 4.184)
                ParseMode = 0
        if ParseMode == 2:
            if len(sline) == 3 and all(isfloat(sline[i]) for i in range(3)):
                if isAtom[iatom]:
                    Force += [float(sline[i]) * 4.184 * 10 for i in range(3)]
                iatom += 1
            if len(Force) == 3*sum(isAtom):
                Forces.append(np.array(Force))
                Force = []
                ParseMode = 0
                iatom = 0
        if line == '0 START of Energies':
            ParseMode = 1
        elif line == '1 Total Force' or line == '2 Total Force':
            ParseMode = 2
    # Obtain energy components
    ParseMode = 0
    Ecomps = OrderedDict()
    for line in open("mdout").readlines():
        if "NSTEP = " in line:
            ParseMode = 1
        if ParseMode == 1:
            if "=" not in line:
                ParseMode = 0
                continue
            else:
                ieq = None
                wkey = []
                # Assume the line is split-able
                for i, w in enumerate(line.split()):
                    if w == '=':
                        ieq = i
                    elif i-1 == ieq:
                        Ecomps.setdefault(' '.join(wkey), []).append(float(w)*4.184)
                        wkey = []
                    else:
                        wkey.append(w)
    Ecomps_Sav = OrderedDict()
    for key in Ecomps:
        if set(Ecomps[key]) == set([0.0]): continue
        elif key.lower() in ['eptot', 'etot', 'volume', 'density']: continue
        else:
            Ecomps_Sav[key] = Ecomps[key][0]
    Ecomps_Sav['EPTOT'] = Ecomps['EPtot'][0]
    # Save just the first frame from the .mdcrd
    Energies = Energies[0]
    Forces = Forces[0]
    return Energies, Forces, Ecomps_Sav

def main():
    # Command line arguments
    gro_file = sys.argv[1]
    top_file = sys.argv[2]
    mdp_file = sys.argv[3]

    # Parse the .mdp file to inform ParmEd
    defines, sysargs, mdp_opts = interpret_mdp(mdp_file)

    # Gromacs calculation
    GMX_Energy, GMX_Force, Ecomps_GMX = Calculate_GMX(gro_file, top_file, mdp_file)
    GMX_Force = GMX_Force.reshape(-1,3)

    # Print Gromacs energy components
    printcool_dictionary(Ecomps_GMX, title="GROMACS energy components")

    # ParmEd-OpenMM calculation
    Structure, OMM_Energy, OMM_Force, Ecomps_OMM = Calculate_ParmEd_OpenMM(gro_file, top_file, sysargs, defines)
    
    if args.amber:
        # AMBER calculation (optional)
        AMBER_Energy, AMBER_Force, Ecomps_AMBER = Calculate_AMBER(Structure, mdp_opts)
        AMBER_Force = AMBER_Force.reshape(-1,3)
        # Print AMBER energy components
        printcool_dictionary(Ecomps_AMBER, title="AMBER energy components")

    # Construct arrays of energy and force differences
    if args.amber:
        Names = ['Gromacs', 'OpenMM', 'AMBER']
        Energies = np.array([GMX_Energy, OMM_Energy.value_in_unit(u.kilojoule_per_mole), AMBER_Energy])
        Forces = np.array([GMX_Force, OMM_Force, AMBER_Force])
    else:
        Names = ['Gromacs', 'OpenMM']
        Energies = np.array([GMX_Energy, OMM_Energy.value_in_unit(u.kilojoule_per_mole)])
        Forces = np.array([GMX_Force, OMM_Force])
    D_Energy = []
    D_FrcRMS = []
    D_FrcMax = []
    D_Names = []
    for i in range(1, len(Names)):
        for j in range(i):
            D_Names.append('%s-%s' % (Names[j],Names[i]))
            D_Energy.append(Energies[j]-Energies[i])
            D_Force = Forces[j]-Forces[i]
            D_FrcRMS.append(np.sqrt(np.mean([sum(k**2) for k in D_Force])))
            D_FrcMax.append(np.sqrt(np.max(np.array([sum(k**2) for k in D_Force]))))

    # Print the net force on the first three atoms (e.g. water molecule)
    # print np.sum(GMX_Force[:3], axis=0)
    # print np.sum(AMBER_Force[:3], axis=0)
    # Final printout
    print "Energy Difference (kJ/mol):"
    for i in range(len(D_Names)):
        print "%-14s % .6e" % (D_Names[i], D_Energy[i])

    print "RMS / Max Force Difference (kJ/mol/nm):"
    for i in range(len(D_Names)):
        print "%-14s % .6e % .6e" % (D_Names[i], D_FrcRMS[i], D_FrcMax[i])

if __name__ == "__main__":
    main()
