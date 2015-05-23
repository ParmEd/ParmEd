#!/usr/bin/env python

import numpy as np
import networkx as nx
import simtk.unit as u
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser(description='Calculate LJ and Coulomb energy for DPPC-cut system')
parser.add_argument('--fudgeLJ', type=float, default=0.0, help='1-4 scale factor for LJ interactions')
parser.add_argument('--fudgeQQ', type=float, default=0.0, help='1-4 scale factor for Coulomb interactions')
parser.add_argument('--nrexcl', type=int, default=3, help='Exclude ')
parser.add_argument('--top', type=str, help='Gromacs .top file (for setting up GMX)')
parser.add_argument('--mdp', type=str, help='Gromacs .mdp file (for setting up GMX)')
parser.add_argument('--pdb', type=str, help='pdb (for setting up OpenMM simulation)')
args = parser.parse_args()

# Coordinates are taken from DPPC_cut.gro
# Force field parameters are taken from DPPC_cut.itp

# Coordinates in Angstroms
xyzs = np.array([[6.009832899999999256e+00, 7.240071079999999881e+00, 1.291002618999999996e+01],
                 [7.760541820000000257e+00, 6.950370019999999371e+00, 1.131817631999999918e+01],
                 [8.069597569999999109e+00, 7.439569630000000267e+00, 1.372942429999999980e+01],
                 [7.349989380000000239e+00, 6.672915480000000343e+00, 1.270214368999999976e+01],
                 [7.371385290000000090e+00, 5.216477320000000084e+00, 1.290020962000000004e+01],
                 [8.758457280000000011e+00, 4.602699010000000257e+00, 1.269969157000000060e+01],
                 [8.891755039999999610e+00, 3.217678310000000153e+00, 1.302961484999999975e+01],
                 [1.027124525999999882e+01, 2.801014829999999733e+00, 1.231163276000000018e+01],
                 [1.039975174000000102e+01, 1.402673599999999965e+00, 1.277910243000000001e+01],
                 [1.029997322999999909e+01, 3.039785810000000144e+00, 1.085130295999999994e+01],
                 [1.125888819000000041e+01, 3.843151119999999921e+00, 1.304006489999999907e+01],
                 [1.203699261999999948e+01, 2.983288120000000099e+00, 1.387677818000000052e+01],
                 [1.313805257000000104e+01, 3.812129639999999764e+00, 1.454129959999999855e+01],
                 [1.413727259999999930e+01, 3.140010959999999685e+00, 1.532171026000000147e+01],
                 [1.535222652000000032e+01, 3.061570319999999956e+00, 1.471562937999999932e+01],
                 [1.559981773999999888e+01, 3.769482379999999910e+00, 1.374071446000000130e+01],
                 [1.227992979999999967e+01, 4.787202320000000455e+00, 1.534985298000000142e+01],
                 [1.129224264999999860e+01, 4.130405409999999833e+00, 1.614859694000000090e+01],
                 [1.017820885999999980e+01, 4.821872579999999964e+00, 1.650970891000000051e+01],
                 [9.810362630000000195e+00, 5.819019909999999740e+00, 1.589060766000000058e+01]])

# Atom names, probably not needed.
atomnames = ['C1', 'C2', 'C3', 'N4', 'C5', 'C6', 'O7', 'P8', 'O9', 'O10', 'O11', 'C12', 'C13', 'O14', 'C15', 'O16', 'C17', 'O18', 'C19', 'O20']

# Dictionary mapping atom names to atom types.
atomnames_to_atomtypes = {'C1':'C3', 'C2':'C3', 'C3':'C3', 'N4':'NL', 'C5':'H2', 'C6':'C2', 'O7':'OS', 'P8':'P' , 'O9':'OM','O10':'OM',
                          'O11':'OS','C12':'C2','C13':'H1','O14':'OS','C15': 'C','O16': 'O','C17':'C2','O18':'OS','C19': 'C','O20': 'O'}

# Permittivity of free space, used for calculating Coulomb interactions.
vacuum_permittivity = 8.854187817620e-12 * u.coulomb ** 2 / u.newton / u.meter ** 2

# Dictionary mapping atom types to Lennard-Jones parameters.
atomtypes_to_lj = {'O':(2.96000e-01, 8.78694e-01),
                   'OM':(2.96000e-01, 8.78694e-01),
                   'NL':(3.24997e-01, 7.11377e-01),
                   'C':(3.74986e-01, 4.39524e-01),
                   'H1':(3.80015e-01, 3.34616e-01),
                   'H2':(3.90504e-01, 4.93637e-01),
                   'P':(3.74004e-01, 8.36713e-01),
                   'OS':(2.99989e-01, 8.79145e-01),
                   'C3':(3.96012e-01, 6.06494e-01),
                   'C2':(3.79979e-01, 4.93950e-01),
                   'OW':(3.16565e-01, 6.50167e-01),
                   'H':(0.00000e+00, 0.00000e+00)}

# Dictionary mapping atom names to charge parameters.
atomnames_to_charge = {'C1':0.4000, 'C2':0.4000, 'C3':0.4000, 'N4':-0.5000, 'C5':0.3000, 'C6':0.4000, 'O7':-0.8000, 'P8':1.7000, 'O9':-0.8000, 'O10':-0.8000, 
                       'O11':-0.7000, 'C12':0.4000, 'C13':0.3000, 'O14':-0.7000, 'C15':0.7000, 'O16':-0.7000, 'C17':0.5000, 'O18':-0.7000, 'C19':0.8000, 'O20':-0.6000}

# Set of dihedral atomtypes
dihedral_atomtypes = set()

# Bonds are indexed from one
bonds=[(1,4), (2, 4), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (8, 10), (8, 11), (11, 12), (12, 13), (13, 14), (13, 17), (14, 15), (15, 16), (17, 18), (18, 19), (19, 20)]

# Dihedrals are indexed from one
dihedrals=[(1, 4, 5, 6), (4, 5, 6, 7), (5, 6, 7, 8), (6, 7, 8, 11), (6, 7, 8, 11), (7, 8, 11, 12), (7, 8, 11, 12), (8, 11, 12, 13), (11, 12, 13, 14), 
           (11, 12, 13, 17), (11, 12, 13, 17), (12, 13, 17, 18), (12, 13, 17, 18), (12, 13, 14, 15), (13, 17, 18, 19), (13, 17, 18, 19)]


G = nx.Graph()
for i in range(len(atomnames)):
    G.add_node(i)
for i in bonds:
    G.add_edge(i[0]-1, i[1]-1)

# Print out missing torsions if desired
print_torsions = 0
if print_torsions:
    print "Finding missing torsions"
    for e in G.edges():
        for n0 in G.neighbors(e[0]):
            for n1 in G.neighbors(e[1]):
                if n0 in [e[0], e[1]]: continue
                if n1 in [e[0], e[1]]: continue
                dihfwd = (n0+1, e[0]+1, e[1]+1, n1+1)
                dihbak = (n1+1, e[1]+1, e[0]+1, n0+1)
                i, j, k, l = dihfwd
                if dihfwd not in dihedrals and dihbak not in dihedrals:
                    print "%2i    %2i    %2i    %2i    1    0.0     0.0     0" % (i, j, k, l)
                dihfwdat = tuple(list((atomnames_to_atomtypes[atomnames[i-1]] for i in dihfwd)))
                dihbakat = tuple(list((atomnames_to_atomtypes[atomnames[i-1]] for i in dihbak)))
                if dihbakat not in dihedral_atomtypes and dihfwdat not in dihedral_atomtypes:
                    dihedral_atomtypes.add(dihfwdat)
    
# for dt in dihedral_atomtypes:
#     i, j, k, l = dt
#     print "%2s    %2s    %2s    %2s    1    0.0     0.0     0" % (i, j, k, l)
            
# raw_input()

# Calculate the Coulomb and Lennard-Jones energy "by hand."
coulomb_energy = 0.0
lj_energy = 0.0
fudgeLJ = args.fudgeLJ
fudgeQQ = args.fudgeQQ
nrexcl = args.nrexcl
excldict = OrderedDict()
for i in range(len(atomnames)):
    for j in range(i):
        if len(nx.shortest_path(G, i, j)) == 4:
            lj_prefactor = fudgeLJ
            qq_prefactor = fudgeQQ
        elif len(nx.shortest_path(G, i, j)) <= (nrexcl+1):
            lj_prefactor = 0.0
            qq_prefactor = 0.0
        else:
            lj_prefactor = 1.0
            qq_prefactor = 1.0
        # Gromacs treats a 1-4 LJ+electrostatics interaction more like a "bonded interaction"
        # so they are technically excluded from the list of "nonbonded interactions"
        if len(nx.shortest_path(G, i, j)) <= (nrexcl+1):
            excldict.setdefault(j, []).append(i)
            excldict.setdefault(i, []).append(j)
        rij = np.linalg.norm(xyzs[i]-xyzs[j])/10
        ani = atomnames[i]
        anj = atomnames[j]
        ati = atomnames_to_atomtypes[ani]
        atj = atomnames_to_atomtypes[anj]
        qi = atomnames_to_charge[ani]
        qj = atomnames_to_charge[anj]
        sigi = atomtypes_to_lj[ati][0]
        sigj = atomtypes_to_lj[atj][0]
        sigij = (sigi+sigj)*0.5
        epsi = atomtypes_to_lj[ati][1]
        epsj = atomtypes_to_lj[atj][1]
        epsij = (epsi*epsj)**0.5
        sig_div_r = sigij/rij
        coulomb_pair = qq_prefactor*qi*qj*u.elementary_charge**2 / (rij*u.nanometers * 4 * np.pi * vacuum_permittivity) * u.AVOGADRO_CONSTANT_NA
        coulomb_energy += coulomb_pair.value_in_unit(u.kilojoule_per_mole)
        lj_pair = lj_prefactor * 4*epsij*(sig_div_r**12 - sig_div_r**6)
        lj_energy += lj_pair

# Gromacs exclusion printout excludes itself.
for i in range(len(atomnames)):
    excldict.setdefault(i, []).append(i)

for k, v in excldict.items():
    excldict[k] = sorted(list(set(v)))

nra = 0
for k, v in excldict.items():
    print "      excls[%i][%i..%i]={%s}" % (k, nra, nra+len(v)-1, ', '.join(['%i' % i for i in v]))
    nra += len(v)
            
print "The Coulomb energy is", coulomb_energy, "kJ/mol"
print "The Lennard-Jones energy is", lj_energy, "kJ/mol"
