#!/usr/bin/env python

import numpy as np
import networkx as nx
import simtk.unit as u

xyzs = np.loadtxt('DPPC_xyz.txt')

atomnames = ['C1', 'C2', 'C3', 'N4', 'C5', 'C6', 'O7', 'P8', 'O9', 'O10', 'O11', 'C12', 'C13', 'O14', 'C15', 'O16', 'C17', 'O18', 'C19', 'O20']

atomnames_to_atomtypes = {'C1':'C3', 'C2':'C3', 'C3':'C3', 'N4':'NL', 'C5':'H2', 'C6':'C2', 'O7':'OS', 'P8':'P' , 'O9':'OM','O10':'OM',
                          'O11':'OS','C12':'C2','C13':'H1','O14':'OS','C15': 'C','O16': 'O','C17':'C2','O18':'OS','C19': 'C','O20': 'O'}

vacuum_permittivity = 8.854187817620e-12 * u.coulomb ** 2 / u.newton / u.meter ** 2

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

atomnames_to_charge = {'C1':0.4000, 'C2':0.4000, 'C3':0.4000, 'N4':-0.5000, 'C5':0.3000, 'C6':0.4000, 'O7':-0.8000, 'P8':1.7000, 'O9':-0.8000, 'O10':-0.8000, 
                       'O11':-0.7000, 'C12':0.4000, 'C13':0.3000, 'O14':-0.7000, 'C15':0.7000, 'O16':-0.7000, 'C17':0.5000, 'O18':-0.7000, 'C19':0.8000, 'O20':-0.6000}

# Bonds are indexed from one
bonds=[(1,4), (2, 4), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (8, 10), (8, 11), (11, 12), (12, 13), (13, 14), (13, 17), (14, 15), (15, 16), (17, 18), (18, 19), (19, 20)]

G = nx.Graph()
for i in range(len(atomnames)):
    G.add_node(i)
for i in bonds:
    G.add_edge(i[0]-1, i[1]-1)

# Calculate the Coulomb and Lennard-Jones energy by hand. :P
exclusions = []
coulomb_energy = 0.0
lj_energy = 0.0
for i in range(len(atomnames)):
    for j in range(i):
        if len(nx.shortest_path(G, i, j)) <= 4:
            exclusions.append((i, j))
        else:
            rij = np.linalg.norm(xyzs[i]-xyzs[j])/10
            ani = atomnames[i]
            anj = atomnames[j]
            ati = atomnames_to_atomtypes[ani]
            atj = atomnames_to_atomtypes[anj]
            qi = atomnames_to_charge[ani]
            qj = atomnames_to_charge[anj]
            coulomb_pair = qi*qj*u.elementary_charge**2 / (rij*u.nanometers * 4 * np.pi * vacuum_permittivity) * u.AVOGADRO_CONSTANT_NA
            coulomb_energy += coulomb_pair.value_in_unit(u.kilojoule_per_mole)
            sigi = atomtypes_to_lj[ati][0]
            sigj = atomtypes_to_lj[atj][0]
            sigij = (sigi+sigj)*0.5
            epsi = atomtypes_to_lj[ati][1]
            epsj = atomtypes_to_lj[atj][1]
            epsij = (epsi*epsj)**0.5
            sig_div_r = sigij/rij
            lj_pair = 4*epsij*(sig_div_r**12 - sig_div_r**6)
            lj_energy += lj_pair
            
print "The hand-calculated Coulomb energy is", coulomb_energy, "kJ/mol"
print "The hand-calculated Lennard-Jones energy is", lj_energy, "kJ/mol"
