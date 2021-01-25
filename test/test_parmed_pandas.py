"""
Tests the pandas functionality of the parmed/structure module
"""
from __future__ import division, print_function, absolute_import
from utils import get_fn

import parmed.structure as structure
from parmed import load_file
from parmed.topologyobjects import Atom
from parmed.utils.six.moves import zip
import numpy as np
import os
import unittest
from utils import create_random_structure
try:
    import pandas as pd
except ImportError:
    pd = None

@unittest.skipIf(pd is None, "Cannot test without pandas")
class TestStructureDataFrame(unittest.TestCase):
    """ Tests the conversion of parmed.Structure to a pd.DataFrame """

    def testPDBPandas(self):
        """ Test 4lzt.pdb conversion to DataFrame (w/ anisou and positions) """
        pdb = load_file(get_fn('4lzt.pdb'))
        df = pdb.to_dataframe()
        self.assertIn('U11', df)
        self.assertIn('U22', df)
        self.assertIn('U33', df)
        self.assertIn('U12', df)
        self.assertIn('U13', df)
        self.assertIn('U23', df)
        # Make sure that when one of the anisotropic B-factors is not set, the
        # NaN gets properly put in its place
        pdb.atoms[0].anisou = None
        df2 = pdb.to_dataframe()
        self.assertIn('U11', df)
        self.assertIn('U22', df)
        self.assertIn('U33', df)
        self.assertIn('U12', df)
        self.assertIn('U13', df)
        self.assertIn('U23', df)
        # Now make sure they are all equal except for the first atom, and that
        # the data matches the PDB information
        anisou = ['U11', 'U22', 'U33', 'U12', 'U13', 'U23']
        for (i, r1), (j, r2) in zip(df.iterrows(), df2.iterrows()):
            if i == 0:
                # Easy NaN check
                self.assertNotEqual(r2.U11, r2.U11)
                self.assertNotEqual(r2.U22, r2.U22)
                self.assertNotEqual(r2.U33, r2.U33)
                self.assertNotEqual(r2.U12, r2.U12)
                self.assertNotEqual(r2.U13, r2.U13)
                self.assertNotEqual(r2.U23, r2.U23)
                # Check against df1
                self.assertNotEqual(r1.U11, r2.U11)
                self.assertNotEqual(r1.U22, r2.U22)
                self.assertNotEqual(r1.U33, r2.U33)
                self.assertNotEqual(r1.U12, r2.U12)
                self.assertNotEqual(r1.U13, r2.U13)
                self.assertNotEqual(r1.U23, r2.U23)
            else:
                self.assertEqual(r1.U11, r2.U11)
                self.assertEqual(r1.U22, r2.U22)
                self.assertEqual(r1.U33, r2.U33)
                self.assertEqual(r1.U12, r2.U12)
                self.assertEqual(r1.U13, r2.U13)
                self.assertEqual(r1.U23, r2.U23)
                self.assertTrue(np.all(np.abs(r1[anisou] - pdb.atoms[i].anisou) < 1e-4))
            self.assertEqual(r1.charge, 0)
            self.assertEqual(r1.charge, pdb.atoms[i].charge)
            self.assertEqual(r1['name'], pdb.atoms[i].name)
            self.assertEqual(r1.solvent_radius, pdb.atoms[i].solvent_radius)
            self.assertEqual(r1.screen, pdb.atoms[i].screen)
            self.assertEqual(r1.type, pdb.atoms[i].type)
            self.assertEqual(r1.occupancy, pdb.atoms[i].occupancy)
            self.assertEqual(r1.bfactor, pdb.atoms[i].bfactor)
            self.assertEqual(r1.altloc, pdb.atoms[i].altloc)
            self.assertEqual(r1.rmin, pdb.atoms[i].rmin)
            self.assertEqual(r1.epsilon, pdb.atoms[i].epsilon)
            self.assertEqual(r1.rmin_14, pdb.atoms[i].rmin_14)
            self.assertEqual(r1.epsilon_14, pdb.atoms[i].epsilon_14)
            self.assertEqual(r1.resname, pdb.atoms[i].residue.name)
            self.assertEqual(r1.resid, pdb.atoms[i].residue.idx)
            self.assertEqual(r1.resnum, pdb.atoms[i].residue.number)
            self.assertEqual(r1.chain, pdb.atoms[i].residue.chain)
            self.assertEqual(r1.join, pdb.atoms[i].join)
            self.assertEqual(r1.nb_idx, pdb.atoms[i].nb_idx)
            self.assertEqual(r1.xx, pdb.atoms[i].xx)
            self.assertEqual(r1.xy, pdb.atoms[i].xy)
            self.assertEqual(r1.xz, pdb.atoms[i].xz)

        self.assertNotIn('vx', df)
        self.assertNotIn('vy', df)
        self.assertNotIn('vz', df)
        self.assertNotIn('type_idx', df)
        self.assertNotIn('class_idx', df)
        for key in df.keys():
            self.assertFalse(key.startswith('multipole'))
        self.assertNotIn('polarizability', df)
        self.assertNotIn('vdw_parent', df)
        self.assertIn('segid', df)

    def testStructureViewPandas(self):
        """ Tests creating a pandas DataFrame from a StructureView """
        parm = load_file(get_fn('tip4p.parm7'))
        parm.load_rst7(get_fn('tip4p.rst7'))
        df = parm.view[:10,:].to_dataframe()
        self.assertEqual(df.shape[0], sum(len(res) for res in parm.residues[:10]))
        for i, r1 in df.iterrows():
            self.assertEqual(r1.charge, parm.atoms[i].charge)
            self.assertEqual(r1['name'], parm.atoms[i].name)
            self.assertEqual(r1.solvent_radius, parm.atoms[i].solvent_radius)
            self.assertEqual(r1.screen, parm.atoms[i].screen)
            self.assertEqual(r1.type, parm.atoms[i].type)
            self.assertEqual(r1.occupancy, parm.atoms[i].occupancy)
            self.assertEqual(r1.bfactor, parm.atoms[i].bfactor)
            self.assertEqual(r1.altloc, parm.atoms[i].altloc)
            self.assertEqual(r1.rmin, parm.atoms[i].rmin)
            self.assertEqual(r1.epsilon, parm.atoms[i].epsilon)
            self.assertEqual(r1.rmin_14, parm.atoms[i].rmin_14)
            self.assertEqual(r1.epsilon_14, parm.atoms[i].epsilon_14)
            self.assertEqual(r1.resname, parm.atoms[i].residue.name)
            self.assertEqual(r1.resid, parm.atoms[i].residue.idx)
            self.assertEqual(r1.resnum, parm.atoms[i].residue.number)
            self.assertEqual(r1.chain, parm.atoms[i].residue.chain)
            self.assertEqual(r1.join, parm.atoms[i].join)
            self.assertEqual(r1.nb_idx, parm.atoms[i].nb_idx)
            self.assertEqual(r1.xx, parm.atoms[i].xx)
            self.assertEqual(r1.xy, parm.atoms[i].xy)
            self.assertEqual(r1.xz, parm.atoms[i].xz)
            self.assertEqual(r1.vx, parm.atoms[i].vx)
            self.assertEqual(r1.vy, parm.atoms[i].vy)
            self.assertEqual(r1.vz, parm.atoms[i].vz)

        self.assertNotIn('type_idx', df)
        self.assertNotIn('class_idx', df)
        for key in df.keys():
            self.assertFalse(key.startswith('multipole'))
        self.assertNotIn('polarizability', df)
        self.assertNotIn('vdw_parent', df)
        self.assertIn('segid', df)

    def testAmberParmPandas(self):
        """ Tests creating a pandas DataFrame from an AmberParm """
        parm = load_file(get_fn('tip4p.parm7'))
        parm.load_rst7(get_fn('tip4p.rst7'))
        df = parm.to_dataframe()
        self.assertEqual(df.shape[0], len(parm.atoms))
        for i, r1 in df.iterrows():
            self.assertEqual(r1.charge, parm.atoms[i].charge)
            self.assertEqual(r1['name'], parm.atoms[i].name)
            self.assertEqual(r1.solvent_radius, parm.atoms[i].solvent_radius)
            self.assertEqual(r1.screen, parm.atoms[i].screen)
            self.assertEqual(r1.type, parm.atoms[i].type)
            self.assertEqual(r1.occupancy, parm.atoms[i].occupancy)
            self.assertEqual(r1.bfactor, parm.atoms[i].bfactor)
            self.assertEqual(r1.altloc, parm.atoms[i].altloc)
            self.assertEqual(r1.rmin, parm.atoms[i].rmin)
            self.assertEqual(r1.epsilon, parm.atoms[i].epsilon)
            self.assertEqual(r1.rmin_14, parm.atoms[i].rmin_14)
            self.assertEqual(r1.epsilon_14, parm.atoms[i].epsilon_14)
            self.assertEqual(r1.resname, parm.atoms[i].residue.name)
            self.assertEqual(r1.resid, parm.atoms[i].residue.idx)
            self.assertEqual(r1.resnum, parm.atoms[i].residue.number)
            self.assertEqual(r1.chain, parm.atoms[i].residue.chain)
            self.assertEqual(r1.join, parm.atoms[i].join)
            self.assertEqual(r1.nb_idx, parm.atoms[i].nb_idx)
            self.assertEqual(r1.xx, parm.atoms[i].xx)
            self.assertEqual(r1.xy, parm.atoms[i].xy)
            self.assertEqual(r1.xz, parm.atoms[i].xz)
            self.assertEqual(r1.vx, parm.atoms[i].vx)
            self.assertEqual(r1.vy, parm.atoms[i].vy)
            self.assertEqual(r1.vz, parm.atoms[i].vz)

        self.assertNotIn('type_idx', df)
        self.assertNotIn('class_idx', df)
        for key in df.keys():
            self.assertFalse(key.startswith('multipole'))
        self.assertNotIn('polarizability', df)
        self.assertNotIn('vdw_parent', df)
        self.assertIn('segid', df)

    def testAmoebaParmPandas(self):
        """ Tests creating a pandas DataFrame from an AmoebaParm """
        parm = load_file(get_fn('nma.parm7'))
        df = parm.to_dataframe()
        self.assertEqual(df.shape[0], parm.ptr('natom'))
        multipoles = ['multipole_111', 'multipole_211', 'multipole_212',
                      'multipole_222', 'multipole_411', 'multipole_412',
                      'multipole_422', 'multipole_413', 'multipole_423',
                      'multipole_433']
        for i, r1 in df.iterrows():
            self.assertEqual(r1.charge, parm.atoms[i].charge)
            self.assertEqual(r1['name'], parm.atoms[i].name)
            self.assertEqual(r1.solvent_radius, parm.atoms[i].solvent_radius)
            self.assertEqual(r1.screen, parm.atoms[i].screen)
            self.assertEqual(r1.type, parm.atoms[i].type)
            self.assertEqual(r1.occupancy, parm.atoms[i].occupancy)
            self.assertEqual(r1.bfactor, parm.atoms[i].bfactor)
            self.assertEqual(r1.altloc, parm.atoms[i].altloc)
            self.assertEqual(r1.rmin, parm.atoms[i].rmin)
            self.assertEqual(r1.epsilon, parm.atoms[i].epsilon)
            self.assertEqual(r1.rmin_14, parm.atoms[i].rmin_14)
            self.assertEqual(r1.epsilon_14, parm.atoms[i].epsilon_14)
            self.assertEqual(r1.resname, parm.atoms[i].residue.name)
            self.assertEqual(r1.resid, parm.atoms[i].residue.idx)
            self.assertEqual(r1.resnum, parm.atoms[i].residue.number)
            self.assertEqual(r1.chain, parm.atoms[i].residue.chain)
            self.assertEqual(r1.join, parm.atoms[i].join)
            self.assertEqual(r1.nb_idx, parm.atoms[i].nb_idx)
            self.assertEqual(r1.type_idx, parm.atoms[i].type_idx)
            self.assertEqual(r1.class_idx, parm.atoms[i].class_idx)
            self.assertTrue(np.all(r1[multipoles] ==
                np.array(parm.atoms[i].multipoles)))
            self.assertEqual(r1.polarizability, parm.atoms[i].polarizability)
            self.assertEqual(r1.vdw_parent, parm.atoms[i].vdw_parent.idx)

        self.assertNotIn('xx', df)
        self.assertNotIn('xy', df)
        self.assertNotIn('xz', df)
        self.assertNotIn('vx', df)
        self.assertNotIn('vy', df)
        self.assertNotIn('vz', df)
        self.assertIn('segid', df)

    def testCharmmPSFPandas(self):
        """ Tests creating a pandas DataFrame from a CharmmPsfFile """
        parm = load_file(get_fn('ala_ala_ala.psf'))
        df = parm.to_dataframe()
        self.assertEqual(df.shape[0], len(parm.atoms))
        for i, r1 in df.iterrows():
            self.assertEqual(r1.charge, parm.atoms[i].charge)
            self.assertEqual(r1['name'], parm.atoms[i].name)
            self.assertEqual(r1.solvent_radius, parm.atoms[i].solvent_radius)
            self.assertEqual(r1.screen, parm.atoms[i].screen)
            self.assertEqual(r1.type, parm.atoms[i].type)
            self.assertEqual(r1.occupancy, parm.atoms[i].occupancy)
            self.assertEqual(r1.bfactor, parm.atoms[i].bfactor)
            self.assertEqual(r1.altloc, parm.atoms[i].altloc)
            self.assertEqual(r1.rmin, parm.atoms[i].rmin)
            self.assertEqual(r1.epsilon, parm.atoms[i].epsilon)
            self.assertEqual(r1.rmin_14, parm.atoms[i].rmin_14)
            self.assertEqual(r1.epsilon_14, parm.atoms[i].epsilon_14)
            self.assertEqual(r1.resname, parm.atoms[i].residue.name)
            self.assertEqual(r1.resid, parm.atoms[i].residue.idx)
            self.assertEqual(r1.resnum, parm.atoms[i].residue.number)
            self.assertEqual(r1.chain, parm.atoms[i].residue.chain)
            self.assertEqual(r1.join, parm.atoms[i].join)
            self.assertEqual(r1.nb_idx, parm.atoms[i].nb_idx)
            self.assertEqual(r1.segid, parm.atoms[i].residue.segid)

        self.assertNotIn('xx', df)
        self.assertNotIn('xy', df)
        self.assertNotIn('xz', df)
        self.assertNotIn('vx', df)
        self.assertNotIn('vy', df)
        self.assertNotIn('vz', df)
        self.assertNotIn('type_idx', df)
        self.assertNotIn('class_idx', df)
        for key in df.keys():
            self.assertFalse(key.startswith('multipole'))
        self.assertNotIn('polarizability', df)
        self.assertNotIn('vdw_parent', df)

    def testLoadDataFrameStructure(self):
        """ Tests the load_dataframe method on Structure """
        struct = create_random_structure(parametrized=True)
        charges = [a.charge for a in struct.atoms]
        self.assertTrue(not all(x == 0 for x in charges))
        df = struct.to_dataframe()
        # First zero-out all of the charges
        struct.load_dataframe(dict(charge=[0 for a in struct.atoms]))
        self.assertTrue(all(a.charge == 0 for a in struct.atoms))
        # Now re-load the dataframe to restore the original charges
        struct.load_dataframe(df)
        self.assertTrue(all(a.charge == x for a, x in zip(struct.atoms, charges)))
        # Change the first atomic properties of *everything* now to
        # make sure that they all get updated
        df_orig = df.copy()
        df.loc[0, 'number'] = 1
        df.loc[0, 'name'] = 'HAHA'
        df.loc[0, 'type'] = 'FUNY'
        df.loc[0, 'atomic_number'] = 92 # uranium
        df.loc[0, 'charge'] *= 2
        df.loc[0, 'mass'] *= 10
        df.loc[0, 'nb_idx'] = 10
        df.loc[0, 'solvent_radius'] *= 2
        df.loc[0, 'screen'] = 0.5
        df.loc[0, 'occupancy'] = 0.1
        df.loc[0, 'bfactor'] = 0.5
        df.loc[0, 'altloc'] = 'X'
        df.loc[0, 'tree'] = 'BLU'
        df.loc[0, 'join'] = 1.0
        df.loc[0, 'irotat'] = 2.0
        df.loc[0, 'rmin'] *= 2
        df.loc[0, 'epsilon'] /= 2
        struct.load_dataframe(df)
        atom = struct.atoms[0]
        self.assertEqual(atom.number, 1)
        self.assertEqual(atom.name, 'HAHA')
        self.assertEqual(atom.type, 'FUNY')
        self.assertEqual(atom.atomic_number, 92)
        self.assertEqual(atom.charge, charges[0]*2)
        self.assertEqual(atom.mass, 10*df_orig.loc[0, 'mass'])
        self.assertEqual(atom.nb_idx, 10)
        self.assertEqual(atom.solvent_radius, 2*df_orig.loc[0, 'solvent_radius'])
        self.assertEqual(atom.screen, 0.5)
        self.assertEqual(atom.occupancy, 0.1)
        self.assertEqual(atom.bfactor, 0.5)
        self.assertEqual(atom.altloc, 'X')
        self.assertEqual(atom.tree, 'BLU')
        self.assertEqual(atom.join, 1.0)
        self.assertEqual(atom.irotat, 2.0)
        self.assertEqual(atom.rmin, 2*df_orig.loc[0, 'rmin'])
        self.assertEqual(atom.epsilon, df_orig.loc[0, 'epsilon']/2)

    def testLoadDataFrameStructureView(self):
        """ Tests the load_dataframe method on StructureView """
        struct = create_random_structure(parametrized=True).view[:10,:]
        charges = [a.charge for a in struct.atoms]
        self.assertTrue(not all(x == 0 for x in charges))
        df = struct.to_dataframe()
        # First zero-out all of the charges
        struct.load_dataframe(dict(charge=[0 for a in struct.atoms]))
        self.assertTrue(all(a.charge == 0 for a in struct.atoms))
        # Now re-load the dataframe to restore the original charges
        struct.load_dataframe(df)
        self.assertTrue(all(a.charge == x for a, x in zip(struct.atoms, charges)))
        # Change the first atomic properties of *everything* now to
        # make sure that they all get updated
        df_orig = df.copy()
        df.loc[0, 'number'] = 1
        df.loc[0, 'name'] = 'HAHA'
        df.loc[0, 'type'] = 'FUNY'
        df.loc[0, 'atomic_number'] = 92 # uranium
        df.loc[0, 'charge'] *= 2
        df.loc[0, 'mass'] *= 10
        df.loc[0, 'nb_idx'] = 10
        df.loc[0, 'solvent_radius'] *= 2
        df.loc[0, 'screen'] = 0.5
        df.loc[0, 'occupancy'] = 0.1
        df.loc[0, 'bfactor'] = 0.5
        df.loc[0, 'altloc'] = 'X'
        df.loc[0, 'tree'] = 'BLU'
        df.loc[0, 'join'] = 1.0
        df.loc[0, 'irotat'] = 2.0
        df.loc[0, 'rmin'] *= 2
        df.loc[0, 'epsilon'] /= 2
        struct.load_dataframe(df)
        atom = struct.atoms[0]
        self.assertEqual(atom.number, 1)
        self.assertEqual(atom.name, 'HAHA')
        self.assertEqual(atom.type, 'FUNY')
        self.assertEqual(atom.atomic_number, 92)
        self.assertEqual(atom.charge, charges[0]*2)
        self.assertEqual(atom.mass, 10*df_orig.loc[0, 'mass'])
        self.assertEqual(atom.nb_idx, 10)
        self.assertEqual(atom.solvent_radius, 2*df_orig.loc[0, 'solvent_radius'])
        self.assertEqual(atom.screen, 0.5)
        self.assertEqual(atom.occupancy, 0.1)
        self.assertEqual(atom.bfactor, 0.5)
        self.assertEqual(atom.altloc, 'X')
        self.assertEqual(atom.tree, 'BLU')
        self.assertEqual(atom.join, 1.0)
        self.assertEqual(atom.irotat, 2.0)
        self.assertEqual(atom.rmin, 2*df_orig.loc[0, 'rmin'])
        self.assertEqual(atom.epsilon, df_orig.loc[0, 'epsilon']/2)
