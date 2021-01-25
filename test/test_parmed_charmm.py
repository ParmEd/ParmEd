"""
Tests for the parmed/charmm subpackage
"""
from __future__ import division, print_function

from collections import OrderedDict, defaultdict
import copy
import numpy as np
import os
import parmed as pmd
from parmed.utils.io import genopen
from parmed.utils.six import iteritems, string_types
from parmed.utils.six.moves import StringIO
from parmed.charmm import charmmcrds, parameters, psf
from parmed.charmm._charmmfile import CharmmFile, CharmmStreamFile
from parmed import exceptions, topologyobjects as to, load_file, ParameterSet
from parmed.topologyobjects import BondType, AngleType, DihedralType, DihedralTypeList
import parmed.unit as u
import random
import unittest
from utils import HAS_GROMACS, FileIOTestCase, get_fn, create_random_structure
import warnings

class TestCharmmBase(FileIOTestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.param22 = parameters.CharmmParameterSet(
            get_fn('top_all22_prot.inp'), get_fn('par_all22_prot.inp')
        )

class TestCharmmCoords(TestCharmmBase):
    """ Test CHARMM coordinate file parsers """

    def test_charmm_crd(self):
        """ Test CHARMM coordinate file parser """
        self.assertTrue(charmmcrds.CharmmCrdFile.id_format(get_fn('1tnm.crd')))
        self._check_crd(charmmcrds.CharmmCrdFile(get_fn('1tnm.crd')))
        # Make sure format ID is good
        # Skipped whitespace
        fn = self.get_fn('test.crd', written=True)
        with open(fn, 'w') as f, open(get_fn('1tnm.crd'), 'r') as f2:
            f.write('\n\n\n\n')
            f.write(f2.read())
        self.assertTrue(charmmcrds.CharmmCrdFile.id_format(fn))
        self._check_crd(charmmcrds.CharmmCrdFile(fn))
        with open(fn, 'w') as f, open(get_fn('1tnm.crd'), 'r') as f2:
            f.write(f2.readline())
            f.write(f2.readline())
            f.write(f2.readline())
            f.write('\n')
            f.write(f2.read())
        self.assertTrue(charmmcrds.CharmmCrdFile.id_format(fn))
        self._check_crd(charmmcrds.CharmmCrdFile(fn))
        # Make sure incomplete files are properly error-detected
        with open(fn, 'w') as f, open(get_fn('1tnm.crd'), 'r') as f2:
            for i in range(100):
                f.write(f2.readline())
        self.assertRaises(exceptions.CharmmError, lambda:
                charmmcrds.CharmmCrdFile(fn))

    def _check_crd(self, crd):
        self.assertEqual(crd.natom, 1414)
        self.assertEqual(max(crd.resno), 91)
        self.assertAlmostEqual(crd.coords.sum(), -218.19346999999757)
        self.assertEqual(crd.coords.shape, (1, crd.natom, 3))
        self.assertEqual(len(crd.atomno), crd.natom)
        self.assertEqual(len(crd.resno), crd.natom)
        self.assertEqual(len(crd.resid), crd.natom)
        self.assertEqual(len(crd.resname), crd.natom)
        self.assertEqual(len(crd.weighting), crd.natom)

    def test_write_crd(self):
        """ Test CHARMM coordinate writing capabilities """
        struct = load_file(get_fn('4lzt.pdb'))
        charmmcrds.CharmmCrdFile.write(struct, self.get_fn('test.crd', written=True))
        crd = charmmcrds.CharmmCrdFile(self.get_fn('test.crd', written=True))
        np.testing.assert_allclose(struct.coordinates,
                                   crd.coordinates.reshape((len(struct.atoms), 3)))
        fd = StringIO()
        charmmcrds.CharmmCrdFile.write(struct, fd)
        fd.seek(0)
        with open(self.get_fn('test2.crd', written=True), 'w') as f:
            f.write(fd.read())
        crd = charmmcrds.CharmmCrdFile(self.get_fn('test2.crd', written=True))
        np.testing.assert_allclose(struct.coordinates,
                                   crd.coordinates.reshape((len(struct.atoms), 3)))

    def test_charmm_rst(self):
        """ Test CHARMM restart file parser """
        crd = charmmcrds.CharmmRstFile(get_fn('sample-charmm.rst'))
        self.assertEqual(crd.natom, 256)
        self.assertEqual(crd.nstep, 100)
        self.assertTrue(hasattr(crd, 'header'))
        self.assertAlmostEqual(crd.coords.sum(), 0.3114525961458884)
        self.assertAlmostEqual(crd.coordinates.sum(), 0.3114525961458884)
        self.assertAlmostEqual(crd.coordsold.sum(), 5439.333671681806)
        self.assertAlmostEqual(crd.coordinatesold.sum(), 5439.333671681806)
        self.assertAlmostEqual(crd.vels.sum(), 42.364377359350534)
        self.assertAlmostEqual(crd.velocities.sum(), 42.364377359350534)
        self.assertEqual(crd.coords.shape, (1, crd.natom, 3))
        self.assertEqual(crd.coordsold.shape, (1, crd.natom, 3))
        self.assertEqual(crd.velocities.shape, (1, crd.natom, 3))
        self.assertTrue(u.is_quantity(crd.positions))
        for xyz, pos in zip(crd.coordinates[0], crd.positions):
            np.testing.assert_equal(xyz, pos.value_in_unit(u.angstroms))
        for xyz, pos in zip(crd.coordinatesold[0], crd.positionsold):
            np.testing.assert_equal(xyz, pos.value_in_unit(u.angstroms))
        # Check variables whose meaning I don't understand
        self.assertEqual(crd.jhstrt, 754200)
        self.assertEqual(crd.npriv, 754200)
        self.assertEqual(crd.nsavc, 100)
        self.assertEqual(crd.enrgstat, [])
        self.assertEqual(crd.nsavv, 10)
        self.assertIs(crd.box, None)
        # Check proper handling of truncated files
        fn = self.get_fn('test.rst', written=True)
        with open(fn, 'w') as f, open(get_fn('sample-charmm.rst'), 'r') as f2:
            for i in range(200):
                f.write(f2.readline())
        self.assertRaises(exceptions.CharmmError, lambda:
                charmmcrds.CharmmRstFile(fn))
        with open(fn, 'w') as f:
            f.write('\n\n\n\n\n\n')
        self.assertRaises(exceptions.CharmmError, lambda:
                charmmcrds.CharmmRstFile(fn))

class TestCharmmPsf(TestCharmmBase):
    """ Test CHARMM PSF file capabilities """

    def test_private_internals(self):
        """ Test private internal functions for CHARMM psf file """
        # _catchindexerror
        func = psf._catchindexerror(lambda: [1, 2, 3][10])
        # _ZeroDict
        self.assertRaises(exceptions.CharmmError, func)
        d1 = psf._ZeroDict()
        d2 = psf._ZeroDict()
        d1['NGRP NST2'] = ([1, 1], [1, 2, 3])
        d1['NUMLP NUMLPH'] = ([3, 3], [1, 2, 3])
        d1['a'] = 0
        d1['b'] = 1
        d1['c'] = 2
        self.assertEqual(d1['NGRP'], ([1, 1], [1, 2, 3]))
        self.assertEqual(d2['NGRP'], ([0, 0], []))
        self.assertEqual(d1['NUMLP'], ([3, 3], [1, 2, 3]))
        self.assertEqual(d2['NUMLP'], ([0, 0], []))
        self.assertEqual(d1['a'], 0)
        self.assertEqual(d1['b'], 1)
        self.assertEqual(d1['c'], 2)
        self.assertEqual(d2['a'], (0, []))
        self.assertEqual(d2['b'], (0, []))
        self.assertEqual(d2['c'], (0, []))
        # CharmmPsfFile._convert staticmethod
        self.assertRaises(exceptions.CharmmError, lambda:
                psf.CharmmPsfFile._convert('bad', int, 'not an integer')
        )
        try:
            psf.CharmmPsfFile._convert('bad', int, 'not an integer')
        except exceptions.CharmmError as e:
            self.assertIn('not an integer', str(e))
        else:
            self.assertTrue(False)

    def test_charmm_psf(self):
        """ Test CHARMM PSF file parsing """
        cpsf = psf.CharmmPsfFile(get_fn('ala_ala_ala.psf'))
        self.assertEqual(len(cpsf.atoms), 33)
        for i, atom in enumerate(cpsf.atoms):
            self.assertEqual(atom.idx, i)
            self.assertEqual(atom.residue.name, 'ALA')
            self.assertTrue(atom in atom.residue) # tests __contains__
        # Check the bond, angle, and torsion partners of the first N atom
        a = cpsf.atoms[0]
        for atom in a.bond_partners:
            self.assertTrue(atom.name in ['HT3', 'HT2', 'CA', 'HT1'])
            self.assertEqual(atom.residue.idx, 0)
            self.assertTrue(atom.type in [2, 22])
        for atom in a.angle_partners:
            self.assertTrue(atom.name in ['HA', 'CB', 'C'])
            self.assertTrue(atom.type in [6, 24, 20])
            self.assertEqual(atom.residue.idx, 0)
        for atom in a.dihedral_partners:
            self.assertTrue(atom.name in ['HB1', 'HB2', 'HB3', 'O', 'N'])
            if atom.name == 'N':
                self.assertEqual(atom.residue.idx, 1)
            else:
                self.assertEqual(atom.residue.idx, 0)
            self.assertTrue(atom.type in [3, 70, 54])
        # Check some atom properties
        self.assertRaises(exceptions.ParameterError, lambda: str(a.atom_type))
        self.assertTrue(all([isinstance(b, to.Bond) for b in a.bonds]))
        self.assertTrue(all([isinstance(an, to.Angle) for an in a.angles]))
        self.assertTrue(all([isinstance(d, to.Dihedral) for d in a.dihedrals]))
        self.assertEqual(len(a.angle_partners), 3)
        self.assertEqual(len(a.angles), 9)
        self.assertEqual(len(a.bond_partners), 4)
        self.assertEqual(len(a.bonds), 4)
        self.assertEqual(len(a.cmaps), 0)
        self.assertEqual(len(a.dihedral_partners), 5)
        self.assertEqual(len(a.dihedrals), 14)
        self.assertEqual(len(a.impropers), 0)
        self.assertEqual(len(a.name), 1)
        self.assertEqual(len(a.props), 3)
        self.assertEqual(len(a.residue), 12)
        self.assertEqual(len(a.residue.segid), 3)
        self.assertEqual(len(a.urey_bradleys), 0)
        # Check attributes of the psf file
        self.assertEqual(len(cpsf.acceptors), 4)
        self.assertEqual(len(cpsf.angles), 57)
        self.assertEqual(len(cpsf.atoms), 33)
        self.assertEqual(len(cpsf.bonds), 32)
        self.assertEqual(len(cpsf.cmaps), 1)
        self.assertEqual(len(cpsf.dihedrals), 74)
        self.assertEqual(len(cpsf.donors), 5)
        self.assertEqual(len(cpsf.flags), 2)
        self.assertEqual(len(cpsf.groups), 9)
        self.assertEqual(len(cpsf.impropers), 5)
        self.assertEqual(len(cpsf.residues), 3)
        self.assertEqual(len(cpsf.title), 2)
        # Check the __contains__ methods of valence terms (make sure the correct
        # number of atoms are in each valence term)
        atoms = cpsf.atoms
        bonds = cpsf.bonds
        for bond in cpsf.bonds:
            self.assertEqual(sum([int(a in bond) for a in atoms]), 2)
        # Other valence terms can also contain bonds
        for i, angle in enumerate(cpsf.angles):
            self.assertEqual(sum([int(a in angle) for a in atoms]), 3)
            self.assertEqual(sum([int(b in angle) for b in bonds]), 2)
        for dih in cpsf.dihedrals:
            self.assertEqual(sum([int(a in dih) for a in atoms]), 4)
            self.assertEqual(sum([int(b in dih) for b in bonds]), 3)
        for imp in cpsf.impropers:
            self.assertEqual(sum([int(a in imp) for a in atoms]), 4)
            self.assertEqual(sum([int(b in imp) for b in bonds]), 3)
        for cmap in cpsf.cmaps:
            self.assertEqual(sum([int(a in cmap) for a in atoms]), 5)
            self.assertEqual(sum([int(b in cmap) for b in bonds]), 4)
        # Test CHARMM groups
        g = to.Group(cpsf.groups[0].atom, cpsf.groups[0].type, cpsf.groups[0].move)
        self.assertEqual(g, cpsf.groups[0])
        g.type = 0
        self.assertNotEqual(g, cpsf.groups[0])
        # Check that copying preserves segid attributes
        psf2 = copy.copy(cpsf)
        for r1, r2 in zip(cpsf.residues, psf2.residues):
            self.assertEqual(r1.chain, r2.chain)
            self.assertEqual(r1.segid, r2.segid)
            self.assertEqual(r1.number, r2.number)
            self.assertEqual(r1.idx, r2.idx)
        # Check that slicing preserves segid attributes as well
        firstres = cpsf[0,:]
        self.assertEqual(cpsf.residues[0].segid, firstres.residues[0].segid)
        for res in (firstres + firstres).residues:
            self.assertEqual(res.segid, firstres.residues[0].segid)
        for res in (firstres * 3).residues:
            self.assertEqual(res.segid, firstres.residues[0].segid)

    def test_xplor_psf(self):
        """ Test Xplor-format CHARMM PSF file parsing """
        # Atom types are strings, not integers like in charmm
        cpsf = psf.CharmmPsfFile(get_fn('ala_ala_ala.psf.xplor'))
        self.assertEqual(len(cpsf.atoms), 33)
        for i, atom in enumerate(cpsf.atoms):
            self.assertEqual(atom.idx, i)
            self.assertEqual(atom.residue.name, 'ALA')
            self.assertTrue(atom in atom.residue) # tests __contains__
        # Check the bond, angle, and torsion partners of the first N atom
        a = cpsf.atoms[0]
        for atom in a.bond_partners:
            self.assertTrue(atom.name in ['HT3', 'HT2', 'CA', 'HT1'])
            self.assertEqual(atom.residue.idx, 0)
            self.assertTrue(atom.type in ['HC', 'CT1'])
        for atom in a.angle_partners:
            self.assertTrue(atom.name in ['HA', 'CB', 'C'])
            self.assertTrue(atom.type in ['HB', 'CT3', 'C'])
            self.assertEqual(atom.residue.idx, 0)
        for atom in a.dihedral_partners:
            self.assertTrue(atom.name in ['HB1', 'HB2', 'HB3', 'O', 'N'])
            if atom.name == 'N':
                self.assertEqual(atom.residue.idx, 1)
            else:
                self.assertEqual(atom.residue.idx, 0)
            self.assertTrue(atom.type in ['HA', 'O', 'NH1'])
        # Check some atom properties
        self.assertRaises(exceptions.ParameterError, lambda: int(a.atom_type))
        self.assertTrue(all([isinstance(b, to.Bond) for b in a.bonds]))
        self.assertTrue(all([isinstance(an, to.Angle) for an in a.angles]))
        self.assertTrue(all([isinstance(d, to.Dihedral) for d in a.dihedrals]))
        self.assertEqual(len(a.angle_partners), 3)
        self.assertEqual(len(a.angles), 9)
        self.assertEqual(len(a.bond_partners), 4)
        self.assertEqual(len(a.bonds), 4)
        self.assertEqual(len(a.cmaps), 0)
        self.assertEqual(len(a.dihedral_partners), 5)
        self.assertEqual(len(a.dihedrals), 14)
        self.assertEqual(len(a.impropers), 0)
        self.assertEqual(len(a.name), 1)
        self.assertEqual(len(a.props), 3)
        self.assertEqual(len(a.residue), 12)
        self.assertEqual(len(a.residue.segid), 3)
        self.assertEqual(len(a.urey_bradleys), 0)
        # Check attributes of the psf file
        self.assertEqual(len(cpsf.acceptors), 4)
        self.assertEqual(len(cpsf.angles), 57)
        self.assertEqual(len(cpsf.atoms), 33)
        self.assertEqual(len(cpsf.bonds), 32)
        self.assertEqual(len(cpsf.cmaps), 1)
        self.assertEqual(len(cpsf.dihedrals), 74)
        self.assertEqual(len(cpsf.donors), 5)
        self.assertEqual(len(cpsf.flags), 2)
        self.assertEqual(len(cpsf.groups), 9)
        self.assertEqual(len(cpsf.impropers), 5)
        self.assertEqual(len(cpsf.residues), 3)
        self.assertEqual(len(cpsf.title), 2)
        # Check the __contains__ methods of valence terms (make sure the correct
        # number of atoms are in each valence term)
        atoms = cpsf.atoms
        bonds = cpsf.bonds
        for bond in cpsf.bonds:
            self.assertEqual(sum([int(a in bond) for a in atoms]), 2)
        # Other valence terms can also contain bonds
        for i, angle in enumerate(cpsf.angles):
            self.assertEqual(sum([int(a in angle) for a in atoms]), 3)
            self.assertEqual(sum([int(b in angle) for b in bonds]), 2)
        for dih in cpsf.dihedrals:
            self.assertEqual(sum([int(a in dih) for a in atoms]), 4)
            self.assertEqual(sum([int(b in dih) for b in bonds]), 3)
        for imp in cpsf.impropers:
            self.assertEqual(sum([int(a in imp) for a in atoms]), 4)
            self.assertEqual(sum([int(b in imp) for b in bonds]), 3)
        for cmap in cpsf.cmaps:
            self.assertEqual(sum([int(a in cmap) for a in atoms]), 5)
            self.assertEqual(sum([int(b in cmap) for b in bonds]), 4)

    def test_charmm_gui_builder(self):
        """ Test parsing of CHARMM PSF from CHARMM-GUI """
        cpsf = psf.CharmmPsfFile(get_fn('parv.psf'))
        self.assertEqual(len(cpsf.acceptors), 0)
        self.assertEqual(len(cpsf.angles), 3004)
        self.assertEqual(len(cpsf.atoms), 1659)
        self.assertEqual(len(cpsf.bonds), 1671)
        self.assertEqual(len(cpsf.cmaps), 107)
        self.assertEqual(len(cpsf.dihedrals), 4377)
        self.assertEqual(len(cpsf.donors), 0)
        self.assertEqual(len(cpsf.flags), 3)
        self.assertEqual(len(cpsf.groups), 1)
        self.assertEqual(len(cpsf.impropers), 295)
        self.assertEqual(len(cpsf.residues), 109)
        self.assertEqual(len(cpsf.title), 3)

    def test_vmd_psf(self):
        """ Test parsing of CHARMM PSF from VMD """
        cpsf = psf.CharmmPsfFile(get_fn('ala_ala_ala_autopsf.psf'))
        # Atom types are strings, not integers like in charmm
        self.assertEqual(len(cpsf.atoms), 33)
        for i, atom in enumerate(cpsf.atoms):
            self.assertEqual(atom.idx, i)
            self.assertEqual(atom.residue.name, 'ALA')
            self.assertTrue(atom in atom.residue) # tests __contains__
        # Check the bond, angle, and torsion partners of the first N atom
        a = cpsf.atoms[0]
        for atom in a.bond_partners:
            self.assertTrue(atom.name in ['HT3', 'HT2', 'CA', 'HT1'])
            self.assertEqual(atom.residue.idx, 0)
            self.assertTrue(atom.type in ['HC', 'CT1'])
        for atom in a.angle_partners:
            self.assertTrue(atom.name in ['HA', 'CB', 'C'])
            self.assertTrue(atom.type in ['HB', 'CT3', 'C'])
            self.assertEqual(atom.residue.idx, 0)
        for atom in a.dihedral_partners:
            self.assertTrue(atom.name in ['HB1', 'HB2', 'HB3', 'O', 'N'])
            if atom.name == 'N':
                self.assertEqual(atom.residue.idx, 1)
            else:
                self.assertEqual(atom.residue.idx, 0)
            self.assertTrue(atom.type in ['HA', 'O', 'NH1'])
        # Check some atom properties
        self.assertRaises(exceptions.ParameterError, lambda: int(a.atom_type))
        self.assertTrue(all([isinstance(b, to.Bond) for b in a.bonds]))
        self.assertTrue(all([isinstance(an, to.Angle) for an in a.angles]))
        self.assertTrue(all([isinstance(d, to.Dihedral) for d in a.dihedrals]))
        self.assertEqual(len(a.angle_partners), 3)
        self.assertEqual(len(a.angles), 9)
        self.assertEqual(len(a.bond_partners), 4)
        self.assertEqual(len(a.bonds), 4)
        self.assertEqual(len(a.cmaps), 0)
        self.assertEqual(len(a.dihedral_partners), 5)
        self.assertEqual(len(a.dihedrals), 14)
        self.assertEqual(len(a.impropers), 0)
        self.assertEqual(len(a.name), 1)
        self.assertEqual(len(a.props), 1)
        self.assertEqual(len(a.residue), 12)
        self.assertEqual(len(a.residue.segid), 2)
        self.assertEqual(len(a.urey_bradleys), 0)
        # Check attributes of the psf file
        self.assertEqual(len(cpsf.acceptors), 0)
        self.assertEqual(len(cpsf.angles), 57)
        self.assertEqual(len(cpsf.atoms), 33)
        self.assertEqual(len(cpsf.bonds), 32)
        self.assertEqual(len(cpsf.cmaps), 1)
        self.assertEqual(len(cpsf.dihedrals), 74)
        self.assertEqual(len(cpsf.donors), 0)
        self.assertEqual(len(cpsf.flags), 1)
        self.assertEqual(len(cpsf.groups), 1)
        self.assertEqual(len(cpsf.impropers), 5)
        self.assertEqual(len(cpsf.residues), 3)
        self.assertEqual(len(cpsf.title), 6)
        # Check the __contains__ methods of valence terms (make sure the correct
        # number of atoms are in each valence term)
        atoms = cpsf.atoms
        bonds = cpsf.bonds
        for bond in cpsf.bonds:
            self.assertEqual(sum([int(a in bond) for a in atoms]), 2)
        # Other valence terms can also contain bonds
        for i, angle in enumerate(cpsf.angles):
            self.assertEqual(sum([int(a in angle) for a in atoms]), 3)
            self.assertEqual(sum([int(b in angle) for b in bonds]), 2)
        for dih in cpsf.dihedrals:
            self.assertEqual(sum([int(a in dih) for a in atoms]), 4)
            self.assertEqual(sum([int(b in dih) for b in bonds]), 3)
        for imp in cpsf.impropers:
            self.assertEqual(sum([int(a in imp) for a in atoms]), 4)
            self.assertEqual(sum([int(b in imp) for b in bonds]), 3)
        for cmap in cpsf.cmaps:
            self.assertEqual(sum([int(a in cmap) for a in atoms]), 5)
            self.assertEqual(sum([int(b in cmap) for b in bonds]), 4)

    def test_inscode_psf(self):
        """ Test PSF with insertion code as part of residue number """
        cpsf = psf.CharmmPsfFile(get_fn('4TVP-dmj_wat-ion.psf'))
        self.assertEqual(len(cpsf.atoms), 66264)
        self.assertEqual(len(cpsf.residues), 20169)
        self.assertEqual(len(cpsf.bonds), 46634)
        self.assertEqual(len(cpsf.angles), 32739)
        self.assertEqual(len(cpsf.dihedrals), 19104)
        self.assertEqual(len(cpsf.impropers), 1257)
        self.assertEqual(len(cpsf.cmaps), 447)
        self.assertEqual(cpsf.residues[281].insertion_code, 'A')

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_from_structure(self):
        """ Tests the CharmmPsfFile.from_structure constructor """
        top1 = load_file(get_fn('benzene_cyclohexane_10_500.prmtop'))
        psf1 = psf.CharmmPsfFile.from_structure(top1)

        top2 = load_file(os.path.join(get_fn('03.AlaGlu'), 'topol.top'))
        psf2 = psf.CharmmPsfFile.from_structure(top2)

        self.assertEqual(len(psf1.atoms), len(top1.atoms))
        self.assertEqual(len(psf2.atoms), len(top2.atoms))
        self.assertEqual(len(psf1.residues), len(top1.residues))
        self.assertEqual(len(psf2.residues), len(top2.residues))

        self.assertEqual(len(psf1.bonds), len(top1.bonds))
        self.assertEqual(len(psf2.bonds), len(top2.bonds))
        self.assertEqual(len(psf1.angles), len(top1.angles))
        self.assertEqual(len(psf2.angles), len(top2.angles))
        self.assertEqual(len(psf1.urey_bradleys), len(top1.urey_bradleys))
        self.assertEqual(len(psf2.urey_bradleys), len(top2.urey_bradleys))
        self.assertEqual(len(psf1.dihedrals), len(top1.dihedrals))
        self.assertEqual(len(psf2.dihedrals), len(top2.dihedrals))
        self.assertEqual(len(psf1.impropers), len(top1.impropers))
        self.assertEqual(len(psf2.impropers), len(top2.impropers))
        self.assertEqual(len(psf1.cmaps), len(top1.cmaps))
        self.assertEqual(len(psf2.cmaps), len(top2.cmaps))
        self.assertEqual(len(psf1.acceptors), len(top1.acceptors))
        self.assertEqual(len(psf2.acceptors), len(top2.acceptors))
        self.assertEqual(len(psf1.donors), len(top1.donors))
        self.assertEqual(len(psf2.donors), len(top2.donors))
        self.assertEqual(len(psf1.groups), len(top1.groups))
        self.assertEqual(len(psf2.groups), len(top2.groups))

        self.assertEqual(len(psf1.bond_types), len(top1.bond_types))
        self.assertEqual(len(psf2.bond_types), len(top2.bond_types))
        self.assertEqual(len(psf1.angle_types), len(top1.angle_types))
        self.assertEqual(len(psf2.angle_types), len(top2.angle_types))
        self.assertEqual(len(psf1.dihedral_types), len(top1.dihedral_types))
        self.assertEqual(len(psf2.dihedral_types), len(top2.dihedral_types))
        self.assertEqual(len(psf1.urey_bradley_types), len(top1.urey_bradley_types))
        self.assertEqual(len(psf2.urey_bradley_types), len(top2.urey_bradley_types))
        self.assertEqual(len(psf1.improper_types), len(top1.improper_types))
        self.assertEqual(len(psf2.improper_types), len(top2.improper_types))
        self.assertEqual(len(psf1.cmap_types), len(top1.cmap_types))
        self.assertEqual(len(psf2.cmap_types), len(top2.cmap_types))

        for atom in psf1.atoms:
            self.assertEqual(atom.type.upper(), atom.type)

        # Test the copy argument
        psf3 = psf.CharmmPsfFile.from_structure(top2, copy=True)
        self.assertIsNot(psf3.atoms, top2.atoms)
        self.assertIsNot(psf3.residues, top2.residues)

        self.assertIsNot(psf3.bonds, top2.bonds)
        self.assertIsNot(psf3.angles, top2.angles)
        self.assertIsNot(psf3.urey_bradleys, top2.urey_bradleys)
        self.assertIsNot(psf3.dihedrals, top2.dihedrals)
        self.assertIsNot(psf3.impropers, top2.impropers)
        self.assertIsNot(psf3.cmaps, top2.cmaps)
        self.assertIsNot(psf3.acceptors, top2.acceptors)
        self.assertIsNot(psf3.donors, top2.donors)
        self.assertIsNot(psf3.groups, top2.groups)

        self.assertIsNot(psf3.bond_types, top2.bond_types)
        self.assertIsNot(psf3.angle_types, top2.angle_types)
        self.assertIsNot(psf3.dihedral_types, top2.dihedral_types)
        self.assertIsNot(psf3.urey_bradley_types, top2.urey_bradley_types)
        self.assertIsNot(psf3.improper_types, top2.improper_types)
        self.assertIsNot(psf3.cmap_types, top2.cmap_types)

    def test_error_handling(self):
        """ Tests error handling of CharmmPsfFile """
        self.assertRaises(exceptions.CharmmError, lambda:
                psf.CharmmPsfFile(get_fn('trx.prmtop'))
        )
        # Print some atoms out-of-order
        with open(get_fn('ala_ala_ala.psf'), 'r') as f, \
                open(self.get_fn('ala_ala_ala2.psf', written=True), 'w') as f2:
            for i in range(15):
                f2.write(f.readline())
            tmp = f.readline()
            f2.write(f.readline())
            f2.write(tmp)
            for line in f:
                f2.write(line)
        with self.assertRaises(exceptions.CharmmError):
            psf.CharmmPsfFile(self.get_fn('ala_ala_ala2.psf', written=True))
        # CHARMM can't handle all potential energy functions
        struct = create_random_structure(True)
        self.assertRaises(ValueError, lambda:
                psf.CharmmPsfFile.from_structure(struct)
        )

    def test_copy_parameters(self):
        """ Tests copy_parameters option in load_parameters """
        top = psf.CharmmPsfFile(get_fn('ala_ala_ala.psf'))
        top.load_parameters(parmset=self.param22, copy_parameters=False)
        b = self.param22.bond_types[(top.atoms[0].type, top.atoms[1].type)]
        b.k = 200
        a = self.param22.angle_types[(top.atoms[1].type, top.atoms[0].type, top.atoms[2].type)]
        a.k = 20
        d = self.param22.dihedral_types[('X', top.atoms[4].type, top.atoms[6].type, 'X')]
        d[0].phi_k = 0.300
        self.assertEqual(top.bonds[0].type, self.param22.bond_types[(top.atoms[0].type, top.atoms[1].type)])
        self.assertEqual(top.angles[0].type, self.param22.angle_types[(top.atoms[1].type, top.atoms[0].type, top.atoms[2].type)])
        self.assertEqual(top.dihedrals[0].type, self.param22.dihedral_types[('X', top.atoms[4].type, top.atoms[6].type, 'X')])

        self.param22.bond_types[(top.atoms[0].type, top.atoms[1].type)] = BondType(300, 1.040)
        self.param22.angle_types[(top.atoms[1].type, top.atoms[0].type, top.atoms[2].type)] = AngleType(k=40, theteq=109.5)

        dtl = DihedralTypeList()
        self.param22.dihedral_types[('X', top.atoms[4].type, top.atoms[6].type, 'X')] = \
            dtl.append(DihedralType(phi_k=0.200, per=3, phase=0.00, scee=1.00, scnb=1.00))
        self.assertNotEqual(top.bonds[0].type, self.param22.bond_types[(top.atoms[0].type, top.atoms[1].type)])
        self.assertNotEqual(top.angles[0].type, self.param22.angle_types[(top.atoms[1].type, top.atoms[0].type, top.atoms[2].type)])
        self.assertNotEqual(top.dihedrals[0].type, self.param22.dihedral_types[('X', top.atoms[4].type, top.atoms[6].type, 'X')])

    def test_psf_with_no_nnb_section(self):
        """ Tests parsing of a PSF file with a truncated NNB section """
        top = psf.CharmmPsfFile(get_fn('nonnb.psf'))
        self.assertEqual(len(top.atoms), 10740)

class TestCharmmParameters(TestCharmmBase):
    """ Test CHARMM Parameter file parsing """

    def test_private_functions(self):
        """ Tests private helper functions for CharmmParameterSet """
        # EmptyStringIterator
        si = parameters._EmptyStringIterator()
        it = iter(si)
        for i in range(random.randint(100, 1000)):
            self.assertEqual(next(it), '')
        self.assertEqual(si[random.randint(0, 10000)], '')
        # _typeconv
        randint = random.randint(0, 100000)
        self.assertEqual(parameters._typeconv(randint), randint)
        self.assertEqual(parameters._typeconv('NOCHNG'), 'NOCHNG')
        self.assertEqual(parameters._typeconv('NoCh'), 'NOCHLT')
        self.assertEqual(parameters._typeconv('Na+'), 'NAPLTU')
        self.assertEqual(parameters._typeconv('NA+'), 'NAP')

    def test_e14_fac(self):
        """ Test reading CHARMM parameter files with 1-4 EEL scaling """
        params = parameters.CharmmParameterSet(
                get_fn('parm14sb_all.prm'),
        )
        for i, tortype in iteritems(params.dihedral_types):
            for typ in tortype:
                self.assertAlmostEqual(typ.scee, 1.2)
        params = parameters.CharmmParameterSet(
                get_fn('parm14sb_all_2.prm')
        )
        for i, tortype in iteritems(params.dihedral_types):
            for typ in tortype:
                self.assertAlmostEqual(typ.scee, 1.2)
        # Now test that adding to the parameter set with a DIFFERENT 1-4 scaling
        # factor is caught
        self.assertRaises(exceptions.CharmmError, lambda:
                params.read_parameter_file(get_fn('par_all36_prot.prm'))
        )
        self.assertRaises(exceptions.CharmmError, lambda:
                parameters.CharmmParameterSet(get_fn('parm14sb_all.prm'),
                                              get_fn('dummy_charmm.str'))
        )

    def test_geometric(self):
        """ Test reading CHARMM parameter file with geometric comb. rule """
        opls = parameters.CharmmParameterSet(get_fn('top_opls_aa.inp'),
                                             get_fn('par_opls_aa.inp'))
        self.assertEqual(opls.combining_rule, 'geometric')
        # Now test error handling corresponding to illegal mixing of
        # incompatible parameter files.
        non_opls = parameters.CharmmParameterSet(get_fn('par_all36_prot.prm'))
        self.assertEqual(non_opls.combining_rule, 'lorentz')
        non_opls.read_topology_file(get_fn('top_opls_aa.inp'))
        self.assertRaises(exceptions.CharmmError, lambda:
                non_opls.read_parameter_file(get_fn('par_geometric_combining.inp'))
        )
        self.assertRaises(exceptions.CharmmError, lambda:
                non_opls.read_parameter_file(get_fn('par_opls_aa.inp'))
        )
        for _, dt in iteritems(opls.dihedral_types):
            for t in dt: t.scee = t.scnb = 1.0
        self.assertRaises(exceptions.CharmmError, lambda:
                opls.read_parameter_file(get_fn('par_all36_prot.prm'))
        )

    def test_single_parameterset(self):
        """ Test reading a single parameter set """
        # Make sure we error if trying to load parameters before topology
        with self.assertWarns(exceptions.ParameterWarning):
            parameters.CharmmParameterSet(get_fn('par_all22_prot.inp'))
        # Test error handling for loading files with unsupported extensions
        with self.assertRaises(ValueError):
            parameters.CharmmParameterSet(get_fn('trx.prmtop'))
        with self.assertRaises(ValueError):
            parameters.CharmmParameterSet('x.inp')
        self._check_single_paramset(
            parameters.CharmmParameterSet(get_fn('top_all22_prot.inp'), get_fn('par_all22_prot.inp'))
        )
        self._check_single_paramset(
            parameters.CharmmParameterSet.load_set(
                tfile=get_fn('top_all22_prot.inp'),
                pfile=get_fn('par_all22_prot.inp'),
            )
        )

    def _check_single_paramset(self, params):
        for i, tup in enumerate(params.atom_types_tuple):
            name, num = tup
            self.assertTrue(params.atom_types_tuple[tup] is params.atom_types_str[name])
            self.assertTrue(params.atom_types_tuple[tup] is params.atom_types_int[num])
        self.assertEqual(i, 94) # 95 types, but i starts from 0
        self.assertEqual(len(params.angle_types), 685)
        self.assertEqual(len(params.atom_types_int), 95)
        self.assertEqual(len(params.atom_types_str), 95)
        self.assertEqual(len(params.atom_types_tuple), 95)
        self.assertEqual(len(params.bond_types), 266)
        self.assertEqual(len(params.cmap_types), 12)
        self.assertEqual(len(params.dihedral_types), 772)
        self.assertEqual(len(params.improper_types), 43)
        self.assertEqual(len(params.nbfix_types), 0)
        self.assertEqual(len(params.parametersets), 1)
        self.assertEqual(len(params.urey_bradley_types), 685)
        # Look at the parsed residue templates and make sure they line up
        self.assertEqual(len(params.residues), 32)
        self.assertEqual(set(params.residues.keys()),
                set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                     'HSD', 'HSE', 'HSP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                     'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ALAD', 'TIP3',
                     'TP3M', 'SOD', 'MG', 'POT', 'CES', 'CAL', 'CLA', 'ZN2'])
        )
        self.assertEqual(len(params.patches), 22)
        for resname, res in iteritems(params.residues):
            if resname in ('TIP3', 'TP3M', 'SOD', 'MG', 'CLA', 'POT', 'CES',
                    'CAL', 'ZN2', 'ALAD'):
                self.assertIs(res.first_patch, None)
                self.assertIs(res.last_patch, None)
                continue
            self.assertIs(res.last_patch, params.patches['CTER'])
            if resname == 'GLY':
                self.assertIs(res.first_patch, params.patches['GLYP'])
            elif resname == 'PRO':
                self.assertIs(res.first_patch, params.patches['PROP'])
            else:
                self.assertIs(res.first_patch, params.patches['NTER'])
        # Look at the number of unique terms
        def uniques(stuff):
            return len({id(value) for value in stuff.values()})
        self.assertEqual(uniques(params.angle_types), 356)
        self.assertEqual(uniques(params.atom_types_int), 95)
        self.assertEqual(uniques(params.atom_types_str), 95)
        self.assertEqual(uniques(params.atom_types_tuple), 95)
        self.assertEqual(uniques(params.bond_types), 140)
        self.assertEqual(uniques(params.cmap_types), 6)
        self.assertEqual(uniques(params.dihedral_types), 396)
        self.assertEqual(uniques(params.improper_types), 43)
        self.assertEqual(uniques(params.urey_bradley_types), 105)
        obj = params.condense()
        self.assertTrue(obj is params)
        # Check that condensing happened
        self.assertEqual(uniques(params.angle_types), 178)
        self.assertEqual(uniques(params.atom_types_int), 95)
        self.assertEqual(uniques(params.atom_types_str), 95)
        self.assertEqual(uniques(params.atom_types_tuple), 95)
        self.assertEqual(uniques(params.bond_types), 103)
        self.assertEqual(uniques(params.cmap_types), 3)
        self.assertEqual(uniques(params.dihedral_types), 81)
        self.assertEqual(uniques(params.improper_types), 20)
        self.assertEqual(uniques(params.urey_bradley_types), 42)
        # Make sure all cmaps have 8 atom type keys
        for key in params.cmap_types:
            self.assertEqual(len(key), 8)

    def test_param_file_only(self):
        """ Test reading only a parameter file with no RTF (CHARMM36) """
        parameters.CharmmParameterSet(get_fn('par_all36_carb.prm')).condense()
        # Make sure read_parameter_file can accept a list of lines *without*
        # comments
        with CharmmFile(get_fn('par_all36_carb.prm'), 'r') as f:
            params = parameters.CharmmParameterSet()
            params.read_parameter_file(f.readlines())

    def test_collection(self):
        """ Test reading a large number of parameter files """
        p = parameters.CharmmParameterSet(
                    get_fn('top_all36_prot.rtf'),
                    get_fn('top_all36_carb.rtf'),
                    get_fn('par_all36_prot.prm'),
                    get_fn('par_all36_carb.prm'),
                    get_fn('toppar_water_ions.str'),
        ).condense()
        # Look at the number of unique terms
        def uniques(stuff):
            myset = set()
            for key in stuff: myset.add(id(stuff[key]))
            return len(myset)
        self.assertEqual(uniques(p.angle_types), 206)
        self.assertEqual(uniques(p.atom_types_int), 123)
        self.assertEqual(uniques(p.atom_types_str), 123)
        self.assertEqual(uniques(p.atom_types_tuple), 123)
        self.assertEqual(uniques(p.bond_types), 114)
        self.assertEqual(uniques(p.cmap_types), 3)
        self.assertEqual(uniques(p.dihedral_types), 257)
        self.assertEqual(uniques(p.improper_types), 15)
        self.assertEqual(uniques(p.nbfix_types), 6)
        self.assertEqual(uniques(p.urey_bradley_types), 45)
        for key in p.cmap_types:
            self.assertEqual(len(key), 8)

    def test_write_params(self):
        """ Tests writing CHARMM RTF/PAR/STR files from parameter sets """
        params = parameters.CharmmParameterSet(get_fn('top_all22_prot.inp'), get_fn('par_all22_prot.inp'))
        params.write(top=self.get_fn('test.rtf', written=True),
                     par=self.get_fn('test.par', written=True))
        params.write(str=self.get_fn('test.str', written=True))
        # Check bad options
        self.assertRaises(ValueError, lambda: params.write())

        params2 = parameters.CharmmParameterSet(
            self.get_fn('test.rtf', written=True), self.get_fn('test.par', written=True)
        )
        params3 = parameters.CharmmParameterSet(self.get_fn('test.str', written=True))
        params4 = parameters.CharmmParameterSet.load_set(sfiles=self.get_fn('test.str', written=True))
        params5 = parameters.CharmmParameterSet.load_set(sfiles=[self.get_fn('test.str', written=True)])

        # Check that all of the params are equal
        self._compare_paramsets(params, params2, copy=True)
        self._compare_paramsets(params, params3, copy=True)
        self._compare_paramsets(params, params4, copy=True)
        self._compare_paramsets(params, params5, copy=True)

    def test_cgenff(self):
        """ Test parsing stream files generated by CGenFF """
        p = parameters.CharmmParameterSet(get_fn('toppar_spin_label_dummy.str'))
        p = p.condense()
        self.assertEqual(len(p.atom_types_str), 4)
        self.assertEqual(p.atom_types_str['CBD'].epsilon, 0)
        self.assertEqual(p.atom_types_str['CBD'].rmin, 0)
        self.assertAlmostEqual(p.atom_types_str['OND'].epsilon, -0.05)
        self.assertAlmostEqual(p.atom_types_str['OND'].rmin, 2.0)

    def test_penalty(self):
        """ Test parsing penalty scores for CGenFF parameters from comments """
        p = parameters.CharmmParameterSet(get_fn('pyrrol.str'))
        # Check bond types
        self.assertEqual(p.bond_types[('CG251O', 'CG2D2')].penalty, 210)
        self.assertEqual(p.bond_types[('CG251O', 'CG3C52')].penalty, 64)
        self.assertEqual(p.bond_types[('CG251O', 'NG3C51')].penalty, 64)
        self.assertEqual(p.bond_types[('CG321', 'NG3C51')].penalty, 20)
        # Check angle types
        self.assertEqual(p.angle_types[('CG2D2','CG251O','CG3C52')].penalty, 268)
        self.assertEqual(p.angle_types[('CG2D2','CG251O','NG3C51')].penalty, 264.5)
        self.assertEqual(p.angle_types[('CG3C52','CG251O','NG3C51')].penalty, 328)
        self.assertEqual(p.angle_types[('CG251O','CG2D2','HGA5')].penalty, 23)
        self.assertEqual(p.angle_types[('CG321','CG321','NG3C51')].penalty, 3.9)
        self.assertEqual(p.angle_types[('NG3C51','CG321','HGA2')].penalty, 3)
        self.assertEqual(p.angle_types[('CG251O','CG3C52','CG2R51')].penalty, 4)
        self.assertEqual(p.angle_types[('CG251O','CG3C52','HGA2')].penalty, 3.5)
        self.assertEqual(p.angle_types[('CG251O','NG3C51','CG2R51')].penalty, 75)
        self.assertEqual(p.angle_types[('CG251O','NG3C51','CG321')].penalty, 123.5)
        self.assertEqual(p.angle_types[('CG2R51','NG3C51','CG321')].penalty, 122.5)
        # Check dihedral types
        self.assertEqual(p.dihedral_types[('CG3C52','CG251O','CG2D2','HGA5')].penalty, 295)
        self.assertEqual(p.dihedral_types[('NG3C51','CG251O','CG2D2','HGA5')].penalty, 294)
        self.assertEqual(p.dihedral_types[('CG2D2','CG251O','CG3C52','CG2R51')].penalty, 136.5)
        self.assertEqual(p.dihedral_types[('CG2D2','CG251O','CG3C52','HGA2')].penalty, 136.5)
        self.assertEqual(p.dihedral_types[('NG3C51','CG251O','CG3C52','CG2R51')].penalty, 165.5)
        self.assertEqual(p.dihedral_types[('NG3C51','CG251O','CG3C52','HGA2')].penalty, 148)
        self.assertEqual(p.dihedral_types[('CG2D2','CG251O','NG3C51','CG2R51')].penalty, 210.5)
        self.assertEqual(p.dihedral_types[('CG2D2','CG251O','NG3C51','CG321')].penalty, 167.5)
        self.assertEqual(p.dihedral_types[('CG3C52','CG251O','NG3C51','CG2R51')].penalty, 232.5)
        self.assertEqual(p.dihedral_types[('CG3C52','CG251O','NG3C51','CG321')].penalty, 189.5)
        self.assertEqual(p.dihedral_types[('CG2R51','CG2R51','CG3C52','CG251O')].penalty, 4)
        self.assertEqual(p.dihedral_types[('HGR51','CG2R51','CG3C52','CG251O')].penalty, 4)
        self.assertEqual(p.dihedral_types[('CG2R51','CG2R51','NG3C51','CG251O')].penalty, 75)
        self.assertEqual(p.dihedral_types[('CG2R51','CG2R51','NG3C51','CG321')].penalty, 31)
        self.assertEqual(p.dihedral_types[('HGR52','CG2R51','NG3C51','CG251O')].penalty, 75)
        self.assertEqual(p.dihedral_types[('HGR52','CG2R51','NG3C51','CG321')].penalty, 31)
        self.assertEqual(p.dihedral_types[('CG331','CG321','CG321','NG3C51')].penalty, 33)
        self.assertEqual(p.dihedral_types[('CG331','CG321','CG321','NG3C51')].penalty, 33)
        self.assertEqual(p.dihedral_types[('NG3C51','CG321','CG321','HGA2')].penalty, 9)
        self.assertEqual(p.dihedral_types[('CG321','CG321','NG3C51','CG251O')].penalty, 89)
        self.assertEqual(p.dihedral_types[('CG321','CG321','NG3C51','CG251O')].penalty, 89)
        self.assertEqual(p.dihedral_types[('CG321','CG321','NG3C51','CG251O')].penalty, 89)
        self.assertEqual(p.dihedral_types[('CG321','CG321','NG3C51','CG2R51')].penalty, 88)
        self.assertEqual(p.dihedral_types[('CG321','CG321','NG3C51','CG2R51')].penalty, 88)
        self.assertEqual(p.dihedral_types[('CG321','CG321','NG3C51','CG2R51')].penalty, 88)
        self.assertEqual(p.dihedral_types[('HGA2','CG321','NG3C51','CG251O')].penalty, 49.5)
        self.assertEqual(p.dihedral_types[('HGA2','CG321','NG3C51','CG2R51')].penalty, 48.5)

    @unittest.skipUnless(HAS_GROMACS, "Cannot run GROMACS tests without GROMACS")
    def test_charmm_parameter_set_conversion(self):
        """ Tests CharmmParameterSet.from_parameterset and from_structure """
        params1 = ParameterSet.from_structure(
                load_file(get_fn('benzene_cyclohexane_10_500.prmtop'))
        )
        params2 = load_file(os.path.join(get_fn('03.AlaGlu'), 'topol.top')).parameterset

        chparams1 = parameters.CharmmParameterSet.from_parameterset(params1)
        chparams2 = parameters.CharmmParameterSet.from_parameterset(params2, copy=True)
        chparams3 = parameters.CharmmParameterSet.from_structure(
                load_file(get_fn('benzene_cyclohexane_10_500.prmtop'))
        )

        self.assertIsInstance(chparams1, parameters.CharmmParameterSet)
        self.assertIsInstance(chparams2, parameters.CharmmParameterSet)
        self.assertIsInstance(chparams3, parameters.CharmmParameterSet)

        self._compare_paramsets(chparams1, params1, copy=False)
        self._compare_paramsets(chparams2, params2, copy=True)
        self._compare_paramsets(chparams1, chparams3, copy=True)

        self._check_uppercase_types(chparams1)
        self._check_uppercase_types(chparams2)
        self._check_uppercase_types(chparams3)

        # GAFF atom types, as in the first parameter set, are all lower-case.
        # Check that name decoration is the established pattern
        for name in chparams1.atom_types:
            self.assertTrue(name.endswith('LTU'))
        for name in chparams3.atom_types:
            self.assertTrue(name.endswith('LTU'))

        # Load a parameter set with NoUreyBradley to make sure it's retained as
        # a singleton. Also build a list of atom type tuples
        for i, (typstr, typ) in enumerate(iteritems(params1.atom_types)):
            params1.atom_types_tuple[(typstr, i+1)] = typ
            params1.atom_types_int[i+1] = typ
        for key in params1.angle_types:
            params1.urey_bradley_types[key] = to.NoUreyBradley
        chparams = parameters.CharmmParameterSet.from_parameterset(params1)
        for _, item in iteritems(chparams.urey_bradley_types):
            self.assertIs(item, to.NoUreyBradley)
        for (typstr1, typ1), (typstr2, typ2) in zip(iteritems(params1.atom_types),
                                                    iteritems(params1.atom_types)):
            self.assertEqual(typstr1, typstr2)
            self.assertEqual(typ1, typ2)
        for (typstr1, typ1), (typstr2, typ2) in zip(iteritems(params1.atom_types_int),
                                                    iteritems(params1.atom_types_int)):
            self.assertEqual(typstr1, typstr2)
            self.assertEqual(typ1, typ2)
        for (typstr1, typ1), (typstr2, typ2) in zip(iteritems(params1.atom_types_tuple),
                                                    iteritems(params1.atom_types_tuple)):
            self.assertEqual(typstr1, typstr2)
            self.assertEqual(typ1, typ2)

        # Convert from the GROMACS topology file parameter set to
        # CharmmParameterSet when loaded from the CHARMM force field
        gmx = pmd.gromacs.GromacsTopologyFile(
                os.path.join(pmd.gromacs.GROMACS_TOPDIR,
                             'charmm27.ff', 'forcefield.itp')
        )
        from_gmx = parameters.CharmmParameterSet.from_parameterset(gmx.parameterset)
        gmx = pmd.gromacs.GromacsTopologyFile(
                os.path.join(pmd.gromacs.GROMACS_TOPDIR,
                             'charmm27.ff', 'forcefield.itp')
        )
        gmx.parameterset.nbfix_types[('X', 'Y')] = (2.0, 3.0)
        from_gmx2 = parameters.CharmmParameterSet.from_parameterset(gmx.parameterset)
        for (key1, typ1), (key2, typ2) in zip(iteritems(from_gmx.cmap_types),
                                              iteritems(from_gmx2.cmap_types)):
            self.assertEqual(key1, key2)
            self.assertEqual(typ1, typ2)
        self.assertEqual(len(from_gmx2.nbfix_types), 1)
        self.assertEqual(from_gmx2.nbfix_types[('X', 'Y')], (2.0, 3.0))

    def test_parameters_from_structure(self):
        """ Test creation of CharmmParameterSet from a Structure """
        top = psf.CharmmPsfFile(get_fn('ala_ala_ala.psf'))
        top.load_parameters(self.param22)
        params = parameters.CharmmParameterSet.from_structure(top)
        self.assertGreater(len(params.urey_bradley_types), 0)
        for key in params.urey_bradley_types:
            self.assertEqual(len(key), 3)

    def test_warning(self):
        """ Tests warning when overwriting parameters"""
        with self.assertWarns(exceptions.ParameterWarning):
            parameters.CharmmParameterSet(
                get_fn('toppar_all36_prot_aldehydes.str'),
                get_fn('toppar_all36_na_modifications.str'),
            )


    def _check_uppercase_types(self, params):
        for aname, atom_type in iteritems(params.atom_types):
            self.assertEqual(aname, aname.upper())
            self.assertEqual(atom_type.name, atom_type.name.upper())
        for key in params.bond_types:
            for k in key:
                self.assertEqual(k.upper(), k)
        for key in params.angle_types:
            for k in key:
                self.assertEqual(k.upper(), k)
        for key in params.dihedral_types:
            for k in key:
                self.assertEqual(k.upper(), k)
        for key in params.cmap_types:
            for k in key:
                self.assertEqual(k.upper(), k)

    def _compare_paramsets(self, set1, set2, copy):
        def get_typeset(set1, set2):
            ids1 = set()
            ids2 = set()
            for _, item in iteritems(set1):
                ids1.add(id(item))
            for _, item in iteritems(set2):
                ids2.add(id(item))
            return ids1, ids2
        def typenames(key):
            if isinstance(key, string_types):
                return parameters._typeconv(key)
            return tuple(typenames(k) for k in key)
        # Bonds
        b1, b2 = get_typeset(set1.bond_types, set2.bond_types)
        self.assertEqual(len(b1), len(b2))
        if copy:
            self.assertFalse(b1 & b2)
        else:
            self.assertEqual(b1, b2)
        for key, item2 in iteritems(set2.bond_types):
            self.assertEqual(set1.bond_types[typenames(key)], item2)
        # Angles
        a1, a2 = get_typeset(set1.angle_types, set2.angle_types)
        self.assertEqual(len(a1), len(a2))
        if copy:
            self.assertFalse(a1 & a2)
        else:
            self.assertEqual(a1, a2)
        for key, item2 in iteritems(set2.angle_types):
            self.assertEqual(set1.angle_types[typenames(key)], item2)
        # Dihedrals
        d1, d2 = get_typeset(set1.dihedral_types, set2.dihedral_types)
        self.assertEqual(len(d1), len(d2))
        if copy:
            self.assertFalse(d1 & d2)
        else:
            self.assertEqual(d1, d2)
        for key, item2 in iteritems(set2.dihedral_types):
            self.assertEqual(set1.dihedral_types[typenames(key)], item2)
        # Impropers
        d1, d2 = get_typeset(set1.improper_types, set2.improper_types)
        self.assertEqual(len(d1), len(d2))
        if copy:
            self.assertFalse(d1 & d2)
        else:
            self.assertEqual(d1, d2)
        for key, item2 in iteritems(set2.improper_types):
            self.assertEqual(set1.improper_types[typenames(key)], item2)
        # Periodic impropers
        d1, d2 = get_typeset(set1.improper_periodic_types, set2.improper_periodic_types)
        self.assertEqual(len(d1), len(d2))
        if copy:
            self.assertFalse(d1 & d2)
        else:
            self.assertEqual(d1, d2)
        for key, item2 in iteritems(set2.improper_periodic_types):
            self.assertEqual(set1.improper_periodic_types[typenames(key)], item2)
        # CMAPs
        d1, d2 = get_typeset(set1.cmap_types, set2.cmap_types)
        self.assertEqual(len(d1), len(d2))
        if copy:
            self.assertFalse(d1 & d2)
        else:
            self.assertEqual(d1, d2)
        for key, item2 in iteritems(set2.cmap_types):
            self.assertEqual(len(key), 8)
            self.assertEqual(set1.cmap_types[typenames(key)], item2)
        # Atom types
        a1, a2 = get_typeset(set1.atom_types, set2.atom_types)
        ndups = 0
        recorded_types = set()
        for typ in set2.atom_types:
            if parameters._typeconv(typ) in recorded_types:
                ndups += 1
            recorded_types.add(parameters._typeconv(typ))
        self.assertEqual(len(a1), len(a2)-ndups)
        if copy:
            self.assertFalse(a1 & a2)
        else:
            self.assertEqual(a1, a2)

    def test_charmm36_rtf(self):
        """Test parsing of CHARMM36 RTF files."""
        # Make sure there are no failures loading CHARMM36 RTF files.
        param36 = parameters.CharmmParameterSet(get_fn('top_all36_prot.rtf'),
                                                get_fn('top_all36_carb.rtf'),
                                                get_fn('top_all36_cgenff.rtf'))

class TestFileWriting(TestCharmmBase):
    """ Tests the various file writing capabilities """

    def test_charmm_file(self):
        """ Test the CharmmFile API and error handling """
        self.assertRaises(ValueError, lambda:
                CharmmFile(get_fn('trx.prmtop'), 'x')
        )
        self.assertRaises(IOError, lambda:
                CharmmFile(get_fn('file_does_not_exist'), 'r')
        )
        with CharmmFile(self.get_fn('newfile.chm', written=True), 'w') as f:
            f.write('abc123\ndef456\nghi789!comment...\n')
        with CharmmFile(self.get_fn('newfile.chm', written=True), 'r') as f:
            self.assertEqual(f.read(), 'abc123\ndef456\nghi789\n')
        with CharmmFile(get_fn('trx.prmtop')) as f1, open(get_fn('trx.prmtop')) as f2:
            firstline = f1.readline()
            self.assertEqual(firstline, f2.readline())
            self.assertEqual(f1.tell(), f2.tell())
            f1.seek(0)
            f2.seek(0)
            firstline2 = f2.readline()
            self.assertEqual(f1.readline(), firstline2)
            self.assertEqual(firstline, firstline2)
            f1.rewind()
            self.assertEqual(f1.readline(), firstline)
        # Now make sure that every way of opening/reading a file in CharmmFile
        # gets rid of ! comments
        with open(self.get_fn('test.chm', written=True), 'w') as f:
            f.write('First line ! first comment\n'
                    'Second line ! second comment\n'
                    'Third line ! third comment\n'
                    'Fourth line ! fourth comment\n')
        with CharmmFile(self.get_fn('test.chm', written=True), 'r') as f:
            lines = []
            comments = []
            line = f.readline()
            while line:
                lines.append(line)
                comments.append(f.comment)
                line = f.readline()
            self.assertEqual(lines, ['First line \n', 'Second line \n',
                                     'Third line \n', 'Fourth line \n'])
            self.assertEqual(comments, ['! first comment', '! second comment',
                                        '! third comment', '! fourth comment'])

    def test_charmm_stream_file(self):
        """ Test the CharmmStreamFile API """
        stream = CharmmStreamFile(get_fn('toppar_spin_label_dummy.str'))
        lines = genopen(get_fn('toppar_spin_label_dummy.str'), 'r').readlines()
        for l1, l2, c in zip(stream, lines, stream.comments):
            if '!' in l2:
                self.assertEqual(l1, l2[:l2.index('!')] + '\n')
                self.assertEqual(c, l2[l2.index('!'):].rstrip())
            else:
                self.assertEqual(l1, l2)
                self.assertEqual(c, '')
        stream.rewind()
        if '!' in lines[0]:
            self.assertEqual(next(iter(stream)), lines[0][:lines[0].index('!')] + '\n')
        else:
            self.assertEqual(next(iter(stream)), lines[0])

    def test_write_simple_psf(self):
        """ Test writing simple PSF files """
        cpsf = psf.CharmmPsfFile(get_fn('ala_ala_ala.psf'))
        cpsf.flags = [f for f in cpsf.flags if f != 'EXT'] # NO EXT!
        fn = self.get_fn('test.psf', written=True)
        cpsf.write_psf(fn)
        cpsf2 = psf.CharmmPsfFile(fn)

    def test_eliminate_duplicate_dihedrals(self):
        """ Test that duplicate torsions are eliminated in PSF writes """
        def count_torsions(parm):
            torsions = defaultdict(int)
            for d in parm.dihedrals:
                if d.improper: continue # Skip impropers
                if d.atom1 > d.atom4:
                    torsions[(d.atom4.idx, d.atom3.idx, d.atom2.idx, d.atom1.idx)] += 1
                else:
                    torsions[(d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)] += 1
            return torsions
        fn = self.get_fn('test.psf', written=True)
        parm = load_file(get_fn('trx.prmtop'))
        ptorsions = count_torsions(parm)
        parm.write_psf(fn)
        cpsf = psf.CharmmPsfFile(fn)
        ctorsions = count_torsions(cpsf)
        self.assertGreater(max(ptorsions.values()), 1)
        self.assertEqual(set(ctorsions.keys()), set(ptorsions.keys()))
        self.assertEqual(max(ctorsions.values()), 1)

    def test_write_charmm(self):
        """ Test writing CHARMM-style PSF files """
        # Test writing CHARMM-style PSFs
        cpsf = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf'))
        cpsf.write_psf(self.get_fn('dhfr_cmap_pbc.psf', written=True))
        cpsf2 = psf.CharmmPsfFile(self.get_fn('dhfr_cmap_pbc.psf', written=True))
        for attr in dir(cpsf):
            if attr.startswith('_'): continue
            # Skip descriptors
            if attr in ('topology', 'positions', 'box_vectors',
                        'velocities', 'name', 'view'):
                continue
            if callable(getattr(cpsf, attr)): continue
            if hasattr(getattr(cpsf, attr), '__len__'):
                self.assertEqual(len(getattr(cpsf, attr)),
                                 len(getattr(cpsf2, attr)))
            else:
                self.assertEqual(getattr(cpsf, attr), getattr(cpsf2, attr))
        f = open(self.get_fn('dhfr_cmap_pbc.psf', written=True), 'r')
        try:
            has_key = False
            for line in f:
                if '!MOLNT' in line:
                    has_key = True
                    break
        finally:
            f.close()
        self.assertTrue(has_key)

    def test_write_vmd(self):
        """ Test writing VMD-style PSF files """
        # Test writing VMD-style PSFs
        cpsf = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf'))
        cpsf.write_psf(self.get_fn('dhfr_cmap_pbc.psf', written=True), vmd=True)
        cpsf2 = psf.CharmmPsfFile(self.get_fn('dhfr_cmap_pbc.psf', written=True))
        for attr in dir(cpsf):
            if attr.startswith('_'): continue
            if attr in ('topology', 'positions', 'box_vectors',
                        'velocities', 'name', 'view'):
                continue
            if callable(getattr(cpsf, attr)): continue
            if hasattr(getattr(cpsf, attr), '__len__'):
                self.assertEqual(len(getattr(cpsf, attr)),
                                 len(getattr(cpsf2, attr)))
            else:
                self.assertEqual(getattr(cpsf, attr), getattr(cpsf2, attr))
        f = open(self.get_fn('dhfr_cmap_pbc.psf', written=True), 'r')
        try:
            has_key = False
            for line in f:
                if '!MOLNT' in line:
                    has_key = True
                    break
        finally:
            f.close()
        self.assertFalse(has_key)

    def test_write_xplor(self):
        """ Test that XPLOR-style CHARMM PSF files have XPLOR flag (#715) """
        parm = pmd.load_file(get_fn('trx.prmtop'))
        fn = self.get_fn('test.psf', written=True)
        parm.save(fn, overwrite=True)
        cpsf = pmd.load_file(fn)
        self.assertIn('XPLOR', cpsf.flags)
