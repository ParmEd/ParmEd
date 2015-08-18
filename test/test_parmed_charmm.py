"""
Tests for the parmed/charmm subpackage
"""
from __future__ import division, print_function

import numpy as np
from parmed.utils.six import iteritems, string_types
from parmed.utils.six.moves import StringIO
from parmed.charmm import charmmcrds, parameters, psf
from parmed import exceptions, topologyobjects as to, load_file, ParameterSet
import os
import unittest
import utils

get_fn = utils.get_fn

class TestCharmmCoords(utils.FileIOTestCase):
    """ Test CHARMM coordinate file parsers """
    
    def testCharmmCrd(self):
        """ Test CHARMM coordinate file parser """
        crd = charmmcrds.CharmmCrdFile(get_fn('1tnm.crd'))
        self.assertEqual(crd.natom, 1414)
        self.assertEqual(max(crd.resno), 91)
        self.assertAlmostEqual(crd.coords.sum(), -218.19346999999757)
        self.assertEqual(crd.coords.shape, (1, crd.natom, 3))
        self.assertEqual(len(crd.atomno), crd.natom)
        self.assertEqual(len(crd.resno), crd.natom)
        self.assertEqual(len(crd.resid), crd.natom)
        self.assertEqual(len(crd.resname), crd.natom)
        self.assertEqual(len(crd.weighting), crd.natom)

    def testWriteCrd(self):
        """ Test CHARMM coordinate writing capabilities """
        struct = load_file(get_fn('4lzt.pdb'))
        charmmcrds.CharmmCrdFile.write(struct, get_fn('test.crd', written=True))
        crd = charmmcrds.CharmmCrdFile(get_fn('test.crd', written=True))
        np.testing.assert_allclose(struct.coordinates,
                                   crd.coordinates.reshape((len(struct.atoms), 3)))

    def testCharmmRst(self):
        """ Test CHARMM restart file parser """
        crd = charmmcrds.CharmmRstFile(get_fn('sample-charmm.rst'))
        self.assertEqual(crd.natom, 256)
        self.assertEqual(crd.nstep, 100)
        self.assertTrue(hasattr(crd, 'header'))
        self.assertAlmostEqual(crd.coords.sum(), 0.3114525961458884)
        self.assertAlmostEqual(crd.coordsold.sum(), 5439.333671681806)
        self.assertAlmostEqual(crd.vels.sum(), 42.364377359350534)
        self.assertEqual(crd.coords.shape, (1, crd.natom, 3))
        self.assertEqual(crd.coordsold.shape, (1, crd.natom, 3))
        self.assertEqual(crd.velocities.shape, (1, crd.natom, 3))
        # Check variables whose meaning I don't understand
        self.assertEqual(crd.jhstrt, 754200)
        self.assertEqual(crd.npriv, 754200)
        self.assertEqual(crd.nsavc, 100)
        self.assertEqual(crd.enrgstat, [])
        self.assertEqual(crd.nsavv, 10)

class TestCharmmPsf(unittest.TestCase):
    """ Test CHARMM PSF file capabilities """
    
    def testCharmmPsf(self):
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
        self.assertEqual(len(a.segid), 3)
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

    def testXplorPsf(self):
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
        self.assertEqual(len(a.segid), 3)
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

    def testCharmmGuiBuilder(self):
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
    
    def testVmdPsf(self):
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
        self.assertEqual(len(a.segid), 2)
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

    def testInscodePSF(self):
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

    def testFromStructure(self):
        """ Tests the CharmmPsfFile.from_structure constructor """
        top1 = load_file(get_fn('benzene_cyclohexane_10_500.prmtop'))
        psf1 = psf.CharmmPsfFile.from_structure(top1)

        top2 = load_file(os.path.join(get_fn('03.AlaGlu'), 'topol.top'))
        psf2 = psf.CharmmPsfFile.from_structure(top2)

        self.assertEqual(len(psf1.atoms), len(top1.atoms))
        self.assertEqual(len(psf2.atoms), len(top2.atoms))

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

class TestCharmmParameters(utils.FileIOTestCase):
    """ Test CHARMM Parameter file parsing """
    
    def testSingleParameterset(self):
        """ Test reading a single parameter set """
        self.assertRaises(RuntimeError, lambda: parameters.CharmmParameterSet(
                                                get_fn('par_all22_prot.inp')))
        params = parameters.CharmmParameterSet(
                                get_fn('top_all22_prot.inp'),
                                get_fn('par_all22_prot.inp'),
        )
        for i, tup in enumerate(params.atom_types_tuple):
            name, num = tup
            self.assertTrue(params.atom_types_tuple[tup] is
                            params.atom_types_str[name])
            self.assertTrue(params.atom_types_tuple[tup] is
                            params.atom_types_int[num])
        self.assertEqual(i, 94) # 95 types, but i starts from 0
        self.assertEqual(len(params.angle_types), 685)
        self.assertEqual(len(params.atom_types_int), 95)
        self.assertEqual(len(params.atom_types_str), 95)
        self.assertEqual(len(params.atom_types_tuple), 95)
        self.assertEqual(len(params.bond_types), 266)
        self.assertEqual(len(params.cmap_types), 12)
        self.assertEqual(len(params.dihedral_types), 772)
        self.assertEqual(len(params.improper_types), 33)
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
            myset = set()
            for key in stuff: myset.add(id(stuff[key]))
            return len(myset)
        self.assertEqual(uniques(params.angle_types), 356)
        self.assertEqual(uniques(params.atom_types_int), 95)
        self.assertEqual(uniques(params.atom_types_str), 95)
        self.assertEqual(uniques(params.atom_types_tuple), 95)
        self.assertEqual(uniques(params.bond_types), 140)
        self.assertEqual(uniques(params.cmap_types), 6)
        self.assertEqual(uniques(params.dihedral_types), 396)
        self.assertEqual(uniques(params.improper_types), 33)
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

    def testParamFileOnly(self):
        """ Test reading only a parameter file with no RTF (CHARMM36) """
        parameters.CharmmParameterSet(get_fn('par_all36_carb.prm')).condense()

    def testCollection(self):
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

    def testWriteParams(self):
        """ Tests writing CHARMM RTF/PAR/STR files from parameter sets """
        params = parameters.CharmmParameterSet(
                                get_fn('top_all22_prot.inp'),
                                get_fn('par_all22_prot.inp'),
        )
        params.write(top=get_fn('test.rtf', written=True),
                     par=get_fn('test.par', written=True))
        params.write(str=get_fn('test.str', written=True))

        params2 = parameters.CharmmParameterSet(
                                get_fn('test.rtf', written=True),
                                get_fn('test.par', written=True)
        )
        params3 = parameters.CharmmParameterSet(get_fn('test.str', written=True))

        # Check that all of the params are equal
        self._compare_paramsets(params, params2, copy=True)
        self._compare_paramsets(params, params3, copy=True)

    def testCGenFF(self):
        """ Test parsing stream files generated by CGenFF """
        p = parameters.CharmmParameterSet(get_fn('toppar_spin_label_dummy.str'))
        p = p.condense()
        self.assertEqual(len(p.atom_types_str), 4)
        self.assertEqual(p.atom_types_str['CBD'].epsilon, 0)
        self.assertEqual(p.atom_types_str['CBD'].rmin, 0)
        self.assertAlmostEqual(p.atom_types_str['OND'].epsilon, -0.05)
        self.assertAlmostEqual(p.atom_types_str['OND'].rmin, 2.0)

    def testPenalty(self):
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

    def testCharmmParameterSetConversion(self):
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
                if key != key.upper():
                    return '%sLTU' % key.upper()
                return key
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
            self.assertEqual(set1.cmap_types[typenames(key)], item2)
        # Atom types
        a1, a2 = get_typeset(set1.atom_types, set2.atom_types)
        self.assertEqual(len(a1), len(a2))
        if copy:
            self.assertFalse(a1 & a2)
        else:
            self.assertEqual(a1, a2)

class TestFileWriting(utils.FileIOTestCase):
    """ Tests the various file writing capabilities """

    def testWriteCharmm(self):
        """ Test writing CHARMM-style PSF files """
        # Test writing CHARMM-style PSFs
        cpsf = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf'))
        cpsf.write_psf(get_fn('dhfr_cmap_pbc.psf', written=True))
        cpsf2 = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf', written=True))
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
        f = open(get_fn('dhfr_cmap_pbc.psf', written=True), 'r')
        try:
            has_key = False
            for line in f:
                if '!MOLNT' in line:
                    has_key = True
                    break
        finally:
            f.close()
        self.assertTrue(has_key)

    def testWriteVmd(self):
        """ Test writing VMD-style PSF files """
        # Test writing VMD-style PSFs
        cpsf = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf'))
        cpsf.write_psf(get_fn('dhfr_cmap_pbc.psf', written=True), vmd=True)
        cpsf2 = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf', written=True))
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
        f = open(get_fn('dhfr_cmap_pbc.psf', written=True), 'r')
        try:
            has_key = False
            for line in f:
                if '!MOLNT' in line:
                    has_key = True
                    break
        finally:
            f.close()
        self.assertFalse(has_key)

if __name__ == '__main__':
    unittest.main()
