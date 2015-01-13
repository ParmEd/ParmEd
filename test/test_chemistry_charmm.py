"""
Tests for the chemistry/charmm subpackage
"""

from chemistry.charmm import charmmcrds, parameters, psf
from chemistry import topologyobjects as to
from chemistry import exceptions
from compat24 import all
import os
import unittest
import utils

get_fn = utils.get_fn

class TestCharmmCoords(unittest.TestCase):
    """ Test CHARMM coordinate file parsers """
    
    def testCharmmCrd(self):
        """ Test CHARMM coordinate file parser """
        crd = charmmcrds.CharmmCrdFile(get_fn('1tnm.crd'))
        self.assertEqual(crd.natom, 1414)
        self.assertEqual(max(crd.resno), 91)
        self.assertAlmostEqual(sum(crd.coords), -218.19346999999757)
        self.assertEqual(len(crd.coords), 3*crd.natom)
        self.assertEqual(len(crd.atomno), crd.natom)
        self.assertEqual(len(crd.resno), crd.natom)
        self.assertEqual(len(crd.resid), crd.natom)
        self.assertEqual(len(crd.resname), crd.natom)
        self.assertEqual(len(crd.weighting), crd.natom)

    def testCharmmRst(self):
        """ Test CHARMM restart file parser """
        crd = charmmcrds.CharmmRstFile(get_fn('sample-charmm.rst'))
        self.assertEqual(crd.natom, 256)
        self.assertEqual(crd.nstep, 100)
        self.assertTrue(hasattr(crd, 'header'))
        self.assertAlmostEqual(sum(crd.coords), 0.3114525961458884)
        self.assertAlmostEqual(sum(crd.coordsold), 5439.333671681806)
        self.assertAlmostEqual(sum(crd.vels), 42.364377359350534)
        self.assertEqual(len(crd.coords), 3*crd.natom)
        self.assertEqual(len(crd.coordsold), 3*crd.natom)
        self.assertEqual(len(crd.vels), 3*crd.natom)
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
        self.assertRaises(exceptions.MissingParameter, lambda: str(a.atom_type))
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
        self.assertEqual(len(a.residue.chain), 3)
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
        self.assertRaises(exceptions.MissingParameter, lambda: int(a.atom_type))
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
        self.assertEqual(len(a.residue.chain), 3)
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
        self.assertRaises(exceptions.MissingParameter, lambda: int(a.atom_type))
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
        self.assertEqual(len(a.residue.chain), 2)
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

class TestCharmmParameters(unittest.TestCase):
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
        p=parameters.CharmmParameterSet(get_fn('par_all36_carb.prm')).condense()

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

    def testCGenFF(self):
        """ Test parsing stream files generated by CGenFF """
        p = parameters.CharmmParameterSet(get_fn('toppar_spin_label_dummy.str'))
        p = p.condense()
        self.assertEqual(len(p.atom_types_str), 4)
        self.assertEqual(p.atom_types_str['CBD'].epsilon, 0)
        self.assertEqual(p.atom_types_str['CBD'].rmin, 0)
        self.assertAlmostEqual(p.atom_types_str['OND'].epsilon, -0.05)
        self.assertAlmostEqual(p.atom_types_str['OND'].rmin, 2.0)

class TestFileWriting(unittest.TestCase):
    """ Tests the various file writing capabilities """

    def setUp(self):
        try:
            os.makedirs(get_fn('writes'))
        except OSError:
            pass

    def tearDown(self):
        try:
            for f in os.listdir(get_fn('writes')):
                os.unlink(get_fn(f, written=True))
            os.rmdir(get_fn('writes'))
        except OSError:
            pass

    def testWriteCharmm(self):
        """ Test writing CHARMM-style PSF files """
        # Test writing CHARMM-style PSFs
        cpsf = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf'))
        cpsf.write_psf(get_fn('dhfr_cmap_pbc.psf', written=True))
        cpsf2 = psf.CharmmPsfFile(get_fn('dhfr_cmap_pbc.psf', written=True))
        for attr in dir(cpsf):
            if attr.startswith('_'): continue
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
