"""
Tests the chemistry/structure module
"""
import cStringIO as StringIO
import chemistry.structure as structure
import unittest
from utils import get_fn

class TestChemistryStructure(unittest.TestCase):
    
    def setUp(self):
        self.pdb = get_fn('4lzt.pdb')
        self.pdbgz = get_fn('4lzt.pdb.gz')
        self.pdbbz2 = get_fn('4lzt.pdb.bz2')
        self.models = get_fn('2koc.pdb')
        self.overflow = get_fn('4lyt_vmd.pdb')
        self.simple = get_fn('ala_ala_ala.pdb')

    def testAscii(self):
        """ Test PDB file parsing """
        self._check4lyt(structure.read_PDB(self.pdb))
        # The PDB file with multiple models
        pdbfile = structure.read_PDB(open(self.models))
        self.assertEqual(len(pdbfile.pdbxyz), 20)
        self.assertEqual(pdbfile.pdbxyz[0][:3], [-8.886, -5.163, 9.647])
        self.assertEqual(pdbfile.pdbxyz[19][-3:], [-12.051, 5.205, -2.146])

    def testGzip(self):
        """ Test Gzipped-PDB file parsing """
        self._check4lyt(structure.read_PDB(self.pdbgz))

    def testBzip(self):
        """ Test Bzipped-PDB file parsing """
        self._check4lyt(structure.read_PDB(self.pdbbz2))

    def testVmdOverflow(self):
        """ Test PDB file where atom and residue numbers overflow """
        pdbfile = structure.read_PDB(self.overflow)
        self.assertEqual(len(pdbfile.atoms), 110237)
        self.assertEqual(len(pdbfile.residues), 35697)
        self.assertEqual(pdbfile.box, [0, 0, 0, 90, 90, 90])

    def testPdbWriteSimple(self):
        """ Test PDB file writing on a very simple input structure """
        pdbfile = structure.read_PDB(self.simple)
        self.assertEqual(len(pdbfile.atoms), 33)
        self.assertEqual(len(pdbfile.residues), 3)
        output = StringIO.StringIO()
        pdbfile.write_pdb(output)
        output.seek(0)
        pdbfile2 = structure.read_PDB(output)
        self.assertEqual(len(pdbfile2.atoms), 33)
        self.assertEqual(len(pdbfile2.residues), 3)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)

    def testPdbWriteModels(self):
        """ Test PDB file writing from NMR structure with models """
        pdbfile = structure.read_PDB(self.models)
        self.assertEqual(len(pdbfile.pdbxyz), 20)
        self.assertEqual(len(pdbfile.atoms), 451)
        output = StringIO.StringIO()
        structure.write_PDB(pdbfile, output)
        output.seek(0)
        pdbfile2 = structure.read_PDB(output)
        self.assertEqual(len(pdbfile2.atoms), 451)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)

    def testPdbWriteXtal(self):
        """ Test PDB file writing from a Xtal structure """
        pdbfile = structure.read_PDB(self.pdb)
        self._check4lyt(pdbfile)
        output = StringIO.StringIO()
        pdbfile.write_pdb(output, renumber=False)
        output.seek(0)
        pdbfile2 = structure.read_PDB(output)
        self._check4lyt(pdbfile)
        self._compareInputOutputPDBs(pdbfile, pdbfile2)

    def _compareInputOutputPDBs(self, pdbfile, pdbfile2):
        # Now go through all atoms and compare their attributes
        for a1, a2 in zip(pdbfile.atoms, pdbfile2.atoms):
            self.assertEqual(a1.atomic_number, a2.atomic_number)
            self.assertEqual(a1.name, a2.name)
            self.assertEqual(a1.type, a2.type)
            self.assertEqual(a1.mass, a2.mass)
            self.assertEqual(a1.charge, a2.charge)
            self.assertEqual(a1.occupancy, a2.occupancy)
            self.assertEqual(a1.bfactor, a2.bfactor)
            self.assertEqual(a1.altloc, a2.altloc)
            self.assertEqual(a1.idx, a2.idx)
            self.assertEqual(set(a1.other_locations.keys()),
                             set(a2.other_locations.keys()))
            self.assertEqual(a1.number, a2.number)
            self.assertEqual(a1.xx, a2.xx)
            self.assertEqual(a1.xy, a2.xy)
            self.assertEqual(a1.xz, a2.xz)
            # Search all alternate locations as well
            for k1, k2 in zip(sorted(a1.other_locations.keys()),
                              sorted(a2.other_locations.keys())):
                self.assertEqual(k1, k2)
                oa1 = a1.other_locations[k1]
                oa2 = a2.other_locations[k2]
                self.assertEqual(oa1.atomic_number, oa2.atomic_number)
                self.assertEqual(oa1.name, oa2.name)
                self.assertEqual(oa1.type, oa2.type)
                self.assertEqual(oa1.mass, oa2.mass)
                self.assertEqual(oa1.charge, oa2.charge)
                self.assertEqual(oa1.occupancy, oa2.occupancy)
                self.assertEqual(oa1.bfactor, oa2.bfactor)
                self.assertEqual(oa1.altloc, oa2.altloc)
                self.assertEqual(oa1.idx, oa2.idx)
        # Now compare all residues
        for r1, r2 in zip(pdbfile.residues, pdbfile2.residues):
            self.assertEqual(r1.name, r2.name)
            self.assertEqual(r1.idx, r2.idx)
            self.assertEqual(r1.ter, r2.ter)
            self.assertEqual(len(r1), len(r2))
            self.assertEqual(r1.insertion_code, r2.insertion_code)
            self.assertEqual(r1.number, r2.number)

    # Private helper test functions
    def _check4lyt(self, obj):
        self.assertEqual(len(obj.pdbxyz), 1)
        self.assertEqual(obj.box,
                         [27.24, 31.87, 34.23, 88.52, 108.53, 111.89])
        self.assertEqual(obj.space_group, 'P 1')
        self.assertEqual(len(obj.atoms), 1164)
        self.assertEqual(len(obj.residues[0]), 9)
        # Check that alternate conformations are taken into account
        total_natoms = 0
        for i, atom in enumerate(obj.atoms):
            total_natoms += 1
            for key in atom.other_locations:
                total_natoms += 1
                atom2 = atom.other_locations[key]
                self.assertEqual(atom.altloc, 'A')
                self.assertEqual(atom2.altloc, 'B')
                if i in [388, 389]:
                    # Something wrong with atoms 388 and 389 -- sum of their
                    # occupancies are 1.02 :)
                    self.assertEqual(atom2.occupancy + atom.occupancy, 1.02)
                else:
                    self.assertEqual(atom2.occupancy + atom.occupancy, 1)
        self.assertEqual(total_natoms, 1183)
        # Check the metadata
        self.assertEqual(obj.experimental, 'X-RAY DIFFRACTION')
        self.assertEqual(len(obj.residues), 274)
        self.assertEqual(obj.pmid, '9761848')
        self.assertEqual(obj.journal_authors, 'M.A.WALSH,T.R.SCHNEIDER,'
                         'L.C.SIEKER,Z.DAUTER,V.S.LAMZIN,K.S.WILSON')
        self.assertEqual(obj.journal, 'ACTA CRYSTALLOGR.,SECT.D')
        self.assertEqual(obj.year, 1998)
        self.assertEqual(obj.keywords, ['HYDROLASE', 'O-GLYCOSYL',
                         'GLYCOSIDASE'])
        self.assertEqual(obj.title, 'REFINEMENT OF TRICLINIC HEN EGG-WHITE '
                         'LYSOZYME AT ATOMIC RESOLUTION.')
        self.assertEqual(obj.doi, '10.1107/S0907444997013656')
        self.assertEqual(obj.volume, '54')
        self.assertEqual(obj.page, '522')
        # Check the TER card is picked up
        for i, residue in enumerate(obj.residues):
            if i == 128:
                self.assertTrue(residue.ter)
            else:
                self.assertFalse(residue.ter)
