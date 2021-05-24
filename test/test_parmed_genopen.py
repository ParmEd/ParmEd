"""
Tests the generalized open routine that is used for virtually all I/O in ParmEd,
and automatically handles compression, etc.
"""

import bz2
import gzip
from contextlib import closing
from parmed.utils.io import genopen
import unittest
from utils import get_fn, FileIOTestCase, is_jenkins

ALPHABET = 'abcdefghijklmnopqrstuvwxyz\n'

class TestGenopen(FileIOTestCase):

    def test_read_normal(self):
        """ Tests genopen reading a normal text file """
        with closing(genopen(get_fn('4lzt.pdb'))) as f:
            self.assertEqual(f.read(), open(get_fn('4lzt.pdb')).read())
        with closing(genopen(get_fn('4lzt.pdb'), 'r')) as f:
            self.assertEqual(f.read(), open(get_fn('4lzt.pdb')).read())

    def test_write_normal(self):
        """ Tests genopen writing a normal text file """
        with closing(genopen(self.get_fn('tmp.txt', written=True), 'w')) as f:
            f.write(ALPHABET)
        self.assertEqual(open(self.get_fn('tmp.txt', written=True), 'r').read(), ALPHABET)

    def test_read_gzipped(self):
        """ Tests genopen reading a gzipped file """
        with closing(genopen(get_fn('4lzt.pdb.gz'))) as f:
            text = gzip.open(get_fn('4lzt.pdb.gz'), 'r').read()
            self.assertEqual(f.read(), text.decode('ascii'))

    def test_write_gzipped(self):
        """ Tests genopen writing a gzipped file """
        with closing(genopen(self.get_fn('test.gz', written=True), 'w')) as f:
            f.write(ALPHABET)
        text = gzip.open(self.get_fn('test.gz', written=True), 'r').read()
        self.assertEqual(text.decode('ascii'), ALPHABET)

    def test_read_bzipped(self):
        """ Tests genopen reading a bzipped file """
        with closing(genopen(get_fn('4lzt.pdb.bz2'), 'r')) as f:
            text = bz2.BZ2File(get_fn('4lzt.pdb.bz2'), 'r').read()
            self.assertEqual(text.decode('ascii'), f.read())

    def test_write_bzipped(self):
        """ Tests genopen writing a bzipped file """
        with closing(genopen(self.get_fn('test.bz2', written=True), 'w')) as f:
            f.write(ALPHABET)
        text = bz2.BZ2File(self.get_fn('test.bz2', written=True), 'r').read()
        self.assertEqual(text.decode('ascii'), ALPHABET)

    def test_read_normal_URL(self):
        """ Tests genopen reading a remote file """
        url = 'https://github.com/ParmEd/ParmEd/raw/master/test/files/tripos1.mol2'
        with closing(genopen(url, 'r')) as f:
            self.assertEqual(f.read(), open(get_fn('tripos1.mol2')).read())

    def test_read_bzipped_URL(self):
        """ Tests genopen reading a bzipped remote file """
        url = 'https://github.com/ParmEd/ParmEd/raw/master/test/files/4lzt.pdb.bz2'
        with closing(genopen(url, 'r')) as f:
            self.assertEqual(f.read(), genopen(get_fn('4lzt.pdb.bz2')).read())

    def test_read_gzipped_URL(self):
        """ Tests genopen reading a gzipped remote file """
        url = 'https://github.com/ParmEd/ParmEd/raw/master/test/files/4lzt.pdb.gz'
        with closing(genopen(url, 'r')) as f:
            self.assertEqual(f.read(), genopen(get_fn('4lzt.pdb.gz')).read())

    @unittest.skipUnless(is_jenkins(), 'Cannot download files from PDB with Travis')
    def test_read_ftp_URL(self):
        """ Tests genopen reading a ftp remote file """
        url = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/05/205l.cif.gz'
        with closing(genopen(url, 'r')) as f:
            self.assertEqual(f.read(), genopen(get_fn('205l.cif.gz')).read())

    def test_append_normal(self):
        """ Tests genopen appending a normal text file """
        with closing(genopen(self.get_fn('test.txt', written=True), 'w')) as f:
            f.write(ALPHABET)
        with closing(genopen(self.get_fn('test.txt', written=True), 'a')) as f:
            f.write(ALPHABET)
        self.assertEqual(open(self.get_fn('test.txt', written=True)).read(), ALPHABET*2)

    def test_append_gzip(self):
        """ Tests genopen appending a gzipped file """
        with closing(genopen(self.get_fn('test.txt.gz', written=True), 'a')) as f:
            f.write(ALPHABET)
        with closing(genopen(self.get_fn('test.txt.gz', written=True), 'a')) as f:
            f.write(ALPHABET)
        text = gzip.open(self.get_fn('test.txt.gz', written=True)).read()
        self.assertEqual(text.decode('ascii'), ALPHABET*2)

    def test_append_bzip(self):
        """ Tests genopen appending a bzipped file """
        with closing(genopen(self.get_fn('test.txt.bz2', written=True), 'a')) as f:
            f.write(ALPHABET)
        with closing(genopen(self.get_fn('test.txt.bz2', written=True), 'a')) as f:
            f.write(ALPHABET)
        text = bz2.BZ2File(self.get_fn('test.txt.bz2', written=True)).read()
        self.assertEqual(text.decode('ascii'), ALPHABET*2)

    def test_append_remote_file(self):
        """ Tests that genopen appending a remote file fails """
        url = 'http://q4md-forcefieldtools.org/REDDB/projects/W-73/tripos1.mol2'
        self.assertRaises(ValueError, lambda: genopen(url, 'a'))
        self.assertRaises(ValueError, lambda: genopen(url, 'w'))

    def test_read_bad_URL(self):
        """ Tests proper exception handling of non-existent URL """
        self.assertRaises(IOError, lambda: genopen('http://asdkfjasdf.lib'))

