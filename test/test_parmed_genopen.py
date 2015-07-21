"""
Tests the generalized open routine that is used for virtually all I/O in ParmEd,
and automatically handles compression, etc.
"""

import bz2
import gzip
from contextlib import closing
from parmed.utils.io import genopen
import unittest
from utils import get_fn, FileIOTestCase

ALPHABET = 'abcdefghijklmnopqrstuvwxyz\n'

class TestGenopen(FileIOTestCase):

    def testReadNormal(self):
        """ Tests genopen reading a normal text file """
        with closing(genopen(get_fn('4lzt.pdb'))) as f:
            self.assertEqual(f.read(), open(get_fn('4lzt.pdb')).read())
        with closing(genopen(get_fn('4lzt.pdb'), 'r')) as f:
            self.assertEqual(f.read(), open(get_fn('4lzt.pdb')).read())

    def testWriteNormal(self):
        """ Tests genopen writing a normal text file """
        with closing(genopen(get_fn('tmp.txt', written=True), 'w')) as f:
            f.write(ALPHABET)
        self.assertEqual(open(get_fn('tmp.txt', written=True), 'r').read(),
                         ALPHABET)

    def testReadGzipped(self):
        """ Tests genopen reading a gzipped file """
        with closing(genopen(get_fn('4lzt.pdb.gz'))) as f:
            text = gzip.open(get_fn('4lzt.pdb.gz'), 'r').read()
            self.assertEqual(f.read(), text.decode('ascii'))

    def testWriteGzipped(self):
        """ Tests genopen writing a gzipped file """
        with closing(genopen(get_fn('test.gz', written=True), 'w')) as f:
            f.write(ALPHABET)
        text = gzip.open(get_fn('test.gz', written=True), 'r').read()
        self.assertEqual(text.decode('ascii'), ALPHABET)

    def testReadBzipped(self):
        """ Tests genopen reading a bzipped file """
        with closing(genopen(get_fn('4lzt.pdb.bz2'), 'r')) as f:
            text = bz2.BZ2File(get_fn('4lzt.pdb.bz2'), 'r').read()
            self.assertEqual(text.decode('ascii'), f.read())

    def testWriteBzipped(self):
        """ Tests genopen writing a bzipped file """
        with closing(genopen(get_fn('test.bz2', written=True), 'w')) as f:
            f.write(ALPHABET)
        text = bz2.BZ2File(get_fn('test.bz2', written=True), 'r').read()
        self.assertEqual(text.decode('ascii'), ALPHABET)

    def testReadNormalURL(self):
        """ Tests genopen reading a remote file """
        url = 'http://q4md-forcefieldtools.org/REDDB/projects/W-73/tripos1.mol2'
        with closing(genopen(url, 'r')) as f:
            self.assertEqual(f.read(), open(get_fn('tripos1.mol2')).read())

    def testReadBzippedURL(self):
        """ Tests genopen reading a bzipped remote file """
        url = 'https://github.com/ParmEd/ParmEd/raw/master/test/files/4lzt.pdb.bz2'
        with closing(genopen(url, 'r')) as f:
            self.assertEqual(f.read(), genopen(get_fn('4lzt.pdb.bz2')).read())

    def testReadGzippedURL(self):
        """ Tests genopen reading a gzipped remote file """
        url = 'https://github.com/ParmEd/ParmEd/raw/master/test/files/4lzt.pdb.gz'
        with closing(genopen(url, 'r')) as f:
            self.assertEqual(f.read(), genopen(get_fn('4lzt.pdb.gz')).read())
