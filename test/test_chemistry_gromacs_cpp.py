"""
Tests the functionality in the chemistry.gromacs package
"""
from chemistry.gromacs._cpp import CPreProcessor
from chemistry.utils.six.moves import range, zip, StringIO
import os
import unittest
from utils import get_fn, has_numpy

class TestGromacsCpp(unittest.TestCase):
    """ Tests the C-preprocessor written for the Gromacs classes """

    def test_ifdef(self):
        """ Tests CPreProcessor #ifdef/#else/#endif preprocessing """
        f = StringIO('#ifdef MY_DEFINE\nPPVAR is set to MY_DEFINE\n'
                     '#else\nMY_DEFINE is not set\n#endif\n')
        pp = CPreProcessor(f) # no defines
        self.assertEqual(pp.read().strip(), 'MY_DEFINE is not set')
        f.seek(0)
        pp = CPreProcessor(f, defines=dict(MY_DEFINE='SUCCESS'))
        self.assertEqual(pp.read().strip(), 'PPVAR is set to SUCCESS')

    def test_ifndef(self):
        """ Tests CPreProcessor #ifndef/#else#endif preprocessing """
        f = StringIO('#ifndef MY_DEFINE\nMY_DEFINE is not set\n'
                     '#else\nPPVAR is set to MY_DEFINE\n#endif\n')
        pp = CPreProcessor(f) # no defines
        self.assertEqual(pp.read().strip(), 'MY_DEFINE is not set')
        f.seek(0)
        pp = CPreProcessor(f, defines=dict(MY_DEFINE='SUCCESS'))
        self.assertEqual(pp.read().strip(), 'PPVAR is set to SUCCESS')


    def test_define(self):
        """ Tests CPreProcessor #define preprocessing """
        f = StringIO('#ifdef MY_DEFINE\nPPVAR is set to MY_DEFINE\n'
                     '#else\nMY_DEFINE is not set\n#endif\n'
                     '#define MY_DEFINE SUCCESS\n#ifdef MY_DEFINE\n'
                     'PPVAR is set to MY_DEFINE\n#else\nMY_DEFINE is not set\n'
                     '#endif\n')
        pp = CPreProcessor(f)
        self.assertEqual(pp.read().strip(),
                         'MY_DEFINE is not set\nPPVAR is set to SUCCESS')

    def test_undef(self):
        """ Tests CPreProcessor #undef preprocessing """
        f = StringIO("MY_DEFINE\n#define MY_DEFINE Changed\nMY_DEFINE\n"
                     "#undef MY_DEFINE\nMY_DEFINE\n")
        pp = CPreProcessor(f)
        self.assertEqual(pp.read().strip(), "MY_DEFINE\nChanged\nMY_DEFINE")

    def test_multiword_define(self):
        """ Tests CPreProcessor #define with multiple words """
        f = StringIO("#define MY_DEFINE Changed to a sentence\nMY_DEFINE")
        self.assertEqual(CPreProcessor(f).read().strip(),
                         "Changed to a sentence")

    def test_pp_regex(self):
        """ Tests CPreProcessor processing with whitespace """
        f = StringIO("#      define MY_DEFINE Something\nMY_DEFINE")
        self.assertEqual(CPreProcessor(f).read().strip(), "Something")
        f = StringIO("#       ifdef      MY_VAR    \nMY_VAR\n#endif\n"
                     "#\t\tdefine MY_VAR       NoSpaces    \n"
                     "#       ifdef      MY_VAR    \n\"MY_VAR\"\n#endif")
        f.seek(0)
        self.assertEqual(CPreProcessor(f).read().strip(), '"NoSpaces"')

    def test_define_inside_ifdef(self):
        """ Tests CPreProcessor #define inside #ifdef/#ifndef """
        f = StringIO("#ifdef MY_DEFINE\n#define MY_DEFINE "
                     "Something\n#endif\nMY_DEFINE")
        pp = CPreProcessor(f)
        self.assertEqual(pp.read().strip(), "MY_DEFINE")
        f = StringIO("#ifdef MY_DEFINE\n#    define MY_DEFINE Something\n"
                     "#endif\nMY_DEFINE\n#ifndef MY_DEFINE\n"
                     "#    define MY_DEFINE Something special\n#endif\n"
                     "MY_DEFINE")
        pp = CPreProcessor(f)
        self.assertEqual(pp.read().strip(), "MY_DEFINE\nSomething special")

    def test_if1_if0(self):
        """ Tests CPreProcessor use of #if 1 and #if 0 to serve as comments """
        f = StringIO('#if 1\n#ifdef MY_DEFINE\nPPVAR is set to MY_DEFINE\n'
                     '#else\nMY_DEFINE is not set\n#endif\n#endif\n')
        pp = CPreProcessor(f) # no defines
        self.assertEqual(pp.read().strip(), 'MY_DEFINE is not set')
        f.seek(0)
        pp = CPreProcessor(f, defines=dict(MY_DEFINE='SUCCESS'))
        self.assertEqual(pp.read().strip(), 'PPVAR is set to SUCCESS')
        f = StringIO('#if 0\n#ifdef MY_DEFINE\nPPVAR is set to MY_DEFINE\n'
                     '#else\nMY_DEFINE is not set\n#endif\n#endif\n')
        pp = CPreProcessor(f) # no defines
        self.assertFalse(pp.read().strip())
        f.seek(0)
        pp = CPreProcessor(f, defines=dict(MY_DEFINE='SUCCESS'))
        self.assertFalse(pp.read().strip())

    def test_nested_ifs(self):
        """ Tests CPreProcessor nested #ifdef and #if """
        f = StringIO("""
#ifdef MY_VAR
line 1
#   if 1
line 2
#   endif /* 1 */
#else /* ! MY_VAR */
#   if 0
line 3
#   else /* ! 0 */
line 4
#   endif /* 0 */
#endif /* MY_VAR */
""")
        pp = CPreProcessor(f)
        self.assertEqual(pp.read().strip(), "line 4")
        f.seek(0)
        pp = CPreProcessor(f, defines=dict(MY_VAR=1))
        self.assertEqual(pp.read().strip(), "line 1\nline 2")

    def test_simple_include(self):
        """ Tests CPreProcessor simple #include directive """
        pp = CPreProcessor(get_fn('pptest1/pptest1.h'))
        self.assertEqual(pp.read().strip(), """\
pptest1 line 1
pptest2 line 1
pptest3 line 1
pptest2 line 2
pptest1 line 2
pptest3 line 1
pptest1 line 3""")

    def test_include_defines(self):
        """ Tests passing defines between included files """
        pp = CPreProcessor(get_fn('pptest2/pptest1.h'))
        self.assertEqual(pp.read().strip(), """\
pptest1 line 1
pptest2 line 1
pptest3 line 1
pptest1 line 2
pptest1 line 3
pptest3 line 1
pptest1 line 4""")
