
from distutils.core import setup, Extension
import os
import sys

# First the ParmedTools packages:
packages = ['ParmedTools', 'ParmedTools.gui', 'ParmedTools.simulations']

# Next the main chemistry packages
packages += ['chemistry', 'chemistry.amber', 'chemistry.modeller',
             'chemistry.tinker', 'chemistry.unit', 'chemistry.amber.mdin',
             'chemistry.charmm', 'chemistry.formats.pdbx', 'chemistry.formats',
             'fortranformat']

# Modules
modules = ['compat24']

# Scripts
scripts = ['parmed.py', 'xparmed.py']

# Optimized readparm
extensions = [Extension('chemistry.amber._rdparm',
                        sources=['src/_rdparm.cpp', 'src/readparm.cpp'],
                        include_dirs=[os.path.join(os.path.abspath('.'),'src')],
                        depends=['src/CompatibilityMacros.h', 'src/readparm.h']
)]

# See if our Python version will support OpenMM. Of the ParmEd-supported
# Pythons, only 2.4 and 2.5 do not work with OpenMM
major, minor = sys.version_info[:2]
if major > 2 or minor >= 6:
    packages.append('chemistry.openmm')

if __name__ == '__main__':

    try:
        from distutils.command.build_py import build_py_2to3 as build_py
        from distutils.command.build_scripts import build_scripts_2to3 as build_scripts
        PY3 = True
    except ImportError:
        from distutils.command.build_py import build_py
        from distutils.command.build_scripts import build_scripts
        PY3 = False

    if PY3:
        # Monkey-patch Mixin2to3 to disable the has_key fix. ParmEd doesn't ever
        # use "has_key" on dicts, and ArgumentList defines that method. So we
        # want to preserve has_key everywhere it appears.
        from distutils.util import Mixin2to3
        from lib2to3.refactor import get_fixers_from_package as get_fixers
        fixers = [x for x in get_fixers('lib2to3.fixes')
                    if not 'has_key' in x and not 'itertools_imports' in x]
        Mixin2to3.fixer_names = fixers

    # See if we have the Python development headers.  If not, don't build the
    # optimized prmtop parser extension
    from distutils import sysconfig
    if not os.path.exists(os.path.join(sysconfig.get_config_vars()['INCLUDEPY'],
                                       'Python.h')):
        extensions = []

    setup(name='ParmEd',
          version='15.0b',
          description='Amber parameter file editor',
          author='Jason Swails',
          author_email='jason.swails -at- gmail.com',
          url='http://jswails.wikidot.com/parmed',
          license='GPL v2 or later',
          packages=packages,
          py_modules=modules,
          ext_modules=extensions,
          cmdclass={'build_py':build_py, 'build_scripts':build_scripts},
          scripts=scripts)
