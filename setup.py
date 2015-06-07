from __future__ import print_function

from distutils.core import setup, Extension
import os
import sys

# First the ParmedTools packages:
packages = ['ParmedTools', 'ParmedTools.gui', 'ParmedTools.simulations']

# Next the main parmed packages
packages += ['parmed', 'parmed.amber', 'parmed.modeller',
             'parmed.tinker', 'parmed.unit', 'parmed.amber.mdin',
             'parmed.charmm', 'parmed.formats.pdbx', 'parmed.rosetta',
             'parmed.formats', 'fortranformat', 'parmed.openmm',
             'parmed.utils', 'parmed.gromacs']

# Scripts
scripts = ['xparmed.py']

# Optimized readparm
extensions = [Extension('parmed.amber._rdparm',
                        sources=['src/_rdparm.cpp', 'src/readparm.cpp'],
                        include_dirs=[os.path.join(os.path.abspath('.'),'src')],
                        depends=['src/CompatibilityMacros.h', 'src/readparm.h']
)]

if __name__ == '__main__':

    from distutils.command.build_py import build_py
    from distutils.command.build_scripts import build_scripts
    import shutil

    # See if we have the Python development headers.  If not, don't build the
    # optimized prmtop parser extension
    from distutils import sysconfig
    if not os.path.exists(os.path.join(sysconfig.get_config_vars()['INCLUDEPY'],
                                       'Python.h')):
        extensions = []

    # Since we changed package names from "chemistry" to "parmed", make sure we
    # delete all of the old versions
    for folder in sys.path:
        folder = os.path.realpath(os.path.abspath(folder))
        if folder == os.path.realpath(os.path.abspath('.')): continue
        chem = os.path.join(folder, 'chemistry')
        pmdtools = os.path.join(folder, 'ParmedTools')
        if os.path.isdir(chem):
            try:
                shutil.rmtree(chem)
                print('Removing %s' % chem, file=sys.stderr)
            except OSError:
                print('Could not remove old chemistry package %s; you should\n'
                      'make sure this is completely removed in order to make\n'
                      'sure you do not accidentally use the old version of '
                      'ParmEd' % chem, file=sys.stderr)
        if os.path.isdir(pmdtools):
            try:
                shutil.rmtree(pmdtools)
                print('Removing %s' % pmdtools, file=sys.stderr)
            except OSError:
                print('Could not remove old ParmedTools package %s; you should\n'
                      'make sure this is completely removed in order to make\n'
                      'sure you do not accidentally use the old version of '
                      'ParmEd' % pmdtools, file=sys.stderr)

    setup(name='ParmEd',
          version='15.0b',
          description='Amber parameter file editor',
          author='Jason Swails',
          author_email='jason.swails -at- gmail.com',
          url='http://jswails.wikidot.com/parmed',
          license='GPL v2 or later',
          packages=packages,
          ext_modules=extensions,
          cmdclass={'build_py':build_py, 'build_scripts':build_scripts},
          scripts=scripts)
