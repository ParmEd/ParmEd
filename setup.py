
from distutils.core import setup, Extension
import os
import sys

# First the ParmedTools packages:
packages = ['ParmedTools', 'ParmedTools.gui', 'ParmedTools.simulations']

# Next the main chemistry packages
packages += ['chemistry', 'chemistry.amber', 'chemistry.modeller',
             'chemistry.tinker', 'chemistry.unit', 'chemistry.amber.mdin',
             'chemistry.charmm', 'chemistry.formats.pdbx', 'chemistry.rosetta',
             'chemistry.formats', 'fortranformat', 'chemistry.openmm',
             'chemistry.utils', 'chemistry.gromacs']

# Scripts
scripts = ['parmed.py', 'xparmed.py']

# Optimized readparm
extensions = [Extension('chemistry.amber._rdparm',
                        sources=['src/_rdparm.cpp', 'src/readparm.cpp'],
                        include_dirs=[os.path.join(os.path.abspath('.'),'src')],
                        depends=['src/CompatibilityMacros.h', 'src/readparm.h']
)]

if __name__ == '__main__':

    from distutils.command.build_py import build_py
    from distutils.command.build_scripts import build_scripts

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
          ext_modules=extensions,
          cmdclass={'build_py':build_py, 'build_scripts':build_scripts},
          scripts=scripts)
