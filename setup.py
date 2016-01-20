from distutils.core import setup, Extension
import os
import sys

if sys.version_info < (2, 7):
    sys.stderr.write('You must have at least Python 2.7 for ParmEd to work '
                     'correctly.\n')
    sys.exit(0)

is_pypy = '__pypy__' in sys.builtin_module_names

if sys.platform == 'darwin' and not is_pypy:
    # You *need* to use clang and clang++ for ParmEd extensions on a Mac;
    # Anaconda does annoying stuff that breaks this, since their distutils
    # automatically tries to use "gcc", which would conflict with the MacPorts
    # gcc... sigh.
    os.environ['CXX'] = 'clang++'
    os.environ['CC'] = 'clang'

# parmed package and all its subpackages
packages = ['parmed', 'parmed.amber', 'parmed.modeller',
            'parmed.tinker', 'parmed.unit', 'parmed.amber.mdin',
            'parmed.charmm', 'parmed.formats.pdbx', 'parmed.rosetta',
            'parmed.formats', 'parmed.utils.fortranformat', 'parmed.openmm',
            'parmed.utils', 'parmed.gromacs', 'parmed.tools', 'parmed.namd',
            'parmed.tools.gui', 'parmed.tools.simulations']

# Optimized readparm
sources = [os.path.join('src', '_rdparm.cpp'),
           os.path.join('src', 'readparm.cpp')]
depends = [os.path.join('src', 'CompatabilityMacros.h'),
           os.path.join('src', 'readparm.h')]
include_dirs = [os.path.join(os.path.abspath('.'), 'src')]

extensions = [Extension('parmed.amber._rdparm',
                        sources=sources,
                        include_dirs=include_dirs,
                        depends=depends)
]

if __name__ == '__main__':

    from distutils.command.build_py import build_py
    from distutils.command.build_scripts import build_scripts
    import shutil
    import parmed

    # See if we have the Python development headers.  If not, don't build the
    # optimized prmtop parser extension
    from distutils import sysconfig
    if not is_pypy and not os.path.exists(
            os.path.join(sysconfig.get_config_vars()['INCLUDEPY'], 'Python.h')):
        extensions = []

    # Delete old versions with old names of scripts and packages (chemistry and
    # ParmedTools for packages, parmed.py and xparmed.py for scripts)
    def deldir(folder):
        try:
            shutil.rmtree(folder)
        except OSError:
            sys.stderr.write(
                    'Could not remove old package %s; you should make sure\n'
                      'this is completely removed in order to make sure you\n'
                      'do not accidentally use the old version of ParmEd\n' %
                      folder
            )
    def delfile(file):
        try:
            os.unlink(file)
        except OSError:
            sys.stderr.write(
                    'Could not remove old script %s; you should make sure\n'
                      'this is completely removed in order to make sure you\n'
                      'do not accidentally use the old version of ParmEd\n' %
                      file
            )

    for folder in sys.path:
        folder = os.path.realpath(os.path.abspath(folder))
        if folder == os.path.realpath(os.path.abspath('.')): continue
        chem = os.path.join(folder, 'chemistry')
        pmdtools = os.path.join(folder, 'ParmedTools')
        pmd = os.path.join(folder, 'parmed.py')
        xpmd = os.path.join(folder, 'xparmed.py')
        pmdc = os.path.join(folder, 'parmed.pyc')
        xpmdc = os.path.join(folder, 'xparmed.pyc')
        if os.path.isdir(chem): deldir(chem)
        if os.path.isdir(pmdtools): deldir(pmdtools)
        if os.path.exists(pmd): delfile(pmd)
        if os.path.exists(xpmd): delfile(xpmd)
        if os.path.exists(pmdc): delfile(pmdc)
        if os.path.exists(xpmdc): delfile(xpmdc)

    for folder in os.getenv('PATH').split(os.pathsep):
        pmd = os.path.join(folder, 'parmed.py')
        xpmd = os.path.join(folder, 'xparmed.py')
        pmdc = os.path.join(folder, 'parmed.pyc')
        xpmdc = os.path.join(folder, 'xparmed.pyc')
        if os.path.exists(pmd): delfile(pmd)
        if os.path.exists(xpmd): delfile(xpmd)
        if os.path.exists(pmdc): delfile(pmdc)
        if os.path.exists(xpmdc): delfile(xpmdc)

    scripts = [os.path.join('scripts', 'parmed'),
               os.path.join('scripts', 'xparmed')]

    setup(name='ParmEd',
          version=parmed.__version__,
          description='Amber parameter file editor',
          author='Jason Swails',
          author_email='jason.swails -at- gmail.com',
          url='http://jswails.wikidot.com/parmed',
          license='GPL v2 or later',
          packages=packages,
          ext_modules=extensions,
          cmdclass={'build_py' : build_py},
          scripts=scripts,
    )
