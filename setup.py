import os
import sys
import versioneer
import struct

if sys.version_info < (2, 7):
    sys.stderr.write('You must have at least Python 2.7 for ParmEd to work '
                     'correctly.\n')
    sys.exit(1)

try:
    if '--no-setuptools' in sys.argv:
        sys.argv.remove('--no-setuptools')
        raise ImportError() # Don't import setuptools...
    from setuptools import setup, Extension
    kws = dict(entry_points={
            'console_scripts' : ['parmed = parmed.scripts:clapp'],
            'gui_scripts' : ['xparmed = parmed.scripts:guiapp']}
    )
except ImportError:
    from distutils.core import setup, Extension
    kws = {'scripts' : [os.path.join('scripts', 'parmed'),
                        os.path.join('scripts', 'xparmed')]
    }

from distutils.command.clean import clean as Clean

def prepare_env_for_osx():
    """ Prepares the environment for OS X building """
    darwin_major_to_osx_map = {
        '11': '10.7',
        '12': '10.8',
        '13': '10.9',
        '14': '10.10',
        '15': '10.11',
        '16': '10.12',
        '17': '10.13',
        '18': '10.14',
    }
    os.environ['CXX'] = 'clang++'
    os.environ['CC'] = 'clang'
    darwin_major = os.uname()[2].split('.')[0]
    if darwin_major in darwin_major_to_osx_map:
        os.environ['MACOSX_DEPLOYMENT_TARGET'] = darwin_major_to_osx_map[darwin_major]

class CleanCommand(Clean):
    """python setup.py clean
    """
    # lightly adapted from scikit-learn package
    description = "Remove build artifacts from the source tree"

    def _clean(self, folder):
        for dirpath, dirnames, filenames in os.walk(folder):
            for filename in filenames:
                if (filename.endswith('.so') or filename.endswith('.pyd')
                        or filename.endswith('.dll')
                        or filename.endswith('.pyc')):
                    os.unlink(os.path.join(dirpath, filename))
            for dirname in dirnames:
                if dirname == '__pycache__':
                    shutil.rmtree(os.path.join(dirpath, dirname))

    def run(self):
        Clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        self._clean('./')


is_pypy = '__pypy__' in sys.builtin_module_names

if sys.platform == 'darwin' and not is_pypy:
    # You *need* to use clang and clang++ for ParmEd extensions on a Mac;
    # Anaconda does annoying stuff that breaks this, since their distutils
    # automatically tries to use "gcc", which would conflict with the MacPorts
    # gcc... sigh.
    prepare_env_for_osx()
elif os.environ.get('CC', '').endswith('pgcc'):
    # PGI compilers don't play nicely with Python extensions. So force GCC
    sys.stderr.write('PGI compilers do not work with Python extensions generally. '
                     'Using GCC instead.\n')
    os.environ['CC'] = 'gcc'
    os.environ['CXX'] = 'g++'

# parmed package and all its subpackages
packages = ['parmed', 'parmed.amber', 'parmed.modeller',
            'parmed.tinker', 'parmed.unit', 'parmed.amber.mdin',
            'parmed.charmm', 'parmed.formats.pdbx', 'parmed.rosetta', 'parmed.rdkit',
            'parmed.formats', 'parmed.utils.fortranformat', 'parmed.openmm',
            'parmed.utils', 'parmed.gromacs', 'parmed.tools', 'parmed.namd',
            'parmed.tools.gui', 'parmed.tools.simulations']

# Optimized readparm
sources = [os.path.join('src', '_rdparm.cpp'),
           os.path.join('src', 'readparm.cpp')]
depends = [os.path.join('src', 'CompatabilityMacros.h'),
           os.path.join('src', 'readparm.h')]
include_dirs = [os.path.join(os.path.abspath('.'), 'src')]

definitions = []

#if using 64 bit python interpreter on Windows, add the MS_WIN64 flag for 64 bit pointers
if sys.platform == 'win32' and (struct.calcsize("P") == 8):
    definitions.append(('MS_WIN64', None))
    definitions.append(('_hypot', 'hypot')) #fix MinGW build -- see http://stackoverflow.com/questions/10660524/error-building-boost-1-49-0-with-gcc-4-7-0/12124708#12124708

extensions = [Extension('parmed.amber._rdparm',
                        sources=sources,
                        include_dirs=include_dirs,
                        depends=depends,
                        define_macros=definitions)
]

if __name__ == '__main__':

    import shutil

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

    cmdclass = dict(clean=CleanCommand)
    cmdclass.update(versioneer.get_cmdclass())
    setup(name='ParmEd',
          version=versioneer.get_version(),
          description='Amber parameter file editor',
          author='Jason Swails',
          author_email='jason.swails@gmail.com',
          url='https://parmed.github.io/ParmEd/html/index.html',
          license='LGPL',
          packages=packages,
          ext_modules=extensions,
          cmdclass=cmdclass,
          test_suite='nose.collector',
          package_data={'parmed.modeller': ['data/*.lib']},
          **kws
    )
