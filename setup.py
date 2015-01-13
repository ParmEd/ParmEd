
from distutils.core import setup, Extension
import os
import sys

# First the ParmedTools packages:
packages = ['ParmedTools', 'ParmedTools.gui', 'ParmedTools.simulations']

# Next the main chemistry packages
packages += ['chemistry', 'chemistry.amber', 'chemistry.tinker',
             'chemistry.unit', 'chemistry.amber.mdin', 'chemistry.charmm',
             'cpinutils', 'fortranformat']

# Modules
modules = ['compat24', 'timer']

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
supports_omm = major > 2 or minor >= 6

def delete_omm_files():
    """ Deletes OpenMM modules in chemistry """
    join = os.path.join
    files = [join('chemistry', 'amber', 'openmmloader.py'),
             join('chemistry', 'charmm', 'openmmloader.py'),
             join('chemistry', 'amber', 'openmmreporters.py'),
             join('chemistry', 'unit', 'doctests.py')]
    ret = dict()
    for f in files:
        ret[f] = open(f, 'r').read()
    for f in files:
        os.unlink(f)
    return ret

def restore_omm_files(filemap):
    """
    Takes the return value of "delete_omm_files" and restores the files to
    their original greatness :)
    """
    # Py3 support not important here, since this code is only executed for
    # Python 2.5 and older
    for f, text in filemap.iteritems():
        open(f, 'w').write(text)

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

    if not supports_omm:
        filemap = delete_omm_files()
        try:
            setup(name='ParmEd',
                  version='SA_1.1.1', # standalone release version
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
        finally:
            restore_omm_files(filemap)
    else:
        setup(name='ParmEd',
              version='SA_1.1.1', # standalone release version
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
