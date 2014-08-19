
from distutils.core import setup

# First the ParmedTools packages:
packages = ['ParmedTools', 'ParmedTools.gui', 'ParmedTools.simulations']

# Next the main chemistry packages
packages += ['chemistry', 'chemistry.amber', 'chemistry.tinker',
             'chemistry.amber.mdin', 'cpinutils', 'chemistry.charmm']

# Modules
modules = ['compat24', 'timer']

# Scripts
scripts = ['parmed.py', 'xparmed.py']

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
        fixers = [x for x in get_fixers('lib2to3.fixes') if not 'has_key' in x]
        Mixin2to3.fixer_names = fixers

    setup(name='ParmEd',
          version='SA_1.1.1', # standalone release version
          description='Amber parameter file editor',
          author='Jason Swails',
          author_email='jason.swails -at- gmail.com',
          url='http://jswails.wikidot.com/parmed',
          license='GPL v2 or later',
          packages=packages,
          py_modules=modules,
          cmdclass={'build_py': build_py, 'build_scripts': build_scripts},
          scripts=scripts)
