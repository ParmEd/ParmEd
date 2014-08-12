
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
    except ImportError:
        from distutils.command.build_py import build_py

    setup(name='ParmEd',
          version='SA_1.1.1', # standalone release version
          description='Amber parameter file editor',
          author='Jason Swails',
          author_email='jason.swails -at- gmail.com',
          url='http://jswails.wikidot.com/parmed',
          license='GPL v2 or later',
          packages=packages,
          py_modules=modules,
          cmdclass={'build_py': build_py},
          scripts=scripts)
