
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

   setup(name='ParmEd',
         version='SA_1.0', # standalone release version
         description='Amber parameter file editor',
         author='Jason Swails',
         author_email='jason.swails -at- gmail.com',
         url='http://jswails.wikidot.com/parmed',
         license='GPL v2 or later',
         packages=packages,
         py_modules=modules,
         scripts=scripts)
