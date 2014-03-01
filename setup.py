
from distutils.core import setup

# First the ParmedTools packages:
packages = ['ParmedTools', 'ParmedTools.gui', 'ParmedTools.simulations']

# Next the main chemistry packages
packages += ['chemistry', 'chemistry.amber', 'chemistry.tinker',
             'chemistry.amber.mdin', 'cpinutils']

# Modules
modules = ['compat24', 'timer']

# Scripts
scripts = ['parmed.py', 'xparmed.py']

if __name__ == '__main__':

   setup(name='ParmEd',
         version='13.5', # Corresponds to Amber release; between Amber 13 and 14
         description='Amber parameter file editor',
         author='Jason Swails',
         author_email='jason.swails -at- gmail.com',
         url='http://jswails.wikidot.com/parmed',
         license='GPL v2 or later',
         packages=packages,
         py_modules=modules,
         scripts=scripts)
