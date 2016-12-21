#!/bin/bash
conda update conda -y
conda install --yes conda-build jinja2 binstar pip
conda config --add channels omnia
conda config --add channels ambermd

conda create -y -n myenv python=$PYTHON_VERSION \
    numpy scipy pandas nose openmm coverage nose-timer \
    python-coveralls ambermini netCDF4
conda update -y -n myenv --all
conda install -y -n myenv pyflakes=1.0.0
conda install -y -n myenv rdkit==2015.09.1 -c omnia
conda install -y -n myenv boost==1.59.0 -c omnia
conda install -y -n myenv nglview -c bioconda
# Do not build dependencies for pysander, since that will pull in an old
# version of ParmEd.
conda install -y --no-deps pysander

source activate myenv

python setup.py install
