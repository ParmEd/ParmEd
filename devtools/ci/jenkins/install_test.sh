#!/bin/bash
CONDAENV=jenkins-parmed-${PYTHON_VERSION}

# First remove the existing conda environment
conda remove -yn ${CONDAENV} --all || true

# Now create the conda environment
conda create -n ${CONDAENV} --no-default-packages
conda config --add channels omnia
conda config --add channels ambermd

# Now add the packages we want
conda install -yn ${CONDAENV} numpy scipy pandas nose openmm coverage nose-timer \
                              python-coveralls ambermini=16.16 netCDF4
conda install -yn ${CONDAENV} pyflakes=1.0.0
conda install -yn ${CONDAENV} rdkit==2015.09.1 -c omnia
conda install -yn ${CONDAENV} boost==1.59.0 -c omnia
conda install -yn ${CONDAENV} nglview -c bioconda

# Make sure we don't install pysander prereqs, since that is just ParmEd!
conda install -yn ${CONDAENV} --no-deps pysander

# Now enter this superamazingawesome environment we just created
source activate ${CONDAENV}

# Lint
echo "Checking the parmed source code with pyflakes"
sh devtools/ci/pyflakes_check.sh

# Now install ParmEd
python setup.py install

# Now run the tests
cd test

export AMBERHOME="`dirname \`which python\``/.."
./run_scripts.sh
coverage run --source=parmed --parallel-mode -m \
    nose -v --with-timer --timer-ok=5s --timer-warning=12s --timer-filter=warning,error .

coverage combine
coverage report -m
coveralls
