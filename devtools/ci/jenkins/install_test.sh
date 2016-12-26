#!/bin/bash
# Make sure we upgrade everything
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

# Now install
python setup.py install

# Now run the tests
do_coverage() {
    echo "Combining coverage data"
    coverage combine
    echo "Reporting..."
    coverage report -m
}

bash devtools/ci/pyflakes_check.sh
echo "Using ParmEd version `parmed --version`"
cd test
echo "Using nosetests...:"
./run_scripts.sh
# Run nose under coverage, since that allows getting the full flexibility of
# the coverage package without sacrificing nose functionality
amber_bin=`dirname \`which tleap\``
export AMBERHOME=`dirname "$amber_bin"`
coverage run --source=parmed --parallel-mode -m \
    nose -vs --with-timer --timer-ok=5s --timer-warning=12s \
         --timer-filter=warning,error .
do_coverage
echo "Done!"
