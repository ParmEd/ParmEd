#!/bin/bash
set -ex

source activate myenv

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
