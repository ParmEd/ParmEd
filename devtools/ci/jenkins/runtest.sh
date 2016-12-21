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
test -z "$MINIMAL_PACKAGES" && export AMBERHOME=$HOME/miniconda/envs/myenv
coverage run --source=parmed --parallel-mode -m \
    nose -vs --with-timer --timer-ok=5s --timer-warning=12s \
         --timer-filter=warning,error .
test -z `which coverage 2>/dev/null` || do_coverage
echo "Running coveralls"
if [ -z "$MINIMAL_PACKAGES" -a "$PYTHON_VERSION" != 'pypy' ]; then
    # Only run coveralls on builds that test everything
    coveralls
fi
echo "Done!"
