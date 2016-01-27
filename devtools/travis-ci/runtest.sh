#!/bin/sh
set -e

do_coverage() {
    echo "Combining coverage data"
    coverage combine
    echo "Reporting..."
    coverage report -m
}

echo "Checking parmed source with pyflakes linter"
if [ "$PYTHON_VERSION" = "pypy" ]; then
    export PYENV_ROOT="${HOME}/.pyenv"
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
fi
sh devtools/travis-ci/pyflakes_check.sh
cd test
echo "Using nosetests...:"
./run_scripts.sh
if [ "$PYTHON_VERSION" = "pypy" ]; then
    # Disable coverage with pypy, since it multiplies the time taken by 6 or
    # something ridiculous like that
    nosetests -vs --with-timer --timer-ok=5s --timer-warning=12s \
              --timer-filter=warning,error .
else
    # Run nose under coverage, since that allows getting the full flexibility of
    # the coverage package without sacrificing nose functionality
    coverage run --source=parmed --parallel-mode -m \
        nose -vs --with-timer --timer-ok=5s --timer-warning=12s \
             --timer-filter=warning,error .
fi
test -z `which coverage 2>/dev/null` || do_coverage
echo "Running coveralls"
test -z `which coveralls` || coveralls
echo "Done!"
