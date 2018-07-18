#!/bin/sh
set -ex

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
sh devtools/ci/pyflakes_check.sh
echo "Using ParmEd version `parmed --version`"
cd test
echo "Using nosetests...:"
if [ "$PYTHON_VERSION" = "pypy" ]; then
    # Disable coverage with pypy, since it multiplies the time taken by 6 or
    # something ridiculous like that
    nosetests -vs --with-timer --timer-ok=5s --timer-warning=12s \
              --timer-filter=warning,error .
else
    ./run_scripts.sh
    # Run nose under coverage, since that allows getting the full flexibility of
    # the coverage package without sacrificing nose functionality
#   test -z "$MINIMAL_PACKAGES" && export AMBERHOME=$HOME/miniconda/envs/myenv
    coverage run --source=parmed --parallel-mode -m \
        nose -vs --with-timer --timer-ok=5s --timer-warning=12s \
             --timer-filter=warning,error .
fi
test -z `which coverage 2>/dev/null` || do_coverage
#echo "Running coveralls"
#if [ -z "$MINIMAL_PACKAGES" -a "$PYTHON_VERSION" != 'pypy' ]; then
#    # Only run coveralls on builds that test everything
#    coveralls
#fi
echo "Done!"

# Workaround for annoying Travis failure on Macs
shell_session_update() {
  echo "Faking shell_session_update"
}
