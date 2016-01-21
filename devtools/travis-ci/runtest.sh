#!/bin/sh

echo "Checking parmed source with pyflakes linter"
if [ "$PYTHON_VERSION" = "pypy" ]; then
    export PYENV_ROOT="${HOME}/.pyenv"
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
fi
sh devtools/travis-ci/pyflakes_check.sh
cd test
echo "Using nosetests...:"
which nosetests
if [ "$PYTHON_VERSION" = "pypy" ]; then
    # Disable coverage with pyflakes, since it multiplies the time taken by 6 or
    # something ridiculous like that
    nosetests -vs --with-timer --timer-ok=5s --timer-warning=12s \
              --timer-filter=warning,error .
else
    nosetests -vs --with-timer --timer-ok=5s --timer-warning=12s \
              --timer-filter=warning,error --with-coverage \
              --cover-package=parmed .
fi
./run_scripts.sh
coverage report -m
