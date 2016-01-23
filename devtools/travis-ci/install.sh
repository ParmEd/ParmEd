if [ "$PYTHON_VERSION" = "pypy" ]; then
    # Upgrade to pypy 4.0.1 -- original recipe taken from google/oauth2client
    git clone https://github.com/yyuu/pyenv.git ${HOME}/.pyenv
    export PYENV_ROOT="${HOME}/.pyenv"
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
    pyenv install pypy-4.0.1
    pyenv global pypy-4.0.1

    pypy -m pip install nose pyflakes nose-timer
    which pyflakes
    pypy -m pip install --user git+https://bitbucket.org/pypy/numpy.git
else # Otherwise, CPython... go through conda
    if [ "$TRAVIS_OS_NAME" = "osx" ]; then
        wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-MacOSX-x86_64.sh -O miniconda.sh;
    else
        wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh;
    fi

    bash miniconda.sh -b

    export PATH=$HOME/miniconda/bin:$PATH
    conda install --yes conda-build jinja2 binstar pip
    conda config --add channels omnia

    if [ -z "$MINIMAL_PACKAGES" ]; then
        conda create -y -n myenv python=$PYTHON_VERSION \
            numpy scipy pandas nose openmm pyflakes coverage nose-timer
        conda update -y -n myenv --all
    else
        # Do not install the full numpy/scipy stack
        conda create -y -n myenv python=$PYTHON_VERSION numpy nose pyflakes \
            coverage nose-timer
    fi

    source activate myenv
fi # CPython
