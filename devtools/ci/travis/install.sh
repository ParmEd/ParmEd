if [ "$PYTHON_VERSION" = "pypy" ]; then
    # Upgrade to pypy 4.0.1 -- original recipe taken from google/oauth2client
    git clone https://github.com/yyuu/pyenv.git ${HOME}/.pyenv
    export PYENV_ROOT="${HOME}/.pyenv"
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
    pyenv install pypy3.6-7.2.0
    pyenv global pypy3.6-7.2.0

    pypy -m pip install nose nose-timer lxml
    pypy -m pip install numpy
else # Otherwise, CPython... go through conda
    if [ "$TRAVIS_OS_NAME" = "osx" ]; then
        wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-MacOSX-x86_64.sh -O miniconda.sh;
    else
        wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh;
    fi

    bash miniconda.sh -b

    export PATH=$HOME/miniconda/bin:$PATH
    conda update -y --all
    conda install --yes conda-build jinja2 anaconda-client pip
    conda config --add channels conda-forge

    if [ -z "$MINIMAL_PACKAGES" ]; then
        # Install all prerequisites
        conda create -y -c conda-forge -n myenv python=$PYTHON_VERSION \
            pandas nose openmm coverage nose-timer \
            netCDF4 rdkit==2020.09.4 nglview ambertools networkx lxml
    else
        # Do not install the full numpy/scipy stack
        conda create -y -n myenv python=$PYTHON_VERSION numpy nose coverage \
            nose-timer
    fi
    source activate myenv
    if [ -z "$MINIMAL_PACKAGES" ]; then
        pip uninstall parmed -y # from ambertools
    fi
    # DEBUG
    conda list
fi # CPython
