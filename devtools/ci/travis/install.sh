if [ "$PYTHON_VERSION" = "pypy" ]; then
    # Upgrade to pypy 4.0.1 -- original recipe taken from google/oauth2client
    git clone https://github.com/yyuu/pyenv.git ${HOME}/.pyenv
    export PYENV_ROOT="${HOME}/.pyenv"
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
    pyenv install pypy2.7-6.0.0
    pyenv global pypy2.7-6.0.0

    pypy -m pip install nose pyflakes==1.0.0 nose-timer lxml
    which pyflakes
    pypy -m pip install numpy==1.15.4 # pin 1.15.4, since the latest version doesn't work on pypy2
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
    # Omnia requires conda-forge
    conda config --add channels omnia --add channels conda-forge

    if [ -z "$MINIMAL_PACKAGES" ]; then
        # Install all prerequisites
        conda create -y -n myenv python=$PYTHON_VERSION \
            numpy scipy pandas nose openmm coverage nose-timer \
            netCDF4
        conda update -y -n myenv --all
        conda install -y -n myenv rdkit==2018.09.1
        conda install -y -n myenv boost==1.69.0
        conda install -y -n myenv nglview
        conda install -y -n myenv ambertools=18 -c ambermd
        conda install -y -n myenv networkx
        conda install -y -n myenv lxml
    else
        # Do not install the full numpy/scipy stack
        conda create -y -n myenv python=$PYTHON_VERSION numpy nose coverage \
            nose-timer
    fi
    source activate myenv
    pip install pyflakes==1.0.0
    if [ -z "$MINIMAL_PACKAGES" ]; then
        pip uninstall parmed -y # from ambertools
    fi
    # DEBUG
    conda list
fi # CPython
