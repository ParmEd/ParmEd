if [ "$PYTHON_VERSION" = "pypy" ]; then
    # Upgrade to pypy 2.6 -- recipe taken from google/oauth2client
    git clone https://github.com/yyuu/pyenv.git ${HOME}/.pyenv
    export PYENV_ROOT="${HOME}/.pyenv"
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
    pyenv install pypy-2.6.0
    pyenv global pypy-2.6.0

    pypy -m pip install --user git+https://bitbucket.org/pypy/numpy.git
    pypy -m pip install nose coverage
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
            $NUMPYVER scipy netcdf4 pandas nose openmm pyflakes
        conda update -y -n myenv --all
    else
        # Do not install the full numpy/scipy stack
        conda create -y -n myenv python=$PYTHON_VERSION $NUMPYVER nose pyflakes
    fi

    source activate myenv

    if [ ! -z "$SCIENTIFIC" ]; then
        # Install ScientificPython
        wget https://sourcesup.renater.fr/frs/download.php/file/4570/ScientificPython-2.9.4.tar.gz
        tar zxvf ScientificPython-2.9.4.tar.gz
        cd ScientificPython-2.9.4
        binprefix=`dirname \`which python\``
        python setup.py install --netcdf_prefix=`dirname $binprefix`
        cd ../
    fi
fi # CPython
