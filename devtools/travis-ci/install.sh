if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
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
