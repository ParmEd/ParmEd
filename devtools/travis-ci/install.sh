sudo apt-get update
sudo apt-get install -qq -y g++ gfortran valgrind csh

MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b


export PATH=$HOME/miniconda/bin:$PATH
conda install --yes conda-build jinja2 binstar pip
conda config --add channels http://conda.binstar.org/omnia

# Python 2.6 does not have argparse or an OpenMM conda package
openmm="openmm"
if [ "$PYTHON_VERSION" = "2.6" ]; then
    argparse="argparse"
    openmm=""
fi

if [ -z "$NO_NUMPY" ]; then
    conda create -y -n myenv python=$PYTHON_VERSION $argparse \
        numpy scipy netcdf4 pandas nose $openmm
else
    # Do not install the full numpy/scipy stack
    conda create -y -n myenv python=$PYTHON_VERSION $argparse nose
fi

source activate myenv
