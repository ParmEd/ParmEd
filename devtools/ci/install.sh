# Remove ParmEd from AmberTools
rm -f `which parmed` `which xparmed`
py_bin_dir=$(dirname `which python`)
rm -fr ${py_bin_dir}/../lib/python${PYTHON_VERSION}/site-packages/ParmEd*
rm -fr ${py_bin_dir}/../lib/python${PYTHON_VERSION}/site-packages/parmed/
conda list

pip install -e .
python -c "import parmed; print(parmed.__version__)"

do_coverage() {
    echo "Combining coverage data"
    coverage combine
    echo "Reporting..."
    coverage report -m
}

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
fi
do_coverage
echo "Done!"
