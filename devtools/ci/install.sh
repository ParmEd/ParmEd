# Remove ParmEd from AmberTools
rm -f `which parmed` `which xparmed`
py_bin_dir=$(dirname `which python`)
rm -fr ${py_bin_dir}/../lib/python${PYTHON_VERSION}/site-packages/ParmEd*
rm -fr ${py_bin_dir}/../lib/python${PYTHON_VERSION}/site-packages/parmed/
conda list

pip install -I .
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
pwd
./run_scripts.sh
# Run pytest under coverage, since that allows getting the full flexibility of
# the coverage package without sacrificing nose functionality
pwd
py.test --cov=parmed --durations=0 --disable-warnings test
do_coverage
echo "Done!"
