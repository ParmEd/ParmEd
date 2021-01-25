# Remove ParmEd from AmberTools
rm -f `which parmed` `which xparmed`
py_bin_dir=$(dirname `which python`)
rm -fr ${py_bin_dir}/../lib/python${PYTHON_VERSION}/site-packages/ParmEd*
rm -fr ${py_bin_dir}/../lib/python${PYTHON_VERSION}/site-packages/parmed/
conda list

pip install -I .
python -c "import parmed; print(parmed.__version__)"

echo "Using ParmEd version `parmed --version`"
cd test
./run_scripts.sh
py.test -n 4 --cov=parmed --durations-min=5 --disable-warnings --cov-append --cov-report=xml .
echo "Done!"
