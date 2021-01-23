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
py.test --cov=parmed --durations=0 --disable-warnings --cov-append .
echo "Done!"
