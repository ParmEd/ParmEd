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
pytest \
    -n 4 \
    --cov-branch \
    --cov=parmed \
    --cov-append \
    --cov-report=xml \
    --durations-min=5 \
    -W ignore::UserWarning \
    -W ignore::parmed.exceptions.ParmedWarning \
    .

echo "Done!"
