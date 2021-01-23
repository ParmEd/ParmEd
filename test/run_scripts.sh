#!/bin/sh

if [ -z "`which coverage 2>/dev/null`" ]; then
    run_cmd="python"
    has_coverage="no"
    echo "No coverage module found..."
else
    run_cmd="coverage run -a"
    has_coverage="yes"
    echo "coverage found and will be used..."
fi

old_pwd=`pwd`
cd `dirname $0`
/bin/rm -fr files/writes
/bin/mkdir -p files/writes
errors=0
failures=0

evaluate_test() {
    if [ $1 -ne 0 ]; then
        echo "ERROR"
        errors=`python -c "print($errors + 1)"`
    elif [ $# -eq 2 ]; then
        diff $DIFFARGS files/saved_scripts/$2 files/writes/$2 >> files/writes/DIFFLOG 2>&1
        if [ $? -ne 0 ]; then
            failures=`python -c "print($failures + 1)"`
            echo "FAILED"
        else
            echo "PASSED"
        fi
    elif [ $# -eq 3 ]; then
        diff $DIFFARGS $3 $2 >> files/writes/DIFFLOG 2>&1
        if [ $? -ne 0 ]; then
            failures=`python -c "print($failures + 1)"`
            echo "FAILED"
        else
            echo "PASSED"
        fi
    fi
}
######   TESTS   ######

# Gromacs CPP tests
printf "Running test python -m parmed.gromacs._cpp test 1..."
$run_cmd -m parmed.gromacs._cpp -i files/pptest1/pptest1.h > files/writes/cpptest1
evaluate_test $? cpptest1

printf "Running test python -m parmed.gromacs._cpp test 2..."
$run_cmd -m parmed.gromacs._cpp -i files/pptest1/pptest1.h \
    -Dline -Dpptest1=REPLACED -o files/writes/cpptest2
evaluate_test $? cpptest2

printf "Running test python -m parmed.gromacs._cpp test 3..."
$run_cmd -m parmed.gromacs._cpp -i - < files/pptest1/pptest1.h \
        -Ifiles/pptest1 > files/writes/cpptest1 || exit
evaluate_test $? cpptest1

# parmed CL tests
# TODO: source, parm, ls, cd
DIFFARGS="-I %VERSION -w"
printf "Running test parmed CL test 1..."
$run_cmd `which parmed` -r > files/writes/parmed1.out 2>&1 << EOF
cd files
cd */
cd nodir
cd ash.parm7
parm ash.parm7
outparm writes/parmed_test1.parm7
EOF
evaluate_test $? files/ash.parm7 files/writes/parmed_test1.parm7

printf "Running test parmed CL test 2..."
/bin/rm -fr files/testbed
mkdir files/testbed
touch files/testbed/file1
touch files/testbed/file2
touch files/testbed/file3
mkdir files/testbed/subdir
mkdir files/testbed/subdir2
touch files/testbed/subdir/file1
touch files/testbed/subdir/file2
touch files/testbed/subdir2/file1
touch files/testbed/subdir2/file2
$run_cmd `which parmed` -nr > files/writes/parmed2.out 2>&1 << EOF
cd files/testbed
ls
ls file*
ls nofile
ls */
ls subdir/file?
EOF
if [ $? -ne 0 ]; then
    echo "FAILED"
else
    echo "PASSED"
fi
#evaluate_test $? parmed2.out

printf "Running test parmed CL test 3..."
cat > files/writes/parmed1.in << EOF
cd files/testbed
ls -C
ls file*
ls nofile
ls */
ls subdir/file?
EOF
$run_cmd `which parmed` -nr 2>&1 << EOF
source files/writes/parmed1.in
EOF
if [ $? -ne 0 ]; then
    echo "FAILED"
else
    echo "PASSED"
fi

###### END TESTS ######

# Clean up if everything passed
if [ $failures -eq 0 -a $errors -eq 0 ]; then
    /bin/rm -fr files/writes files/testbed
else
    test -f files/writes/DIFFLOG && cat files/writes/DIFFLOG
fi
cd "$old_pwd"

if [ $failures -gt 0 -o $errors -gt 0 ]; then
    exit 1
fi
exit 0
