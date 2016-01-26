#!/bin/sh

if [ -z "`which coverage 2>/dev/null`" ]; then
    run_cmd="python"
    has_coverage="no"
    echo "No coverage module found..."
else
    run_cmd="coverage run --parallel-mode --source=parmed"
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
    else
        diff files/saved_scripts/$2 files/writes/$2 > /dev/null 2>&1
        if [ $? -ne 0 ]; then
            failures=`python -c "print($failures + 1)"`
            echo "FAILED"
        else
            echo "PASSED"
        fi
    fi
}
######   TESTS   ######
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
###### END TESTS ######

# Clean up if everything passed
if [ $failures -eq 0 ]; then
    /bin/rm -fr files/writes
fi
cd "$old_pwd"

if [ $failures -gt 0 -o $errors -gt 0 ]; then
    exit 1
fi
exit 0
