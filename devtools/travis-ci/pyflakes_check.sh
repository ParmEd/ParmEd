#!/bin/sh

# This script runs pyflakes and filters out all of the known violations to
# provide an output of what the linters find.

which pyflakes > /dev/null

if [ $? -ne 0 ]; then
    echo "pyflakes not found... cannot run check."
    exit 1
fi

notfound() {
    echo "parmed package directory not found!"
    exit 1
}

test -d parmed || notfound

pyflakes parmed | \
    grep -v -e "from parmed.topologyobjects import \*" \
            -e "from parmed.modeller.residue import \*" \
            -e "from parmed.tools.actions import \*" \
            -e "from parmed.utils.six.moves.tkinter import \*" \
            -e "^parmed\/unit" \
            -e "^parmed\/utils\/six.py" \
            -e "^parmed\/utils\/fortranformat" | tee pyflakes.log

nfail=`cat pyflakes.log | wc -l`

if [ $nfail -gt 0 ]; then
    echo "Detected pyflakes failures"
    exit 1
fi

echo "pyflakes reported clean exit"
/bin/rm -f pyflakes.log
exit 0
