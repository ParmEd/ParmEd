#!/bin/sh

# This script processes the HTML files to remove references to _directory and
# replace it with directory... Github-pages seems to require this, but I'm not
# sure how to get it to behave this way automatically from the beginning, so
# this is a cleanser

for html in `find html/ -name "*.html"`; do
    sed -i -e "s/_modules/modules/g" \
           -e "s/_sources/sources/g" \
           -e "s/_static/static/g" \
        $html
done
