#!/bin/sh

# autotoolize local project
mkdir -p config m4
autoreconf --force --install

# autotoolize third party
(cd htslib && autoreconf)
