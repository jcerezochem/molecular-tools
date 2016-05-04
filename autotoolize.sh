#!/bin/bash

rm configure; aclocal; autoconf; automake -a
#rm configure; aclocal; autoconf; sed -i "s/-lmkl_intel/-lmkl_intel_lp64 -lmkl_sequential -lmkl_core/g" configure; automake -a


