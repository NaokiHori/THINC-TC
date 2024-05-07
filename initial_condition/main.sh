#!/bin/bash

## domain
# domain lengths
export lx=0.4e+0
export ly=0.4e+0
export lz=0.8e+0
# number of cell centers
export glisize=32
export gljsize=32
export glksize=64

## volume fraction
export vfrac=0.1

## where to write resulting NPY files
export dirname="output"

python3 main.py ${dirname}
