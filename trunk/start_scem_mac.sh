#!/bin/bash

export TMPDIR=$HOME/workspace/TMP_SCEM

export MODELDIR=$HOME/workspace/ParSCEM/Model

mpirun   octave --eval master

