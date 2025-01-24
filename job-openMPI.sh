#!/bin/bash

# Source project configuration if it exists
if [ -f ".project_config" ]; then
    source .project_config
fi

# Check if filename argument is provided
if [ -z "$1" ]; then
    echo "Error: No filename provided"
    echo "Usage: $0 <filename>"
    exit 1
fi

MAXlevel="12"
tmax="2.5"
We="4.0"
Ohd="0.034"
Bo="0.5"
Ohf="0.670"
hf="0.10"
rhof="1.0"
sigma23="1e0"

CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions $1.c -o $1 -lm
mpirun -np 4 ./$1 $MAXlevel $tmax $We $Ohd $Bo $Ohf $hf $rhof $sigma23
