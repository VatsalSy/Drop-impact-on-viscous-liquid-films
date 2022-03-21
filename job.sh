#!/bin/bash

MAXlevel="12"
tmax="2.5"
We="4.0"
Ohd="0.034"
Bo="0.5"
Ohf="0.670"
hf="0.10"

export OMP_NUM_THREADS=8
qcc -fopenmp -Wall -O2 dropFilm.c -o dropFilm -lm
./dropFilm $MAXlevel $tmax $We $Ohd $Bo $Ohf $hf

