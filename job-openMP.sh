#!/bin/bash

MAXlevel="12"
tmax="2.5"
We="4.0"
Ohd="0.034"
Bo="0.5"
Ohf="0.670"
hf="0.10"
rhof="1.0"
sigma23="1e0"

export OMP_NUM_THREADS=8
qcc -O2 -Wall -disable-dimensions dropFilm.c -o dropFilm -lm
./dropFilm $MAXlevel $tmax $We $Ohd $Bo $Ohf $hf $rhof $sigma23

