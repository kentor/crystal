#!/bin/sh
cmd='gcc -lm -lgsl -lgslcblas -O2 new_kmc.c new_lattice.c auxilary.c -o kmc.o'
echo $cmd
$cmd
