kmc.o: new_kmc.c new_lattice.c auxiliary.c new_lattice.h
	gcc -lm -lgsl -lgslcblas -O3 -o kmc.o new_kmc.c new_lattice.c auxiliary.c
