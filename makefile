kmc: kmc.c lattice.c lattice.h ll.h parse_arg.h
	gcc -o kmc -lm -lgsl -lgslcblas -O3 -std=gnu99 kmc.c lattice.c