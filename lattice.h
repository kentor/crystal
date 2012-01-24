#ifndef LATTICE_H
#define LATTICE_H

typedef enum state {
   vacuum, silver, surface, pvp
} state_t;

typedef struct site {
   int id;
   int pos[3];
   int nn[12];
   int neighbors;
   int nn_count;
   state_t state;
   double energy;
   double rate;
} site_t;

site_t *new_lattice(int m);

#endif
