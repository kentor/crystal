#ifndef __lattice_h__
#define __lattice_h__

typedef enum state {
   vacuum, silver, surface, pvp
} state_t;

typedef struct site {
   int id;
   int pos[3];
   int nn[12];
   int nn_count;
   int nnn[42];
   int nnn_count;
   int neighbors;
   state_t state;
   double energy;
   double rate;
} site_t;

site_t *new_lattice(int m);

#endif
