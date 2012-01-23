#ifndef LATTICE_H
#define LATTICE_H

typedef enum state {
   _vacuum, _silver, _surface, _pvp
} state_t;

typedef struct site {
   int id;
   int pos[3];
   int neighbors;
   int nn_count;
   state_t state;
   double energy;
   double rate;
   struct site *nn[12];
} site_t;

site_t *new_lattice(int m);
int find_id_by_pos(int x, int y, int z, int m);

#endif
