#ifndef LATTICE_H
#define LATTICE_H

#define _site(ptr) ((site_t *)(ptr))

typedef enum state {
   _vacuum, _silver, _surface, _pvp
} state;

typedef struct site {
   int id;
   int pos[3];
   int neighbors;
   int nn_count;
   state state;
   double energy;
   double rate;
   struct site *nn[12];
} site_t;

typedef struct lattice {
   int m;
   int nsites;
   site_t *site;
} lattice_t;

lattice_t new_lattice(int m);
int find_id_by_pos(int x, int y, int z, int m);

#endif
