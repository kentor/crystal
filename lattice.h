typedef enum state {
   _vacuum, _silver, _surface
} state;

typedef struct site {
   int id;
   int pos[3];
   int neighbors;
   state state;
   double energy;
   double rate;
   int nn_count;
   struct site *nn[12];
} site;

typedef struct lattice {
   int m;
   int nsites;
   site *site;
} lattice;

lattice lattice_new(int _m);
int find_id_by_pos(int _x, int _y, int _z, int m);
