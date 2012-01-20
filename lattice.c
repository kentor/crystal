#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "parse_arg.h"
#include "set.c"

enum states {
   _vacuum, _silver, _surface
};

typedef struct site {
   int id;
   int pos[3];
   int neighbors;
   enum states state;
   double energy;
   double rate;
   struct site *nn[12];
   int nn_count;
} site;

typedef struct lattice {
   int m;
   site site[];
} lattice;

void lattice_initialize(int argc, char **argv);
int find_id_by_pos(int _x, int _y, int _z, int m);
void run_unit_tests(void);

int _m = 2, _nsites, _seed;
site *lat;
char *_file = "movie.xyz";
bool _draw;

void lattice_initialize(int argc, char **argv)
{
   for (int i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
         char *flag = argv[i] + 1;
         char *arg = argv[i+1];
         parse_arg(flag, "m", _m, atoi(arg));
         parse_arg(flag, "d", _draw, !_draw);
         parse_arg(flag, "o", _file, strdup(arg));
      }
   }

   _nsites = 4*_m*_m*_m;

   lat = malloc(_nsites*sizeof(site));
   for (int id = 0; id < _nsites; id++) {
      lat[id].id = id;
      lat[id].rate = 1.0;
   }

   for (int n = 0; n < _nsites; n += 4) {
      int x, y, z, r;
      z = (n/4) / (_m*_m); r = (n/4) % (_m*_m);
      y = r / _m; x = r % _m;

      for (int i = 0; i < 4; i++) {
         lat[n+i].pos[0] = 2*x;
         lat[n+i].pos[1] = 2*y;
         lat[n+i].pos[2] = 2*z;
      }

      lat[n+1].pos[0] += 1; lat[n+1].pos[1] += 1;
      lat[n+2].pos[0] += 1; lat[n+2].pos[2] += 1;
      lat[n+3].pos[1] += 1; lat[n+3].pos[2] += 1;
   }

   for (int n = 0; n < _nsites; n++) {
      int x, y, z;
      x = lat[n].pos[0]; y = lat[n].pos[1]; z = lat[n].pos[2];

      #ifndef make_neighbors_of
      #define make_neighbors_of(id1, _x, _y, _z, _m) \
      do { \
         int id2 = find_id_by_pos((_x), (_y), (_z), (_m)); \
         if (id2 != -1) { \
            lat[id1].nn[lat[id1].nn_count++] = &lat[id2]; \
            lat[id2].nn[lat[id2].nn_count++] = &lat[id1]; \
         } \
      } while (0)
      #endif

      make_neighbors_of(n, x+1, y+1, z  , _m);
      make_neighbors_of(n, x+1, y-1, z  , _m);
      make_neighbors_of(n, x+1, y  , z+1, _m);
      make_neighbors_of(n, x+1, y  , z-1, _m);
      make_neighbors_of(n, x  , y+1, z+1, _m);
      make_neighbors_of(n, x  , y+1, z-1, _m);
   }
}

int main(int argc, char **argv)
{
   run_unit_tests();
   // lattice_initialize(argc, argv);
}

int find_id_by_pos(int _x, int _y, int _z, int m)
{
   if (_x < 0 || _x >= 2*m || _y < 0 || _y >= 2*m || _z < 0 || _z >= 2*m) {
      return -1;
   }

   int x, y, z, yr, zr;
   x = _x / 2;
   y = _y / 2; yr = _y % 2;
   z = _z / 2; zr = _z % 2;

   return 4*(x + y*m + z*m*m) + (zr ? (yr ? 3 : 2) : (yr ? 1 : 0));
}

void run_unit_tests(void)
{
   describe_find_id_by_pos: {
      it_should_return_id_when_given_pos_and_m: {
         assert(find_id_by_pos(0,0,0,2) == 0);
         assert(find_id_by_pos(1,1,0,2) == 1);
         assert(find_id_by_pos(1,0,1,2) == 2);
         assert(find_id_by_pos(0,1,1,2) == 3);
         assert(find_id_by_pos(2,0,0,2) == 4);
         assert(find_id_by_pos(3,1,0,2) == 5);
         assert(find_id_by_pos(3,0,1,2) == 6);
         assert(find_id_by_pos(2,1,1,2) == 7);
         assert(find_id_by_pos(0,2,0,2) == 8);
         assert(find_id_by_pos(2,2,0,2) == 12);
         assert(find_id_by_pos(0,0,2,2) == 16);
         assert(find_id_by_pos(2,0,2,2) == 20);
         assert(find_id_by_pos(2,2,2,2) == 28);

         assert(find_id_by_pos(0,0,0,30) == 0);
         assert(find_id_by_pos(1,1,0,30) == 1);
         assert(find_id_by_pos(1,0,1,30) == 2);
         assert(find_id_by_pos(0,1,1,30) == 3);
         assert(find_id_by_pos(0,2,0,30)== 120);
         assert(find_id_by_pos(58,58,58,30)== 107996);
         assert(find_id_by_pos(59,59,58,30)== 107997);
         assert(find_id_by_pos(59,58,59,30)== 107998);
         assert(find_id_by_pos(58,59,59,30)== 107999);

         puts(".");
      }

      it_should_return_neg_1_when_given_oob_pos: {
         assert(find_id_by_pos(-1,0,0,2) == -1); 
         assert(find_id_by_pos(4,0,0,2) == -1); 
         assert(find_id_by_pos(0,-1,0,30) == -1); 
         assert(find_id_by_pos(0,0,60,30) == -1); 

         puts(".");
      }
   }
}