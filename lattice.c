#include <stdio.h>
#include <stdlib.h>
#include "lattice.h"

lattice_t new_lattice(int _m)
{
   int _nsites = 4*_m*_m*_m;
   lattice_t lat;
   lat.m = _m;
   lat.nsites = _nsites;
   lat.site = malloc(_nsites*sizeof(site_t));

   for (int n = 0; n < _nsites; n++) {
      lat.site[n].id = n;
      lat.site[n].neighbors = 0;
      lat.site[n].nn_count = 0;
      lat.site[n].state = _vacuum;
      lat.site[n].energy = 0.0;
      lat.site[n].rate = 1.0;
      for (int i = 0; i < 12; lat.site[n].nn[i++] = NULL);
   }

   for (int n = 0; n < _nsites; n += 4) {
      int x, y, z, r;
      z = (n/4) / (_m*_m); r = (n/4) % (_m*_m);
      y = r / _m; x = r % _m;

      for (int i = 0; i < 4; i++) {
         lat.site[n+i].pos[0] = 2*x;
         lat.site[n+i].pos[1] = 2*y;
         lat.site[n+i].pos[2] = 2*z;
      }

      lat.site[n+1].pos[0] += 1; lat.site[n+1].pos[1] += 1;
      lat.site[n+2].pos[0] += 1; lat.site[n+2].pos[2] += 1;
      lat.site[n+3].pos[1] += 1; lat.site[n+3].pos[2] += 1;
   }

   for (int n = 0; n < _nsites; n++) {
      int x, y, z;
      x = lat.site[n].pos[0]; y = lat.site[n].pos[1]; z = lat.site[n].pos[2];

      #ifndef make_neighbors_of
      #define make_neighbors_of(id1, _x, _y, _z, _m) \
      do { \
         int id2 = find_id_by_pos((_x), (_y), (_z), (_m)); \
         if (id2 != -1) { \
            lat.site[id1].nn[lat.site[id1].nn_count++] = &lat.site[id2]; \
            lat.site[id2].nn[lat.site[id2].nn_count++] = &lat.site[id1]; \
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

   return lat;
}

int find_id_by_pos(int _x, int _y, int _z, int m)
{
   if (_x < 0 || _x >= 2*m || _y < 0 || _y >= 2*m || _z < 0 || _z >= 2*m)
      return -1;

   int x, y, z, yr, zr;
   x = _x / 2;
   y = _y / 2; yr = _y % 2;
   z = _z / 2; zr = _z % 2;

   return 4*(x + y*m + z*m*m) + (zr ? (yr ? 3 : 2) : (yr ? 1 : 0));
}
