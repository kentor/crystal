#include <stdio.h>
#include <stdlib.h>
#include "lattice.h"

static int find_id_by_pos(int _x, int _y, int _z, int m);

site_t *new_lattice(int _m)
{
   int _nsites = 4*_m*_m*_m;
   site_t *site = malloc(_nsites * sizeof(site_t));

   for (int n = 0; n < _nsites; n++) {
      site[n].id = n;
      site[n].neighbors = 0;
      site[n].nn_count = 0;
      site[n].nnn_count = 0;
      site[n].state = vacuum;
      site[n].energy = 0.0;
      site[n].rate = 1.0;
      for (int i = 0; i < 12; site[n].nn[i++] = -1);
      for (int i = 0; i < 42; site[n].nnn[i++] = -1);
   }

   for (int n = 0; n < _nsites; n+=4) {
      int x, y, z, r;
      z = (n/4) / (_m*_m); r = (n/4) % (_m*_m);
      y = r / _m; x = r % _m;

      for (int i = 0; i < 4; i++) {
         site[n+i].pos[0] = 2*x;
         site[n+i].pos[1] = 2*y;
         site[n+i].pos[2] = 2*z;
      }

      site[n+1].pos[0] += 1; site[n+1].pos[1] += 1;
      site[n+2].pos[0] += 1; site[n+2].pos[2] += 1;
      site[n+3].pos[1] += 1; site[n+3].pos[2] += 1;
   }

   for (int n = 0; n < _nsites; n++) {
      int x, y, z;
      x = site[n].pos[0]; y = site[n].pos[1]; z = site[n].pos[2];

      #define make_neighbors_of(id1, _x, _y, _z, _m) \
      do { \
         int id2 = find_id_by_pos((_x), (_y), (_z), (_m)); \
         if (id2 != -1) { \
            site[id1].nn[site[id1].nn_count++] = id2; \
            site[id2].nn[site[id2].nn_count++] = id1; \
         } \
      } while (0)

      #define make_neighbors_neighbors_of(id1, _x, _y, _z, _m) \
      do { \
         int id2 = find_id_by_pos((_x), (_y), (_z), (_m)); \
         if (id2 != -1) { \
            site[id1].nnn[site[id1].nnn_count++] = id2; \
            site[id2].nnn[site[id2].nnn_count++] = id1; \
         } \
      } while (0)

      make_neighbors_of(n, x+1, y+1, z  , _m);
      make_neighbors_of(n, x+1, y-1, z  , _m);
      make_neighbors_of(n, x+1, y  , z+1, _m);
      make_neighbors_of(n, x+1, y  , z-1, _m);
      make_neighbors_of(n, x  , y+1, z+1, _m);
      make_neighbors_of(n, x  , y+1, z-1, _m);

      make_neighbors_neighbors_of(n, x+2, y+2, z  , _m);
      make_neighbors_neighbors_of(n, x+2, y  , z  , _m);
      make_neighbors_neighbors_of(n, x  , y+2, z  , _m);
      make_neighbors_neighbors_of(n, x+2, y+1, z+1, _m);
      make_neighbors_neighbors_of(n, x+2, y+1, z-1, _m);
      make_neighbors_neighbors_of(n, x+1, y+2, z+1, _m);
      make_neighbors_neighbors_of(n, x+1, y+2, z-1, _m);
      make_neighbors_neighbors_of(n, x+2, y-2, z  , _m);
      make_neighbors_neighbors_of(n, x+2, y-1, z+1, _m);
      make_neighbors_neighbors_of(n, x+2, y-1, z-1, _m);
      make_neighbors_neighbors_of(n, x+1, y-2, z-1, _m);
      make_neighbors_neighbors_of(n, x+1, y-2, z+1, _m);
      make_neighbors_neighbors_of(n, x+2, y  , z+2, _m);
      make_neighbors_neighbors_of(n, x  , y  , z+2, _m);
      make_neighbors_neighbors_of(n, x+1, y+1, z+2, _m);
      make_neighbors_neighbors_of(n, x+1, y-1, z+2, _m);
      make_neighbors_neighbors_of(n, x+2, y  , z-2, _m);
      make_neighbors_neighbors_of(n, x+1, y-1, z-2, _m);
      make_neighbors_neighbors_of(n, x+1, y+1, z-2, _m);
      make_neighbors_neighbors_of(n, x  , y+2, z+2, _m);
      make_neighbors_neighbors_of(n, x  , y+2, z-2, _m);
   }

   return site;
}

static int find_id_by_pos(int _x, int _y, int _z, int m)
{
   if (_x < 0 || _x >= 2*m || _y < 0 || _y >= 2*m || _z < 0 || _z >= 2*m)
      return -1;

   int x, y, z, yr, zr;
   x = _x / 2;
   y = _y / 2; yr = _y % 2;
   z = _z / 2; zr = _z % 2;

   return 4*(x + y*m + z*m*m) + (zr ? (yr ? 3 : 2) : (yr ? 1 : 0));
}
