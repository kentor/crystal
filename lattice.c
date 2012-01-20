#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
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
   int nn_count;
   struct site *nn[12];
} site;

typedef struct lattice {
   int m;
   int nsites;
   site site[];
} lattice;

lattice *lattice_new(int _m);
void lattice_free(lattice *lat);
int find_id_by_pos(int _x, int _y, int _z, int m);
void run_unit_tests(void);

lattice *lattice_new(int _m)
{
   int _nsites = 4*_m*_m*_m;
   lattice *lat = malloc(sizeof(lattice) + _nsites*sizeof(site));
   lat->m = _m;
   lat->nsites = _nsites;

   for (int id = 0; id < _nsites; id++) {
      lat->site[id].id = id;
      lat->site[id].rate = 1.0;
   }

   for (int n = 0; n < _nsites; n += 4) {
      int x, y, z, r;
      z = (n/4) / (_m*_m); r = (n/4) % (_m*_m);
      y = r / _m; x = r % _m;

      for (int i = 0; i < 4; i++) {
         lat->site[n+i].pos[0] = 2*x;
         lat->site[n+i].pos[1] = 2*y;
         lat->site[n+i].pos[2] = 2*z;
      }

      lat->site[n+1].pos[0] += 1; lat->site[n+1].pos[1] += 1;
      lat->site[n+2].pos[0] += 1; lat->site[n+2].pos[2] += 1;
      lat->site[n+3].pos[1] += 1; lat->site[n+3].pos[2] += 1;
   }

   for (int n = 0; n < _nsites; n++) {
      int x, y, z;
      x = lat->site[n].pos[0]; y = lat->site[n].pos[1]; z = lat->site[n].pos[2];

      #ifndef make_neighbors_of
      #define make_neighbors_of(id1, _x, _y, _z, _m) \
      do { \
         int id2 = find_id_by_pos((_x), (_y), (_z), (_m)); \
         if (id2 != -1) { \
            lat->site[id1].nn[lat->site[id1].nn_count++] = &lat->site[id2]; \
            lat->site[id2].nn[lat->site[id2].nn_count++] = &lat->site[id1]; \
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

void lattice_free(lattice *lat)
{
   free(lat);
}

void add_silver(site *s)
{
   s->state = _silver;
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

int main(int argc, char **argv)
{
   run_unit_tests();
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
         assert(find_id_by_pos(0,2,0,30) == 120);
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

   describe_lattice_new: {
      lattice *lat2 = lattice_new(2);
      lattice *lat30 = lattice_new(30);

      it_should_have_correct_m_and_nsites: {
         assert(lat2->m == 2);
         assert(lat2->nsites == 32);
         assert(lat30->m == 30);
         assert(lat30->nsites == 108000);
         puts(".");
      }

      it_should_create_nsites_site_objects: {
         for (int i = 0; i < lat2->nsites; i++) {
            assert(lat2->site[i].id == i);
         }
         for (int i = 0; i < lat30->nsites; i++) {
            assert(lat30->site[i].id == i);
         }
         puts(".");
      }

      it_should_determine_position_of_sites: {
         assert(find_id_by_pos(lat2->site[28].pos[0], lat2->site[28].pos[1],
            lat2->site[28].pos[2], 2) == 28);
         assert(find_id_by_pos(lat30->site[54346].pos[0], lat30->site[54346].pos[1],
            lat30->site[54346].pos[2], 30) == 54346);
         puts(".");
      }

      it_should_determine_neighbor_count_for_sites: {
         assert(lat2->site[0].nn_count == 3);
         assert(lat2->site[28].nn_count == 12);
         puts(".");
      }

      it_should_make_correct_neighbors_of_sites: {
         assert(lat2->site[0].nn[0] == &lat2->site[1]);
         assert(lat2->site[1].nn[0] == &lat2->site[0]);
         assert(lat2->site[0].nn[4] == NULL);
         puts(".");
      }

      lattice_free(lat2);
      lattice_free(lat30);
   }
}
