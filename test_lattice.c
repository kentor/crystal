#include <stdio.h>
#include <assert.h>
#include "lattice.c"

int main()
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
         assert(find_id_by_pos(0,0,2,30) == 3600);
         assert(find_id_by_pos(58,58,58,30) == 107996);
         assert(find_id_by_pos(59,59,58,30) == 107997);
         assert(find_id_by_pos(59,58,59,30) == 107998);
         assert(find_id_by_pos(58,59,59,30) == 107999);
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
      lattice lat2 = new_lattice(2);
      lattice lat30 = new_lattice(30);

      it_should_have_correct_m_and_nsites: {
         assert(lat2.m == 2);
         assert(lat2.nsites == 32);
         assert(lat30.m == 30);
         assert(lat30.nsites == 108000);
         puts(".");
      }

      it_should_create_nsites_site_objects: {
         for (int i = 0; i < lat2.nsites; i++) {
            assert(lat2.site[i].id == i);
         }
         for (int i = 0; i < lat30.nsites; i++) {
            assert(lat30.site[i].id == i);
         }
         puts(".");
      }

      it_should_determine_position_of_sites: {
         assert(find_id_by_pos(lat2.site[28].pos[0], lat2.site[28].pos[1],
            lat2.site[28].pos[2], 2) == 28);
         assert(find_id_by_pos(lat30.site[54346].pos[0], lat30.site[54346].pos[1],
            lat30.site[54346].pos[2], 30) == 54346);
         puts(".");
      }

      it_should_determine_neighbor_count_for_sites: {
         assert(lat2.site[0].nn_count == 3);
         assert(lat2.site[28].nn_count == 12);
         puts(".");
      }

      it_should_make_correct_neighbors_of_sites: {
         assert(lat2.site[0].nn[0] == &lat2.site[1]);
         assert(lat2.site[1].nn[0] == &lat2.site[0]);
         assert(lat2.site[0].nn[4] == NULL);
         puts(".");
      }
   }
   return 0;
}
