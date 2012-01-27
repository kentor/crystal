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
      site_t *site2 = new_lattice(2);
      site_t *site30 = new_lattice(30);

      it_should_create_nsites_site_objects: {
         for (int i = 0; i < 4*2*2*2; i++) {
            assert(site2[i].id == i);
         }
         for (int i = 0; i < 4*30*30*30; i++) {
            assert(site30[i].id == i);
         }
         puts(".");
      }

      it_should_determine_position_of_sites: {
         assert(find_id_by_pos(site2[28].pos[0], site2[28].pos[1],
            site2[28].pos[2], 2) == 28);
         assert(find_id_by_pos(site30[54346].pos[0], site30[54346].pos[1],
            site30[54346].pos[2], 30) == 54346);
         puts(".");
      }

      it_should_determine_neighbor_count_for_sites: {
         assert(site2[0].nn_count == 3);
         assert(site2[28].nn_count == 12);
         puts(".");
      }

      it_should_make_correct_neighbors_of_sites: {
         assert(site2[0].nn[0] == 1);
         assert(site2[1].nn[0] == 0);
         assert(site2[0].nn[4] == -1);
         puts(".");
      }
   }
   return 0;
}
