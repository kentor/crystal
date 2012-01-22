#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include "lattice.h"
#include "set.h"

void initialize_kmc(int m);
void destroy_kmc(void);
void add_silver(site *site);
void rm_silver(site *site);
void site_to_surface(site *site);
void run_unit_tests(void);

lattice lat;
double _T = 0.0375;
double _V = 1e34;
int _rad = 4;
int _nsilver = 2000;
int _npvp = 200;
int _nsteps = 1000000;
int _interval = 1000;
bool _draw = true;
set *total_silver_set;
set *silver_set;
set *surface_set;
set *pvp_set;
double silver_energy[13] = { 0, 1, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4, 2.55, 2.7, 2.85, 3 };
double pvp_energy[13] = { 0, 0, 2.4, 2.5, 2.6, 2.5, 2.4, 2, 0.5, 0, 0, 0, 0 };

int main(int argc, char **argv)
{
   run_unit_tests();
   // initialize_kmc(30);
   /*
   draw();
   for (int step = 1; step <= _nsteps; step++) {
      kmc();
      if (step % _interval == 0) {
         print_progress();
         draw();
      }
   }
   */
   return 0;
}
/*
void kmc(void)
{
   int index = 0;
   double partial_sum[2*set_size(surface_set)+set_size(silver_set)+set_size(pvp_set)][2];
   double ktot = 0.0;
   set_enum *surface_set_enum = set_to_enum(surface_set);
   set_enum *silver_set_enum = set_to_enum(silver_set);
   set_enum *pvp_set_enum = set_to_enum(pvp_set);

   double rate_addition_silver = (_nsilver - set_size(total_silver_set)) / _V;
   double rate_addition_pvp = (_npvp - set_size(pvp_set)) / _V;

   while (set_enum_next(surface_set_enum)) {
      int id = surface_set_enum->value;
      partial_sum[index++][0] = ktot += rate_addition_silver;
      partial_sum[index][1] = &lat.site[id];
   }
   set_enum_rewind(surface_set_enum);
   while (set_enum_next(surface_set_enum)) {
      int id = surface_set_enum->value;
      partial_sum[index++][0] = ktot += rate_addition_pvp;
      partial_sum[index][1] = &lat.site[id];
   }
   set_enum_free(surface_set_enum);

   while (set_enum_next(silver_set_enum)) {
      int id = silver_set_enum->value;
      partial_sum[index++][0] = ktot += lat.site[id].rate;
      partial_sum[index][1] = &lat.site[id];
   }
   set_enum_free(silver_set_enum);
   
   while (set_enum_next(pvp_set_enum)) {
      int id = pvp_set_enum->value;
      psum[psumsize++] = ktot += lat.site[id].rate;
      object_ary[psumsize] = &lat.site[id];
   }
   set_enum_free(pvp_set_enum);

   r = ktot * gsl_rng_uniform(_rng);

   for (int i = 0; i < psumsize; i++) {
      if (r <= psum[i]) {
         int a = 0;
         if (a <= i && i < (a += set_size(surface_set)))
            add_silver(psum[i][1]);
         else if (a <= i && i < (a += set_size(surface_set)))
            add_pvp(psum[i][1]);
         else if (a <= i && i < (a += set_size(silver_set)))
            rm_silver(psum[i][1]);
         else if (a <= i && i < (a += set_size(pvp_set)))
            rm_pvp(psum[i][1]);
         return;
      }
   }
   exit(0);
}

void draw(void)
{
   if (!_bdraw) return;
}
*/
void initialize_kmc(int m)
{
   lat = new_lattice(m);
   total_silver_set = new_set(lat.nsites);
   silver_set = new_set(lat.nsites);
   surface_set = new_set(lat.nsites);
   pvp_set = new_set(lat.nsites);
}

void destroy_kmc(void)
{
   free(total_silver_set);
   free(silver_set);
   free(surface_set);
   free(pvp_set);
}

void add_silver(site *site)
{
   if (site->state != _surface) return;
   set_delete(surface_set, site->id);
   set_insert(total_silver_set, site->id);
   if (site->neighbors < 11)
      set_insert(silver_set, site->id);
   site->state = _silver;

   for (int i = 0; i < site->nn_count; i++) {
      site->nn[i]->neighbors++;
      if (site->nn[i]->state == _vacuum) {
         set_insert(surface_set, site->nn[i]->id);
         site->nn[i]->state = _surface;
      }
      else if (site->nn[i]->state == _silver && site->nn[i]->neighbors == site->nn[i]->nn_count) {
         set_delete(silver_set, site->nn[i]->id);
      }
   }
}

void rm_silver(site *site)
{
   if (site->state != _silver) return;
   set_delete(silver_set, site->id);
   set_delete(total_silver_set, site->id);
   if (site->neighbors > 0) {
      set_insert(surface_set, site->id);
      site->state = _surface;
   }
   else {
      site->state = _vacuum;
   }

   for (int i = 0; i < site->nn_count; i++) {
      site->nn[i]->neighbors--;
      if (site->nn[i]->state == _surface && site->nn[i]->neighbors == 0) {
         set_delete(surface_set, site->nn[i]->id);
         site->nn[i]->state = _vacuum;
      }
      else if (site->nn[i]->state == _silver && !set_include(silver_set, site->nn[i]->id)) {
         set_insert(silver_set, site->nn[i]->id);
      }
   }
}

void add_pvp(site *site)
{
   if (site->state != _surface) return;
   set_delete(surface_set, site->id);
   set_insert(pvp_set, site->id);
   site->state = _pvp;
}

void rm_pvp(site *site)
{
   if (site->state != _pvp) return;
   set_delete(pvp_set, site->id);
   if (site->neighbors > 0) {
      set_insert(surface_set, site->id);
      site->state = _surface;
   }
   else {
      site->state = _vacuum;
   }
}

void site_to_surface(site *site)
{
   if (site->state != _vacuum) return;
   set_insert(surface_set, site->id);
   site->state = _surface;
}

#include <assert.h>
void run_unit_tests(void)
{
   describe_site_to_surface: {
      initialize_kmc(2);

      it_should_make_site_to_surface: {
         assert(lat.site[0].state == _vacuum);
         site_to_surface(&lat.site[0]);
         site_to_surface(&lat.site[29]);
         assert(lat.site[0].state == _surface);
         assert(lat.site[29].state == _surface);
         puts(".");
      }

      it_should_add_site_to_surface_set: {
         assert(set_include(surface_set, 0));
         assert(set_include(surface_set, 29));
         puts(".");
      }

      puts("-");
      destroy_kmc();
   }

   describe_add_silver: {
      initialize_kmc(2);

      it_should_not_do_anything_if_site_is_not_surface: {
         assert(lat.site[28].state == _vacuum);
         add_silver(&lat.site[28]);
         assert(lat.site[28].state == _vacuum);
         puts(".");
      }

      it_should_make_site_into_silver: {
         site_to_surface(&lat.site[28]);
         assert(!set_include(silver_set, 28));
         assert(!set_include(total_silver_set, 28));
         add_silver(&lat.site[28]);
         assert(lat.site[28].state == _silver);
         puts(".");
      }

      it_should_make_neighbors_into_surface_sites: {
         for (int i = 0; i < lat.site[28].nn_count; i++) {
            assert(lat.site[28].nn[i]->state == _surface);
         }
         puts(".");
      }

      it_should_increase_silver_neighbor_count_of_neighbors: {
         for (int i = 0; i < lat.site[28].nn_count; i++) {
            assert(lat.site[28].nn[i]->neighbors == 1);
         }
         puts(".");
      }

      it_should_add_site_to_silver_and_total_silver_set: {
         assert(set_include(silver_set, 28));
         assert(set_include(total_silver_set, 28));
         puts(".");
      }

      it_should_remove_site_from_surface_set: {
         assert(!set_include(surface_set, 28));
         puts(".");
      }

      it_should_add_neighbors_to_surface_set: {
         for (int i = 0; i < lat.site[28].nn_count; i++) {
            assert(set_include(surface_set, lat.site[28].nn[i]->id));
         }
         puts(".");
      }

      it_should_remove_site_from_silver_set_when_it_becomes_bulk: {
         for (int i = 0; i < lat.site[28].nn_count; i++) {
            add_silver(lat.site[28].nn[i]);
         }
         assert(!set_include(silver_set, 28));
         assert(set_include(total_silver_set, 28));
         puts(".");
      }

      puts("-");
      destroy_kmc();
   }

   describe_rm_silver: {
      initialize_kmc(2);

      it_should_not_do_anything_if_site_is_not_silver: {
         assert(lat.site[28].state == _vacuum);
         rm_silver(&lat.site[28]);
         assert(lat.site[28].state == _vacuum);
         puts(".");
      }

      it_should_make_site_into_surface_or_vacuum_depending_on_neighbors: {
         site_to_surface(&lat.site[0]);
         site_to_surface(&lat.site[1]);
         add_silver(&lat.site[0]);
         add_silver(&lat.site[1]);
         assert(lat.site[0].state == _silver);
         assert(lat.site[1].state == _silver);
         assert(lat.site[4].state == _surface);
         rm_silver(&lat.site[1]);
         assert(lat.site[0].state == _silver);
         assert(lat.site[1].state == _surface);
         assert(lat.site[4].state == _vacuum);
         puts(".");
      }

      it_should_remove_silver_from_silver_set: {
         assert(!set_include(silver_set, 1));
         assert(!set_include(total_silver_set, 1));
         puts(".");
      }

      it_should_add_site_to_surface_set_if_its_not_vacuum: {
         assert(set_include(surface_set, 1));
         assert(!set_include(surface_set, 4));
         assert(!set_include(surface_set, 0));
         puts(".");
      }

      it_should_decrease_neighbor_count_of_neighbors: {
         assert(lat.site[4].neighbors == 0);
         assert(lat.site[0].neighbors == 0);
         puts(".");
      }

      it_should_insert_previously_bulk_crystals_to_silver_set: {
         destroy_kmc();
         initialize_kmc(2);
         site_to_surface(&lat.site[28]);
         add_silver(&lat.site[28]);
         for (int i = 0; i < lat.site[28].nn_count; i++) {
            add_silver(lat.site[28].nn[i]);
         }
         assert(lat.site[28].neighbors == 12);
         assert(!set_include(silver_set, 28));
         rm_silver(&lat.site[29]);
         assert(set_include(silver_set, 28));
         puts(".");
      }

      puts("-");
      destroy_kmc();
   }

   describe_add_pvp: {
      initialize_kmc(2);

      it_should_make_site_into_pvp_and_insert_into_pvp_set: {
         site_to_surface(&lat.site[0]);
         add_pvp(&lat.site[0]);
         assert(lat.site[0].state == _pvp);
         assert(set_include(pvp_set, 0));
         assert(!set_include(surface_set, 0));
         puts(".");
      }

      puts("-");
      destroy_kmc();
   }

   describe_rm_pvp: {
      initialize_kmc(2);

      it_should_remove_pvp_and_and_remove_from_sets_appropriately: {
         site_to_surface(&lat.site[0]);
         site_to_surface(&lat.site[28]);
         add_silver(&lat.site[28]);
         add_pvp(&lat.site[0]);
         add_pvp(&lat.site[29]);
         rm_pvp(&lat.site[0]);
         rm_pvp(&lat.site[29]);
         assert(lat.site[0].state == _vacuum);
         assert(lat.site[29].state == _surface);
         assert(!set_include(surface_set, 0));
         assert(set_include(surface_set, 29));
         puts(".");
      }

      destroy_kmc();
   }
}
