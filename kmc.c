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
void add_pvp(site *site);
void rm_pvp(site *site);

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
gsl_rng *_rng;
int _seed;
set *total_silver_set, *silver_set, *surface_set, *pvp_set;
double silver_energy[13] = { 0, 1, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4, 2.55, 2.7, 2.85, 3 };
double pvp_energy[13] = { 0, 0, 2.4, 2.5, 2.6, 2.5, 2.4, 2, 0.5, 0, 0, 0, 0 };

int main(int argc, char **argv)
{
   _rng = gsl_rng_alloc(gsl_rng_mt19937);
   gsl_rng_set(_rng, _seed);
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

void kmc(void)
{
   int index = 0;
   double p_sum[2*set_size(surface_set)+set_size(silver_set)+set_size(pvp_set)];
   site *obj_ary[2*set_size(surface_set)+set_size(silver_set)+set_size(pvp_set)];
   double ktot = 0.0;
   double rate_addition_silver = (_nsilver - set_size(total_silver_set)) / _V;
   double rate_addition_pvp = (_npvp - set_size(pvp_set)) / _V;
   set_enum *surface_set_enum = set_to_enum(surface_set);
   set_enum *silver_set_enum = set_to_enum(silver_set);
   set_enum *pvp_set_enum = set_to_enum(pvp_set);

   #ifndef fill_partial_sum
   #define fill_partial_sum(enumerator, rate) \
   do { \
      while (set_enum_next(enumerator)) { \
         int id = enumerator->value; \
         p_sum[index++] = ktot += rate; \
         obj_ary[index] = &lat.site[id]; \
      } \
   } while (0) 
   #endif

   fill_partial_sum(surface_set_enum, rate_addition_silver);
   set_enum_rewind(surface_set_enum);

   fill_partial_sum(surface_set_enum, rate_addition_pvp);
   set_enum_free(surface_set_enum);

   fill_partial_sum(silver_set_enum, lat.site[id].rate);
   set_enum_free(silver_set_enum);
   
   fill_partial_sum(pvp_set_enum, lat.site[id].rate);
   set_enum_free(pvp_set_enum);

   double r = ktot * gsl_rng_uniform(_rng);
   for (int i = 0; i < index; i++) {
      if (r <= p_sum[i]) {
         int a = 0;
         if (a <= i && i < (a += set_size(surface_set)))
            add_silver(obj_ary[i]);
         else if (a <= i && i < (a += set_size(surface_set)))
            add_pvp(obj_ary[i]);
         else if (a <= i && i < (a += set_size(silver_set)))
            rm_silver(obj_ary[i]);
         else if (a <= i && i < (a += set_size(pvp_set)))
            rm_pvp(obj_ary[i]);
         return;
      }
   }
   exit(0);
}

// void draw(void)
// {
//    if (!_draw) return;
// }

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
