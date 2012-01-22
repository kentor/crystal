#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
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

lattice lat;
double _T = 0.0375;
double _V = 1e34;
int _rad = 4;
int _nslvr = 2000;
int _npvp = 200;
int _nsteps = 1000000;
int _interval = 1000;
int _seed;
gsl_rng *_rng;
set *ttl_slvr_tbl, *slvr_tbl, *surf_tbl, *pvp_tbl;
bool _draw = true;
double slvr_nrg[13] = { 0, 1, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4, 2.55, 2.7, 2.85, 3 };
double pvp_nrg[13] = { 0, 0, 2.4, 2.5, 2.6, 2.5, 2.4, 2, 0.5, 0, 0, 0, 0 };
double slvr_nrg_d[13] = { 0 }, pvp_nrg_d[13] = { 0 };

int main(int argc, char **argv)
{

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
   int _size = 2*set_size(surf_tbl) + set_size(slvr_tbl) + set_size(pvp_tbl);
   site *obj_ary[_size];
   double p_sum[_size];
   double ktot = 0.0;
   double rate_addition_silver = (_nslvr - set_size(ttl_slvr_tbl)) / _V;
   double rate_addition_pvp = (_npvp - set_size(pvp_tbl)) / _V;

   #ifndef fill_p_sum
   #define fill_p_sum(table, rate) \
   do { \
      set_iter iter; \
      gpointer key, member; \
      set_iter_init(iter, table) \
      while (set_iter_next(iter, key, member)) { \
         p_sum[index++] = ktot += rate; \
         obj_ary[index] = _site(member); \
      } \
   } while (0) 
   #endif

   fill_p_sum(surf_tbl, rate_addition_silver);
   fill_p_sum(surf_tbl, rate_addition_pvp);
   fill_p_sum(slvr_tbl, _site(member)->rate);
   fill_p_sum(pvp_tbl, _site(member)->rate);

   if (_size != index) exit(0);

   double r = ktot * gsl_rng_uniform(_rng);
   for (int i = 0; i < index; i++) {
      if (r <= p_sum[i]) {
         int a = 0;
         if (a <= i && i < (a += set_size(surf_tbl)))
            add_silver(obj_ary[i]);
         else if (a <= i && i < (a += set_size(surf_tbl)))
            add_pvp(obj_ary[i]);
         else if (a <= i && i < (a += set_size(slvr_tbl)))
            rm_silver(obj_ary[i]);
         else if (a <= i && i < (a += set_size(pvp_tbl)))
            rm_pvp(obj_ary[i]);
         return;
      }
   }
   exit(0);
}

void initialize_kmc(int m)
{
   lat = new_lattice(m);
   ttl_slvr_tbl = set_new();
   slvr_tbl = set_new();
   surf_tbl = set_new();
   pvp_tbl = set_new();

   // _rng = gsl_rng_alloc(gsl_rng_mt19937);
   // gsl_rng_set(_rng, _seed);

   for (int i = 1; i < 13; i++) {
      slvr_nrg_d[i] = slvr_nrg[i] - slvr_nrg[i-1];
      pvp_nrg_d[i] = pvp_nrg[i] - pvp_nrg[i-1];
   }
}

void destroy_kmc(void)
{
   set_destroy(ttl_slvr_tbl);
   set_destroy(slvr_tbl);
   set_destroy(surf_tbl);
   set_destroy(pvp_tbl);
}

void add_silver(site *site)
{
   if (site->state != _surface) return;
   set_remove(surf_tbl, site);
   set_insert(ttl_slvr_tbl, site);
   if (site->neighbors < 11)
      set_insert(slvr_tbl, site);
   site->state = _silver;

   for (int i = 0; i < site->nn_count; i++) {
      site->nn[i]->neighbors++;
      if (site->nn[i]->state == _vacuum) {
         set_insert(surf_tbl, site->nn[i]);
         site->nn[i]->state = _surface;
      }
      else if (site->nn[i]->state == _silver && site->nn[i]->neighbors == site->nn[i]->nn_count) {
         set_remove(slvr_tbl, site->nn[i]);
      }
   }
}

void rm_silver(site *site)
{
   if (site->state != _silver) return;
   set_remove(slvr_tbl, site);
   set_remove(ttl_slvr_tbl, site);
   if (site->neighbors > 0) {
      set_insert(surf_tbl, site);
      site->state = _surface;
   }
   else {
      site->state = _vacuum;
   }

   for (int i = 0; i < site->nn_count; i++) {
      site->nn[i]->neighbors--;
      if (site->nn[i]->state == _surface && site->nn[i]->neighbors == 0) {
         set_remove(surf_tbl, site->nn[i]);
         site->nn[i]->state = _vacuum;
      }
      else if (site->nn[i]->state == _silver && !set_include(slvr_tbl, site->nn[i])) {
         set_insert(slvr_tbl, site->nn[i]);
      }
   }
}

void add_pvp(site *site)
{
   if (site->state != _surface) return;
   set_remove(surf_tbl, site);
   set_insert(pvp_tbl, site);
   site->state = _pvp;
}

void rm_pvp(site *site)
{
   if (site->state != _pvp) return;
   set_remove(pvp_tbl, site);
   if (site->neighbors > 0) {
      set_insert(surf_tbl, site);
      site->state = _surface;
   }
   else {
      site->state = _vacuum;
   }
}

void site_to_surface(site *site)
{
   if (site->state != _vacuum) return;
   set_insert(surf_tbl, site);
   site->state = _surface;
}

void update_nrg_around(site *site)
{
   set *nnn = set_new();

   for (int L1 = 0; L1 < site->nn_count; L1++) {
      for (int L2 = 0; L2 < site->nn[L1]->nn_count; L2++) {
         if (site->nn[L1]->nn[L2]->state == _silver) {
            set_insert(nnn, site->nn[L1]->nn[L2]);
         }
      }
   }
}
