#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "lattice.h"
#include "parse_arg.h"
#include "ll.h"

void parse_and_print_vars(int argc, char **argv);
void initialize_kmc(int m, int rad);
void kmc(void);

void add_silver(int site_id);
void rm_silver(int site_id);
void add_pvp(int site_id);
void rm_pvp(int site_id);
void update_nrg_around(int site_id);
void draw(int max_slvr_d, int max_pvp_d, char *fn);

site_t *site;
double _T = 0.0350;
double _V = 1e34;
int _m = 30;
int _rad = 4;
int _max_slvr = 10000;
int _max_pvp = 10000;
int _nsteps = 10000000;
int _interval = 10000;
int _seed;
int _max_slvr_d, _max_pvp_d;
int _center[3];
gsl_rng *_rng;
bool _draw = true;
double slvr_nrg[13] = { 0, 1, 1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4, 2.55, 2.7, 2.85, 3 };
double pvp_nrg[13] = { 0, 0, 2.4, 2.5, 2.6, 2.5, 2.4, 2, 0.5, 0, 0, 0, 0 };
double slvr_nrg_d[13], pvp_nrg_d[13];
char *fn = "movie.xyz";

int *ll;
list_t _slvr_l = { .head = null, .size = 0 },
   _surf_l = { .head = null, .size = 0 },
   _pvp_l = { .head = null, .size = 0 };
int _ttlslvr = 0;

#define not(_state, _id) site[_id].state != _state
#define is(_state, _id) site[_id].state == _state
#define to(_state, _id) site[_id].state = _state

int main(int argc, char **argv)
{
   parse_and_print_vars(argc, argv);
   initialize_kmc(_m, _rad);

   draw(3000, 3000, fn);
   for (int step = 1; step <= _nsteps; step++) {
      kmc();
      if (step % _interval == 0) {
         draw(3000, 3000, fn);
      }
   }
   return 0;
}

void kmc(void)
{
   int index = 0;
   int size = 2*_surf_l.size + _slvr_l.size + _pvp_l.size;
   int ids[size];
   double p_sum[size];
   double ktot = 0.0;
   double rate_addition_silver = (_max_slvr - _ttlslvr) / _V;
   double rate_addition_pvp = (_max_pvp - _pvp_l.size) / _V;

   for (int i = _surf_l.head; i != null; i = ll[i]) {
      ids[index] = i;
      p_sum[index++] = ktot += rate_addition_silver;
   }
   for (int i = _surf_l.head; i != null; i = ll[i]) {
      ids[index] = i;
      p_sum[index++] = ktot += rate_addition_pvp;
   }
   for (int i = _slvr_l.head; i != null; i = ll[i]) {
      ids[index] = i;
      p_sum[index++] = ktot += site[i].rate;
   }
   for (int i = _pvp_l.head; i != null; i = ll[i]) {
      ids[index] = i;
      p_sum[index++] = ktot += site[i].rate;
   }

   double r = ktot * gsl_rng_uniform(_rng);
   for (int i = 0; i < index; i++) {
      if (r <= p_sum[i]) {
         int a = 0;
         if (a <= i && i < (a += _surf_l.size))
            add_silver(ids[i]);
         else if (a <= i && i < (a += _surf_l.size))
            add_pvp(ids[i]);
         else if (a <= i && i < (a += _slvr_l.size))
            rm_silver(ids[i]);
         else if (a <= i && i < (a += _pvp_l.size))
            rm_pvp(ids[i]);
         return;
      }
   }

   abort();
}

void initialize_kmc(int m, int rad)
{
   int nsites = 4*m*m*m;
   site = new_lattice(m);
   new_ll(ll, nsites);

   _rng = gsl_rng_alloc(gsl_rng_mt19937);
   gsl_rng_set(_rng, _seed);

   for (int i = 1; i < 13; i++) {
      slvr_nrg_d[i] = slvr_nrg[i] - slvr_nrg[i-1];
      pvp_nrg_d[i] = pvp_nrg[i] - pvp_nrg[i-1];
   }

   int id = m/2;
   id = 4*(id + id*m + id*m*m);
   for (int i = 0; i < 3; i++) {
      _center[i] = site[id].pos[i];
   }
   for (int i = 0; i < nsites; i++) {
      int dx, dy, dz;
      dx = site[i].pos[0] - _center[0];
      dy = site[i].pos[1] - _center[1];
      dz = site[i].pos[2] - _center[2];
      if (dx*dx + dy*dy + dz*dz < rad*rad) {
         add_silver(i);
      }
   }
}

void add_silver(int id)
{
   if (is(surface, id)) {
      ll_remove(ll, id, _surf_l);
   }
   else if (not(vacuum, id)) {
      return;
   }
   if (site[id].neighbors < site[id].nn_count) {
      ll_insert(ll, id, _slvr_l);
   }
   to(silver, id);

   for (int i = 0; i < site[id].nn_count; i++) {
      int neigh = site[id].nn[i];
      site[neigh].neighbors++;
      if (is(vacuum, neigh)) {
         ll_insert(ll, neigh, _surf_l);
         to(surface, neigh);
      }
      else if (is(silver, neigh) && site[neigh].neighbors == site[neigh].nn_count) {
         ll_remove(ll, neigh, _slvr_l);
      }
   }
   _ttlslvr++;
   update_nrg_around(id);
}

void rm_silver(int id)
{
   if (not(silver, id)) return;
   if (site[id].neighbors < site[id].nn_count) {
      ll_remove(ll, id, _slvr_l);
   }
   if (site[id].neighbors > 0) {
      ll_insert(ll, id, _surf_l);
      to(surface, id);
   }
   else {
      to(vacuum, id);
   }

   for (int i = 0; i < site[id].nn_count; i++) {
      int neigh = site[id].nn[i];
      site[neigh].neighbors--;
      if (is(surface, neigh) && site[neigh].neighbors == 0) {
         ll_remove(ll, neigh, _surf_l);
         to(vacuum, neigh);
      }
      else if (is(silver, neigh) && site[neigh].neighbors+1 == site[neigh].nn_count) {
         ll_insert(ll, neigh, _slvr_l);
      }
   }
   _ttlslvr--;
   update_nrg_around(id);
}

void add_pvp(int id)
{
   if (is(surface, id)) {
      ll_remove(ll, id, _surf_l);
   }
   else if (not(vacuum, id)) {
      return;
   }
   ll_insert(ll, id, _pvp_l);
   to(pvp, id);
   update_nrg_around(id);
}

void rm_pvp(int id)
{
   if (not(pvp, id)) return;
   ll_remove(ll, id, _pvp_l);
   if (site[id].neighbors > 0) {
      ll_insert(ll, id, _surf_l);
      to(surface, id);
   }
   else {
      to(vacuum, id);
   }
   update_nrg_around(id);
}

void update_nrg_around(int id)
{
   // #ifndef update_nrg_silver
   // #define update_nrg_silver(_id) \
   // do { \
   //    site[_id].energy = slvr_nrg[site[_id].neighbors]; \
   //    for (int i = 0; i < site[_id].nn_count; i++) { \
   //       int _neigh = site[_id].nn[i]; \
   //       if (is(silver, _neigh)) { \
   //          site[_id].energy += slvr_nrg_d[site[_neigh].neighbors]; \
   //       } \
   //       else if (is(pvp, _neigh)) { \
   //          site[_id].energy += pvp_nrg_d[site[_neigh].neighbors]; \
   //       } \
   //    } \
   //    site[_id].rate = exp(-site[_id].energy / _T); \
   // } while(0)
   // #endif

   // #ifndef update_nrg_pvp
   // #define update_nrg_pvp(_id) \
   // do { \
   //    site[_id].energy = pvp_nrg[site[_id].neighbors]; \
   //    site[_id].rate = exp(-site[_id].energy / _T); \
   // } while(0)
   // #endif

   // if (is(silver, id)) update_nrg_silver(id);
   // else if (is(pvp, id)) update_nrg_pvp(id);

   // for (int neigh = 0; neigh < site[id].nn_count; neigh++) {
   //    if (is(silver, neigh)) update_nrg_silver(neigh);
   //    else if (is(pvp, neigh)) update_nrg_pvp(neigh);
   // }

   // for (int nneigh = 0; nneigh < site[id].nnn_count; nneigh++) {
   //    if (is(silver, nneigh)) update_nrg_silver(nneigh);
   // }
   int cell, size = 0, cell_list[27];
   cell = id/4;
   int x, y, z, r, m = 30;
   z = cell / (m*m); r = cell % (m*m);
   y = r / m; x = r % m;

   for (int k = z-1; k <= z+1; k++) {
      if (k < 0 || k >= m) continue;
      for (int j = y-1; j <= y+1; j++) {
         if (j < 0 || j >= m) continue;
         for (int i = x-1; i <= x+1; i++) {
            if (i < 0 || i >= m) continue;
            cell_list[size++] = i + j*m + k*m*m;
         }
      }
   }

   for (int i = 0; i < size; i++) {
      for (int j = cell_list[i]*4; j < cell_list[i]*4+4; j++) {
         if (is(silver, j)) {
            site[j].energy = slvr_nrg[site[j].neighbors];

            for (int k = 0; k < site[j].nn_count; k++) {
               int neigh = site[j].nn[k];
               if (is(silver, neigh)) {
                  site[j].energy += slvr_nrg_d[site[neigh].neighbors];
               }
               else if (is(pvp, neigh)) {
                  site[j].energy += pvp_nrg_d[site[neigh].neighbors];
               }
            }
            
            site[j].rate = exp(-site[j].energy / _T);
         }
         else if (is(pvp, j)) {
            site[j].energy = pvp_nrg[site[j].neighbors];
            site[j].rate = exp(-site[j].energy / _T);
         }
      }
   }
}

void draw(int max_slvr_d, int max_pvp_d, char *fn)
{
   if (!_draw) return;

   static FILE *fp = NULL;
   int count = 0;

   if (fp == NULL) fp = fopen(fn, "w");

   fprintf(fp, "%d\n\n", max_slvr_d + max_pvp_d);

   for (int i = _slvr_l.head; i != null; i = ll[i]) {
      fprintf(fp, "C");
      for (int j = 0; j < 3; fprintf(fp, " %d", site[i].pos[j++]));
      fprintf(fp, "\n");
      count++;

      if (count == max_slvr_d) break;
   }

   for (int i = 0; i < max_slvr_d - count; i++) {
      fprintf(fp, "C");
      for (int j = 0; j < 3; fprintf(fp, " %d", _center[j++]));
      fprintf(fp, "\n");
   }

   count = 0;
   for (int i = _pvp_l.head; i != null; i = ll[i]) {
      fprintf(fp, "F");
      for (int j = 0; j < 3; fprintf(fp, " %d", site[i].pos[j++]));
      fprintf(fp, "\n");
      count++;

      if (count == max_pvp_d) break;
   }

   for (int i = 0; i < max_pvp_d - count; i++) {
      fprintf(fp, "F");
      for (int j = 0; j < 3; fprintf(fp, " %d", _center[j++]));
      fprintf(fp, "\n");
   }
}

void parse_and_print_vars(int argc, char **argv)
{
   for (int i = 0; i < argc; i++) {
      if (argv[i][0] == '-') {
         char *flag = argv[i]+1;
         char *arg  = argv[i+1];
         parse_arg(flag, "m", _m, atoi(arg));
         parse_arg(flag, "T", _T, atof(arg));
         parse_arg(flag, "V", _V, atof(arg));
         parse_arg(flag, "n", _nsteps, atoi(arg));
         parse_arg(flag, "i", _interval, atoi(arg));
      }
   }

   printf("\n");
   printf("  Unit cells per edge (-m): %d\n", _m);
   printf("  Temperature (-T): %0.4lf\n", _T);
   printf("  Volume (-V): %0.4le\n", _V);
   printf("  Seed radius (-r): %d\n", _rad);
   printf("  Simulation steps (-n): %d\n", _nsteps);
   printf("  Drawing interval (-i): %d\n", _interval);
   printf("\n");
}