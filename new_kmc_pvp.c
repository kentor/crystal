#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "new_lattice.h"
#include "flags.h"

void seed_crystal(int radius);
void grow(void);
void init(void);
void add_crystal_en(int site);
void rm_crystal_en(int site);
void add_pvp_en(int site);
void rm_pvp_en(int site);
void update_energy(int site);
void check(void);
void print_progress(int current, int final);
int is_on_face(int site);
int is_crystal_id(int x, int y, int z);

extern int m;
extern int Nsites;
extern int *ll;
extern struct site *lattice;
extern int center[3];
extern int crystals_hoc;
extern int sites_hoc;
extern int pvp_hoc;
extern int ncrystals;
extern int nsites;
extern int npvp;
extern gsl_rng *rng;

// default values
double T = 0.036;
double V = 1e34;
int radius = 4;
int tot_crystals = 2000;
int tot_pvp = 200;
int nsteps = 1000000;
int interval = 1000;
int max_crystals_drawn, max_pvps_drawn;
double e8 = 2.6;
double e9 = 2.7;
double b100 = 0.5;

double energy[13], energy_diff[13];
double penergy[13], penergy_diff[13];

int main(int argc, char **argv)
{
    int i;

    // initiate lattice / memory allocation
    init_lattice(argc, argv);

    // init random number generator
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(0));

    // override defaults
    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-')
        {
            char *flag = argv[i] + 1;
            char *arg = argv[i+1];

            if (flag_match(flag, 1, "T"))
                T = atof(arg);
            else if (flag_match(flag, 1, "V"))
                V = atof(arg);
            else if (flag_match(flag, 1, "n"))
                nsteps = atoi(arg);
            else if (flag_match(flag, 1, "r"))
                radius = atoi(arg);
            else if (flag_match(flag, 1, "t"))
                tot_crystals = atoi(arg);
            else if (flag_match(flag, 1, "tp"))
                tot_pvp = atoi(arg);
            else if (flag_match(flag, 1, "dm"))
                max_crystals_drawn = atoi(arg);
            else if (flag_match(flag, 1, "dpm"))
                max_pvps_drawn = atoi(arg);
            else if (flag_match(flag, 1, "i"))
                interval = atoi(arg);
            else if (flag_match(flag, 1, "e8"))
                e8 = atof(arg);
            else if (flag_match(flag, 1, "e9"))
                e9 = atof(arg);
            else if (flag_match(flag, 1, "b100"))
                b100 = atof(arg);
        }

    if (max_crystals_drawn == 0)
        max_crystals_drawn = tot_crystals;
    if (max_pvps_drawn == 0)
        max_pvps_drawn = tot_pvp;

    // init seed and energy scale
    init();
    
    printf("  Configuration:\n");
    printf("    Temperature (-T): %.4lf\n", T);
    printf("    Volume (-V): %.4le\n", V);
    printf("    Seed radius (-r): %d\n", radius);
    printf("    Max crystals (-t): %d\n", tot_crystals);
    printf("    Max pvp (-tp): %d\n", tot_pvp);
    printf("    Max crystals drawn (-dm): %d\n", max_crystals_drawn);
    printf("    Max PVPs drawn (-dpm): %d\n", max_pvps_drawn);
    printf("    Drawing interval (-i): %d\n", interval);
    printf("    Simulation steps (-n): %d\n", nsteps);
    printf("\n");

    printf("  Energy scale:\n");
    for (i = 0; i < 13; i++)
        printf("    energy[%d]: %.2lf\n", i, energy[i]);
    printf("\n\n");

    printf("  PVP energy scale:\n");
    for (i = 0; i < 13; i++)
        printf("    penergy[%d]: %.2lf\n", i, penergy[i]);
    printf("\n\n");

    printf("    PVP (100) bonus (-b100): %.2lf\n", b100);
    printf("\n\n");

    draw(max_crystals_drawn, max_pvps_drawn);

    int n;

    for (n = 1; n <= nsteps; n++)
    {
        grow();
        if (n % interval == 0)
        {
            print_progress(n, nsteps);
            draw(max_crystals_drawn, max_pvps_drawn);
        }
    }

    printf("\n\n");

    printf("Crystals adsorbed: %d, crystals floating: %d\n", ncrystals, tot_crystals - ncrystals);
    printf("PVPs adsorbed: %d, pvps floating: %d\n", npvp, tot_pvp - npvp);
}

void check(void)
{
    check_lattice();

    int i, j;

    for (i = 0; i < Nsites; i++)
        if (is_crystal(i))
        {
            double en = energy[lattice[i].neighbors];

            for (j = 0; j < 12; j++)
            {
                int neigh = lattice[i].nn[j];

                if (neigh == EMPTY || !is_crystal(neigh)) continue;

                en += energy_diff[lattice[neigh].neighbors];
            }

            assert(en == lattice[i].energy);
        }
}

void seed_crystal(int radius)
{
    int i, site;
    int dx, dy, dz;

    site = m/2;
    site = 4*(site + site*m + site*m*m);

    // save the coordinates of the center of the seed
    for (i = 0; i < 3; i++)
        center[i] = lattice[site].pos[i];

    for (i = 0; i < Nsites; i++)
    {
        dx = lattice[i].pos[0] - center[0];
        dy = lattice[i].pos[1] - center[1];
        dz = lattice[i].pos[2] - center[2];

        if (dx*dx + dy*dy + dz*dz < radius*radius)
            add_crystal_en(i);
    }
}

void print_progress(int current, int final)
{
    float f;
    int i, bars;
    int max_bars = 50;
    
    f = (float) current / final;

    bars = current / (final / max_bars);

    printf("  ");
    printf("%3d%% ", (int)(100*f));
    printf("[");
    for (i = 0; i < bars; i++)
        printf("=");
    for (i = bars; i < max_bars; i++)
        printf(" ");
    printf("] ");
    printf("%d/%d ", current, final);
    printf("\r");
    fflush(stdout);
}

void grow(void)
{
    int i, j;
    double rate_add, rate_add_pvp;
    double ktot = 0;
    double psum[2+ncrystals+npvp];
    int counter = 0;
    int crystals_ar[ncrystals];
    int pvp_ar[npvp];
    double r;

    // number of free crystals (total crystals - adsorbed crystals) divided by the volume
    rate_add = (tot_crystals - ncrystals)/V;
    rate_add_pvp = (tot_pvp - npvp)/V;

    // first entry in partial sum is the rate of addition multiplied by number of sites
    psum[counter++] = ktot = rate_add*nsites;

    // second entry is rate of addition of pvp multiplied by number of sites
    psum[counter++] = ktot += rate_add_pvp*nsites;

    // the rest of the partial sum contains rates for each individual crystal. loop through the linked list to get the crystals
    for (i = 0, j = crystals_hoc; i < ncrystals; i++, j = ll[j])
    {
        psum[counter++] = ktot += lattice[j].rate;
        crystals_ar[i] = j;
    }

    // loop through pvps and get their removal rates
    for (i = 0, j = pvp_hoc; i < npvp; i++, j = ll[j])
    {
        psum[counter++] = ktot += lattice[j].rate;
        pvp_ar[i] = j;
    }

    // for (i = 0; i < counter; i++)
    //     printf("psum[%d]\t%.4le\n", i, psum[i]);

    if (ktot == 0.0)
    {
        printf("!\n");
        exit(1);
    }

    // pick random double in the interval (0, ktot)
    r = ktot*gsl_rng_uniform(rng);

    for (i = 0; i < counter; i++) 
        if (r <= psum[i])
        {
            if (i == 0)
                add_crystal_en(rand_ss());
            else if (i == 1)
                add_pvp_en(rand_ss());
            else if (i < ncrystals + 2)
                rm_crystal_en(crystals_ar[i-2]);
            else
                rm_pvp_en(pvp_ar[i-ncrystals-2]);
            return;
        }

    printf("ERROR!!!!!!!!!\n");
    exit(2);
}

void init(void)
{
    int i;

    energy[0] = 0;
    energy[1] = 1.2;
    energy[2] = 2.0;
    energy[3] = 2.1;
    energy[4] = 2.2;
    energy[5] = 2.3;
    energy[6] = 2.4;
    energy[7] = 2.5;
    energy[8] = e8;
    energy[9] = e9;
    energy[10] = 2.8;
    energy[11] = 2.9;
    energy[12] = 3.0;

    energy_diff[0] = 0;
    for (i = 1; i < 13; i++)
        energy_diff[i] = energy[i] - energy[i-1];

    penergy[0] = 0;
    penergy[1] = 0;
    penergy[2] = 0;
    penergy[3] = 2.7;
    penergy[4] = 2.8;
    penergy[5] = 2.8;
    penergy[6] = 2.6;
    penergy[7] = 2.0;
    penergy[8] = 1.0;
    penergy[9] = 0;
    penergy[10] = 0;
    penergy[11] = 0;
    penergy[12] = 0;

    penergy_diff[0] = 0;
    for (i = 1; i < 13; i++)
        penergy_diff[i] = penergy[i] - penergy[i-1];

    seed_crystal(radius);
}

void add_crystal_en(int site)
{
    add_crystal(site);
    update_energy(site);
}

void rm_crystal_en(int site)
{
    rm_crystal(site);
    update_energy(site);
}

void add_pvp_en(int site)
{
    add_pvp(site);
    update_energy(site);
}

void rm_pvp_en(int site)
{
    rm_pvp(site);
    update_energy(site);
}

void update_energy(int site)
{
    // get the cell and neighboring cells of site, and store it into cell_list[]
    int cell, cell_list_size = 0, cell_list[27];
    int i, j, k, l;
    int x, y, z;

    cell = site/4;
    z = cell/(m*m); 
    y = (cell-z*m*m)/m;
    x = cell-z*m*m-y*m;

    for (k = z-1; k <= z+1; k++)
    {
        if (k < 0 || k >= m) continue;
        for (j = y-1; j <= y+1; j++)
        {
            if (j < 0 || j >= m) continue;
            for (i = x-1; i <= x+1; i++)
            {
                if (i < 0 || i >= m) continue;
                cell_list[cell_list_size++] = i + j*m + k*m*m;
            }
        }
    }

    // loop through every crystal of the cell list and update their energies and rates
    for (i = 0; i < cell_list_size; i++)
        for (j = cell_list[i]*4; j < cell_list[i]*4+4; j++)
        {
            lattice[j].energy = lattice[j].bonus_energy = 0;
            lattice[j].rate = 1;

            if (is_crystal(j))
            {
                // coordination energy of current crystal
                lattice[j].energy = energy[lattice[j].neighbors];

                // contributions to the energy of all neighboring crystals
                for (k = 0; k < 12; k++)
                {
                    int neigh = lattice[j].nn[k];

                    if (neigh == EMPTY) continue;

                    if (is_crystal(neigh))
                        lattice[j].energy += energy_diff[lattice[neigh].neighbors];
                    else if (is_pvp(neigh))
                        lattice[j].energy += penergy_diff[lattice[neigh].neighbors];
                }

                lattice[j].rate = exp(-(lattice[j].energy + lattice[j].bonus_energy)/T);
            }
            // else if (is_pvp(j) && lattice[j].neighbors <= 6 && is_on_face(j))
            // {
            //     lattice[j].energy = b100;
            //     lattice[j].rate = exp(-lattice[j].energy/T);
            // }
            else if (is_pvp(j))
            {
                // if (lattice[j].neighbors == 4 && is_on_face(j))
                //     lattice[j].energy = b100;
                // else
                    lattice[j].energy = penergy[lattice[j].neighbors];

                lattice[j].rate = exp(-lattice[j].energy/T);
            }
        }
}

int is_on_face(int site)
{
    int x = lattice[site].pos[0];
    int y = lattice[site].pos[1];
    int z = lattice[site].pos[2];

    // check bottom 4 sites
    if (x != 0 && x != 2*m-1 && y != 0 && y != 2*m-1 && z != 0)
        if (is_crystal_id(x-1,y,z-1) && is_crystal_id(x,y-1,z-1) && is_crystal_id(x+1,y,z-1) && is_crystal_id(x,y+1,z-1))
            return 1;

    // check left 4 sites
    if (x != 0 && y != 0 && y != 2*m-1 && z != 0 && z != 2*m-1)
        if (is_crystal_id(x-1,y-1,z) && is_crystal_id(x-1,y,z-1) && is_crystal_id(x-1,y,z+1) && is_crystal_id(x-1,y+1,z))
            return 1;

    // check top 4 sites
    if (x != 0 && x != 2*m-1 && y != 0 && y != 2*m-1 && z != 2*m-1)
        if (is_crystal_id(x-1,y,z+1) && is_crystal_id(x,y-1,z+1) && is_crystal_id(x,y+1,z+1) && is_crystal_id(x+1,y,z+1))
            return 1;

    // check right 4 sites
    if (x != 2*m-1 && y != 0 && y != 2*m-1 && z != 0 && z != 2*m-1)
        if (is_crystal_id(x+1,y-1,z) && is_crystal_id(x+1,y,z-1) && is_crystal_id(x+1,y,z+1) && is_crystal_id(x+1,y+1,z))
            return 1;

    // check front 4 sites
    if (x != 0 && x != 2*m-1 && y != 0 && z != 0 && z != 2*m-1)
        if (is_crystal_id(x-1,y-1,z) && is_crystal_id(x,y-1,z-1) && is_crystal_id(x,y-1,z+1) && is_crystal_id(x+1,y-1,z))
            return 1;

    // check back 4 sites
    if (x != 0 && x != 2*m-1 && y != 2*m-1 && z != 0 && z != 2*m-1)
        if (is_crystal_id(x-1,y+1,z) && is_crystal_id(x,y+1,z-1) && is_crystal_id(x+1,y+1,z) && is_crystal_id(x,y+1,z+1))
            return 1;

    return 0;
}

int is_crystal_id(int x, int y, int z)
{
    extern int ***id;

    return is_crystal(id[x][y][z]);
}
