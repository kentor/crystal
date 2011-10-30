#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "new_lattice.h"
#include "flags.h"

int m;                      // number of cells per edge length
int Nsites;                 // number of total sites
const int A = 2;            // spacing between cells
struct site *lattice;       // lattice, an array of sites, of course
int center[3];              // coord of crystal seed center

int *ll;
int ***id;
int crystals_hoc = EMPTY;
int sites_hoc = EMPTY;
int pvp_hoc = EMPTY;
int ncrystals = 0;
int nsites = 0;
int npvp = 0;

int bdraw;                  // booleans

FILE *fcoords;
char fncoords[50];

gsl_rng *rng;               // random number generator

void init_lattice(int argc, char **argv)
{
    int i, j, k, n;
    int seed;

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    
    // default values for lattice
    m = 30;
    bdraw = 1;
    sprintf(fncoords, "coords%d.xyz", time(NULL));
    seed = time(0);

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-')
        {
            char *flag = argv[i] + 1;
            char *arg = argv[i+1];

            if (flag_match(flag, 1, "m"))
                m = atoi(arg);
            else if (flag_match(flag, 1, "d"))
                bdraw = !bdraw;
            else if (flag_match(flag, 1, "o"))
                sprintf(fncoords, "%s", arg);
            else if (flag_match(flag, 1, "se"))
                seed = atoi(arg);
        }

    Nsites = 4*m*m*m;
    gsl_rng_set(rng, seed);

    printf("\n");
    printf(" FCC lattice information:\n");
    printf("   # of unit cells per edge length (-m): %d\n", m);
    printf("   # of total sites: %d\n", Nsites);
    printf("   Drawing (-d): %s\n", bdraw ? "true" : "false");
    if (bdraw)
    {
        printf("     Output file (-o): %s\n", fncoords);
    }
    printf("   Seed (-se): %d\n", seed);
    printf("\n");

    if (bdraw) fcoords = fopen(fncoords, "w");

    // initialize lattice, an array of structs
    lattice = malloc(Nsites*sizeof(struct site));

    // initialize state for each site
    for (i = 0; i < Nsites; i++)
    {
        lattice[i].state = VACUUM;
        lattice[i].neighbors = 0;
        for (j = 0; j < 12; j++)
            lattice[i].nn[j] = EMPTY;
        lattice[i].edge = 0;
        lattice[i].energy = 0;
        lattice[i].bonus_energy = 0;
        lattice[i].rate = 1;
    }

    // calculate coordinate of each lattice site arranged in FCC. do this cell by cell
    for (n = 0; n < Nsites/4; n++)
    {
        int site;
        int x, y, z;
        const int disp = A/2;

        // get coordinates of the lowest corner of cell
        z = n/(m*m);
        y = (n-z*m*m)/m;
        x = n-z*m*m-y*m;
        
        // smallest site id in the cell
        site = 4*n;

        // initialize all 4 particles in the corner of the cell
        for (i = 0; i < 4; i++)
        {
            lattice[site+i].pos[0] = x*A;
            lattice[site+i].pos[1] = y*A;
            lattice[site+i].pos[2] = z*A;
        }

        // displace 3 of them to get the FCC unit cell
        lattice[site+1].pos[0] += disp;
        lattice[site+1].pos[1] += disp;
        lattice[site+2].pos[0] += disp;
        lattice[site+2].pos[2] += disp;
        lattice[site+3].pos[1] += disp;
        lattice[site+3].pos[2] += disp;
    }     

    // determine nearest neighbors of each site
    for (n = 0; n < Nsites; n++)
    {
        int count = 0;
        int cell, site;
        int x, y, z;
        int dx, dy, dz;

        // determine cell of site
        cell = n/4;

        // decompose cell into x y and z
        z = cell/(m*m); 
        y = (cell-z*m*m)/m;
        x = cell-z*m*m-y*m;

        // loop through neighboring cells
        for (k = z-1; k <= z+1; k++)
        {
            if (k < 0 || k >= m)
            {
                lattice[n].edge = 1;
                continue;
            }

            for (j = y-1; j <= y+1; j++)
            {
                if (j < 0 || j >= m)
                {
                    lattice[n].edge = 1;
                    continue;
                }

                for (i = x-1; i <= x+1; i++)
                {
                    if (i < 0 || i >= m)
                    {
                        lattice[n].edge = 1;
                        continue;
                    }

                    // determine current neighboring cell number
                    cell = i + j*m + k*m*m;

                    // loop through the 4 sites in the cell and see if it's a neighbor site
                    for (site = cell*4; site < cell*4+4; site++)
                    {
                        if (site == n) continue;

                        dx = lattice[n].pos[0] - lattice[site].pos[0];
                        dy = lattice[n].pos[1] - lattice[site].pos[1];
                        dz = lattice[n].pos[2] - lattice[site].pos[2];

                        // if r^2 is equal to 0.5*A^2 then it is a neighbor
                        if ((dx*dx + dy*dy + dz*dz) == A*A/2)
                            lattice[n].nn[count++] = site;
                    }
                }
            }
        }
    }

    // allocate memory for linked list
    ll = malloc(Nsites*sizeof(int));
    for (i = 0; i < Nsites; i++)
        ll[i] = EMPTY;

    // initialize id lookup table by coordinates
    id = malloc(2*m*sizeof(int **));
    for (i = 0; i < 2*m; i++)
        id[i] = malloc(2*m*sizeof(int *));
    for (i = 0; i < 2*m; i++)
        for (j = 0; j < 2*m; j++)
            id[i][j] = malloc(2*m*sizeof(int));
        
    for (i = 0; i < Nsites; i++)
    {
        int x = lattice[i].pos[0];
        int y = lattice[i].pos[1];
        int z = lattice[i].pos[2];

        id[x][y][z] = i;
    }
}

// void seed_crystal(int radius)
// {
//     int i, site;
//     int dx, dy, dz;
// 
//     site = m/2;
//     site = 4*(site + site*m + site*m*m);
// 
//     // save the coordinates of the center of the seed
//     for (i = 0; i < 3; i++)
//         center[i] = lattice[site].pos[i];
// 
//     for (i = 0; i < Nsites; i++)
//     {
//         dx = lattice[i].pos[0] - center[0];
//         dy = lattice[i].pos[1] - center[1];
//         dz = lattice[i].pos[2] - center[2];
// 
//         if (dx*dx + dy*dy + dz*dz < radius*radius)
//             add_crystal(i);
//     }
// }

void seed_crystal_cube(int radius)
{
    int i, j, k;
    int cent, site, cell;

    // center cell index, rounded down
    cent = m/2;
    
    // allocate memory for linked list
    ll = malloc(Nsites*sizeof(int));
    for (i = 0; i < Nsites; i++)
        ll[i] = EMPTY;
    site = 4*(cent + cent*m + cent*m*m);

    // save coordinate of the center
    for (i = 0; i < 3; i++)
        center[i] = lattice[site].pos[i];

    // make all lattice sites into sites within the cube of specific radius
    for (k = cent - radius; k < cent + radius; k++)
        for (j = cent - radius; j < cent + radius; j++)
            for (i = cent - radius; i < cent + radius; i++)
            {
                cell = i + j*m + k*m*m;

                for (site = cell*4; site < cell*4+4; site++)
                    add_crystal(site);
            }
}

void add_crystal(int site)
{
    int i;

    assert(is_vacuum(site) || is_surface_site(site));

    if (is_surface_site(site))
        rm_from_ss_ll(site);

    add_to_sc_ll(site);

    lattice[site].state = CRYSTAL;

    // increment the crystal neighbor count on all the neighboring sites
    for (i = 0; i < 12; i++)
    {
        int neigh = lattice[site].nn[i];
        
        // quit program if nearest neighbor of current site is outside the sim. box
        if (neigh == EMPTY)
        {
            printf("Reached end of box. Exiting!\n");
            exit(1);
        }
        
        lattice[neigh].neighbors++;
       
        // if neighbor is vacuum, turn it into a surface site
        if (is_vacuum(neigh))
        {
            add_to_ss_ll(neigh);
            lattice[neigh].state = SURFACE_SITE;
        }
    }
}

void rm_crystal(int site)
{
    int i;
    
    assert(is_crystal(site));

    rm_from_sc_ll(site);

    if (lattice[site].neighbors == 0)
        lattice[site].state = VACUUM;
    else
    {
        add_to_ss_ll(site);
        lattice[site].state = SURFACE_SITE;
    }

    // decrement the # of neighbors for each neighboring site
    for (i = 0; i < 12; i++)
    {
        int neigh = lattice[site].nn[i];

        lattice[neigh].neighbors--;

        if (is_surface_site(neigh) && lattice[neigh].neighbors == 0)
        {
            rm_from_ss_ll(neigh);
            lattice[neigh].state = VACUUM;
        }
    }
}

void add_pvp(int site)
{
    assert(is_vacuum(site) || is_surface_site(site));

    if (is_surface_site(site))
        rm_from_ss_ll(site);

    lattice[site].state = PVP;

    add_to_pvp_ll(site);
}

void rm_pvp(int site)
{
    assert(is_pvp(site));

    rm_from_pvp_ll(site);

    if (lattice[site].neighbors == 0)
        lattice[site].state = VACUUM;
    else
    {
        lattice[site].state = SURFACE_SITE;
        add_to_ss_ll(site);
    }
}

void mv_crystal(int crystal, int site)
{
    rm_crystal(crystal);
    add_crystal(site);
}

int is_crystal(int site)
{
    return lattice[site].state == CRYSTAL;
}

int is_surface_site(int site)
{
    return lattice[site].state == SURFACE_SITE;
}

int is_vacuum(int site)
{
    return lattice[site].state == VACUUM;
}

int is_pvp(int site)
{
    return lattice[site].state == PVP;
}

int is_neigh(int site1, int site2)
{
    int i;

    for (i = 0; i < 12; i++)
        if (lattice[site1].nn[i] == site2)
            return 1;

    return 0;
}

int is_edge(int site)
{
    return lattice[site].edge;
}

void add_to_sc_ll(int site)
{
    ll[site] = crystals_hoc;
    crystals_hoc = site;
    ncrystals++;
}

void rm_from_sc_ll(int site)
{
    int i = crystals_hoc;

    if (i == site)
        crystals_hoc = ll[site];
    else while (ll[i] != site)
        i = ll[i];

    ll[i] = ll[site];
    ll[site] = EMPTY;
    ncrystals--;
}

void add_to_ss_ll(int site)
{
    ll[site] = sites_hoc;
    sites_hoc = site;
    nsites++;
}

void rm_from_ss_ll(int site)
{
    int i = sites_hoc;

    if (i == site)
        sites_hoc = ll[site];
    else while (ll[i] != site)
        i = ll[i];

    ll[i] = ll[site];
    ll[site] = EMPTY;
    nsites--;
}

void add_to_pvp_ll(int site)
{
    ll[site] = pvp_hoc;
    pvp_hoc = site;
    npvp++;
}

void rm_from_pvp_ll(int site)
{
    int i = pvp_hoc;

    if (i == site)
        pvp_hoc = ll[site];
    else while (ll[i] != site)
        i = ll[i];

    ll[i] = ll[site];
    ll[site] = EMPTY;
    npvp--;
}

int rand_ss(void)
{
    assert(nsites != 0);

    int site;
    int i, loops = gsl_rng_uniform_int(rng, nsites);

    for (i = 0, site = sites_hoc; i < loops; i++)
        site = ll[site];

    return site;
}

int nu_crystal_neigh(int site)
{
    int i, count;

    for (i = 0, count = 0; i < 12; i++)
        if (lattice[site].nn[i] != EMPTY && is_crystal(lattice[site].nn[i]))
            count++;

    return count;
}

void check_lattice(void)
{
    int i;

    for (i = 0; i < Nsites; i++)
        assert(nu_crystal_neigh(i) == lattice[i].neighbors);
}

void draw(int max_crystals, int max_pvp)
{
    int i, j;
    int count;

    if (bdraw)
    {
        fprintf(fcoords, "%d\n\n", max_crystals+max_pvp);

        for (i = 0, count = 0; i < Nsites; i++)
            if (is_crystal(i) && lattice[i].neighbors != 12)
            {
                fprintf(fcoords, "C");
                for (j = 0; j < 3; j++)
                    fprintf(fcoords, " %d", lattice[i].pos[j]);
                fprintf(fcoords, "\n");
                count++;

                if (count == max_crystals)
                    break;
            }

        // then stick the rest at the reservoir
        for (i = 0; i < max_crystals - count; i++)
        {
            fprintf(fcoords, "C");
            for (j = 0; j < 3; j++)
                fprintf(fcoords, " %d", center[j]);
            fprintf(fcoords, "\n");
        }

        // draw pvps
        for (i = 0, count = 0; i < Nsites; i++)
        {
            if (is_pvp(i))
            {
                fprintf(fcoords, "F");
                for (j = 0; j < 3; j++)
                    fprintf(fcoords, " %d", lattice[i].pos[j]);
                fprintf(fcoords, "\n");
                count++;
            }
        }

        for (i = 0; i < max_pvp - count; i++)
        {
            fprintf(fcoords, "F");
            for (j = 0; j < 3; j++)
                fprintf(fcoords, " %d", center[j]);
            fprintf(fcoords, "\n");
        }
    }
}

void print_site_info(int site)
{
    int i;

    printf("|| Site #%d Info:\n", site);
    printf("||   State: %d\n", lattice[site].state);
    printf("||   Coord: (%d, %d, %d)\n", lattice[site].pos[0], lattice[site].pos[1], lattice[site].pos[2]);
    printf("||   Crystal Neighbors: %d\n", lattice[site].neighbors);
    printf("||   Neighbors: ");
    for (i = 0; i < 12; i++)
    {
        printf("%d", lattice[site].nn[i]);
        if (i != 11)
            printf(", ");
    }
    printf("\n||   Energy: %.4lf", lattice[site].energy);
    printf("\n||   Rate: %.4le", lattice[site].rate);
    printf("\n");
}

void done(void)
{
    fclose(fcoords);
}
