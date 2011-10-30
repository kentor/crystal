#define EMPTY -1

enum states {
    VACUUM, CRYSTAL, SURFACE_SITE, PVP
};

struct site {
    int nn[12];             // id of the 12 nearest neighbors or -1 if not a valid neighbor
    int pos[3];             // the (x, y, z) coordinates of the site
    enum states state;      // possible states defined in the enumerator
    int neighbors;          // how many crystal neighbors
    int edge;
    double energy;
    double bonus_energy;
    double rate;
};

void init_lattice(int argc, char **argv);

void seed_crystal_cube(int radius);

void add_crystal(int site);
void rm_crystal(int site);
void mv_crystal(int crystal, int site);
void add_pvp(int site);
void rm_pvp(int site);

int is_crystal(int site);
int is_surface_site(int site);
int is_vacuum(int site);
int is_neigh(int site1, int site2);
int is_edge(int site);

void add_to_sc_ll(int site);
void rm_from_sc_ll(int site);
void add_to_ss_ll(int site);
void rm_from_ss_ll(int site);
void add_to_pvp_ll(int site);
void rm_from_pvp_ll(int site);

int rand_ss(void);

int rand_crystal_neigh(int site);
int nu_crystal_neigh(int site);
int nu_pvp_neigh(int site);

void draw(int max_crystals, int max_pvp);
void print_site_info(int site);
void done(void);
