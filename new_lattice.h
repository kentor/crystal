#define EMPTY -1

enum states {
    VACUUM, CRYSTAL, SURFACE_SITE, PVP
};

struct site {
    int nn[12];             // id of the 12 nearest neighbors. id of -1 is OOB.
    int pos[3];             // the (x,y,z) coordinates of the site
    enum states state;      // possible states defined in the enum
    int neighbors;          // number of crystal neighbors
    int edge;               // boolean for crystal at the edge of the box
    double energy;          // energy required to remove the particle
    double rate;            // rate of removing the particle
};