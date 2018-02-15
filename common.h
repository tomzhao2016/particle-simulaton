#include <set>

#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;


//
// bins
//
typedef struct
{
  std::set<int> particle_idx;
  std::set<int> neighbour_idx;
} bin_t;



//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;

  int cur_bin;

} particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );

int bin_length(int n);

void init_bins(bin_t *b, int n, particle_t *p );
void find_neighbours(bin_t *bins, int cur_bin, int len_bin);
void update_bin(particle_t &p, bin_t *bins, int p_idx);
//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
