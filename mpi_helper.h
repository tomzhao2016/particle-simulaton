#include <set>
#include <map>
#include "common2.h"

// #ifndef __CS267_MPI_HELPER_H__
// #define __CS267_MPI_HELPER_H__

// inline int min( int a, int b ) { return a < b ? a : b; }
// inline int max( int a, int b ) { return a > b ? a : b; }



//void get_dest_bin(int curpos, int newpos);          


int get_proc_x(double pos_x, int num_proc_x);
int get_proc_y(double pos_y, int num_proc_y);
int* get_procs(double pos_x, double pos_y, int num_proc_x, int num_proc_y);

int* get_bin_size(int num_proc_x, int num_proc_y, int rank, int bin_len);
int glob2loc_row(int global_row, int idx_row, int num_proc_x, int num_bin_row);
int glob2loc_col(int global_col, int idx_col, int num_proc_y, int num_bin_col);
void init_local_bins(bin_t* local_bins, particle_t* local_particles,int local_size, int *local_bin_size, 
	int num_proc_x, int num_proc_y, int rank, int bin_len);
void find_local_neighbors(bin_t *bins, int cur_bin, int len_row, int len_col);

void clean_local_bins(bin_t *local_bins, int local_bin_size);
int update_local_bins(bin_t *local_bins, std::map<double,particle_t>local_particles_native,
	int *local_bin_size, int num_proc_x, int num_proc_y, int rank, int bin_len);