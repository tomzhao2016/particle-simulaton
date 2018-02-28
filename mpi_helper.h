#include <set>
#include <map>
#include "common2.h"

// #ifndef __CS267_MPI_HELPER_H__
// #define __CS267_MPI_HELPER_H__

// inline int min( int a, int b ) { return a < b ? a : b; }
// inline int max( int a, int b ) { return a > b ? a : b; }



// typedef struct{
// 	std::map<int, mbin_t> native_bins;
// 	std::map<int, mbin_t> neighbor_bins;
// 	std::map<int, mbin_t> edge_bins;
// } processor_t;

//void get_dest_bin(int curpos, int newpos);

//void init_mbins(mbin_t mbins, int n, particle_t* particles);

int get_proc_x(double pos_x, int num_proc_x);
int get_proc_y(double pos_y, int num_proc_y);
int* get_procs(double pos_x, double pos_y, int num_proc_x, int num_proc_y);

int* get_bin_size(int num_proc_x, int num_proc_y, int rank, int bin_len);
int glob2loc_row(int global_row, int idx_row, int num_proc_x, int num_bin_row);
int glob2loc_col(int global_col, int idx_col, int num_proc_y, int num_bin_col);
void init_local_bins(bin_t* local_bins, particle_t* local_particles,int local_size, int *local_bin_size,int *offsets, int num_proc_x, int num_proc_y, int rank, int bin_len);


