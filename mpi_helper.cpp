#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common2.h"



int get_proc_x(double pos_x, int num_proc_x)
{
	// Returns the position of the particle processor along the x-direction
	double len = get_size() / num_proc_x;

	return (int) floor(pos_x / len);
}

int get_proc_y(double pos_y, int num_proc_y)
{
	// Returns the position of the particle processor along the y-direction
	double len = get_size() / num_proc_y;
	return (int) floor(pos_y / len);
}

int* get_bin_size(int num_proc_x, int num_proc_y, int rank, int bin_len){
	int idx_col = rank/num_proc_x;
	int idx_row = rank%num_proc_y;
	int *num_bin = new int[2];
	num_bin[0] = (bin_len + num_proc_x - 1) /num_proc_x;
	num_bin[1] = (bin_len + num_proc_y - 1)/num_proc_y;
	if (idx_col == num_proc_y - 1)
		num_bin[1] = bin_len - idx_col*num_bin[1] + 1；
	else if (idx_col == 0)
		num_bin[1] += 1；
	else
		num_bin[1] += 2;

	if (idx_row == num_proc_x - 1)
		num_bin[0] = bin_len - idx_row*num_bin[0] + 1；
	else if (idx_row == 0)
		num_bin[0] += 1；
	else
		num_bin[0] += 2;
	return num_bin;

}
//
// initialize bins locally
// 1. the global index of bins
//
// bin range: 
// - - - - -
//
init_local_bins(bin_t* local_bins, particle_t* local_particles, int *local_bin_size ,int rank, int bin_len){
	for (int idx = 0; idx < sizeof(local_particles)/sizeof(particle_t); idx++){
		(int)floor(local_particle[idx].x/cutoff)
		local_particle[idx].y
	}
}