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
	num_bin[0] = (bin_len + num_proc_x - 1)/num_proc_x;
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
//

void init_local_bins(bin_t* local_bins, particle_t* local_particles, int *local_bin_size, int num_proc_x, int num_proc_y, int rank, int bin_len){
	//
	// num of bins in each processor
	//
	int idx_col = rank/num_proc_x;
	int idx_row = rank%num_proc_y;
	int *num_bin = new int[2];
	num_bin[0] = (bin_len + num_proc_x - 1) /num_proc_x;
	num_bin[1] = (bin_len + num_proc_y - 1)/num_proc_y;

	//
	// assign each particle to bins
	//
	for (int idx = 0; idx < sizeof(local_particles)/sizeof(particle_t); idx++){
		//
		// global bin index in row and col
		//
		int bin_row = (int)floor(local_particle[idx].x/cutoff);
		int bin_col = (int)floor(local_particle[idx].y/cutoff);

		//
		// local bin index in row and col
		//
		int local_row = bin_row%num_bin[0];
		int local_col = bin_col%num_bin[1];

		//
		// bin idx in 1D array
		//
		int cur_bin = local_row * num_bin[0] + local_col;

		// 
		// insert particle into bins, set flag.
		//
		local_bins[cur_bin].particle.insert(local_particles[idx]);
	}

	int local_bin_col = local_bin_size[1];
	int local_bin_row = local_bin_size[0];
	for (int i = 0; i<local_bin_col*local_bin_row; i++){
		local_bins[i].flag = 0;
	}
	if (idx_col == num_proc_y - 1){
		for (int i= 0 ; i< local_bin_row;i++){
			local_bins[i].flag = 2;
			local_bins[i + local_bin_col].flag = 1;
		}
	}
	else if (idx_col == 0)
		for (int i= 0 ; i < local_bin_row;i++){
			local_bins[i + (local_bin_col-1)*local_bin_row].flag = 2;
			local_bins[i + (local_bin_col-2)*local_bin_row].flag = 1;
		}
	else
		for (int i= 0 ; i < local_bin_row;i++){
			local_bins[i + (local_bin_col-1)*local_bin_row].flag = 2;
			local_bins[i].flag = 2;
			local_bins[i + loacl_bin_col].flag = 1;
			local_bins[i + (local_bin_col-2)*local_bin_row].flag = 1;
		}
	


	if (idx_row == num_proc_x - 1)
		for (int i= 0 ; i< local_bin_col;i++){
			local_bins[local_bin_row*i].flag = 2;
			local_bins[1+local_bin_row*i].flag = 1;
		}
	else if (idx_row == 0)
		for (int i= 0 ; i< local_bin_col;i++){
			local_bins[local_bin_row*i+local_bin_col-1].flag = 2;
			local_bins[local_bin_row*i+local_bin_col-2].flag = 1;
		}
	else
		for (int i= 0 ; i< local_bin_col;i++){
			local_bins[local_bin_row*i+local_bin_col-1].flag = 2;
			local_bins[local_bin_row*i].flag = 2;
			local_bins[local_bin_row*i+local_bin_col-2].flag = 1;
			ocal_bins[1+local_bin_row*i].flag = 1;
		}
	find_neighbors(local_bins, idx, len_bin);
}