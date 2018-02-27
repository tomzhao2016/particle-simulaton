#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common2.h"
#include <iostream>



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

// int* get_bin_size(int num_proc_x, int num_proc_y, int rank, int bin_len){
// 	int idx_col = rank/num_proc_x;
// 	int idx_row = rank%num_proc_y;
// 	int *num_bin = new int[2];
// 	num_bin[0] = (bin_len + num_proc_x - 1)/num_proc_x;
// 	num_bin[1] = (bin_len + num_proc_y - 1)/num_proc_y;
// 	if (idx_col == num_proc_y - 1)
// 		num_bin[1] = bin_len - idx_col*num_bin[1] + 1；
// 	else if (idx_col == 0)
// 		num_bin[1] += 1；
// 	else
// 		num_bin[1] += 2;

// 	if (idx_row == num_proc_x - 1)
// 		num_bin[0] = bin_len - idx_row*num_bin[0] + 1；
// 	else if (idx_row == 0)
// 		num_bin[0] += 1；
// 	else
// 		num_bin[0] += 2;
// 	return num_bin;

// }


int* get_procs(double pos_x, double pos_y, int num_proc_x, int num_proc_y)
{
	// This functions takes in the positions of a particle
	// and the total number of processes along each direction
	// It returns an array of length 9, with all the processes ID
	// that might potentially contain this particle.

	// Note: The returned array might contain elements that are -1,
	// and those elements are just empty. 

	int *process_ids = (int *) malloc(9 * sizeof(int));
	for (int i = 0; i < 9; i++)
	{
		// Initialize the array with all -1s. 
		process_ids[i] = -1;
	}

	int index = 0;

	int x_proc = get_proc_x(pos_x, num_proc_x);
	int y_proc = get_proc_y(pos_y, num_proc_y);
	int proc_x = x_proc;
	int proc_y = y_proc;
	int native_proc = (y_proc * num_proc_x) + x_proc;

	process_ids[index++] = native_proc;

	double bin_len = get_cut_off();

	// Up?
	int up_proc = get_proc_y(pos_y - bin_len, num_proc_y);
	if ((pos_y - bin_len > 0) && (up_proc != proc_y))
	{
		process_ids[index++] = (up_proc * num_proc_x) + x_proc;
	}

	// Down?
	int down_proc = get_proc_y(pos_y + bin_len, num_proc_y);
	if ((pos_y + bin_len < get_size()) && (down_proc != proc_y))
	{
		process_ids[index++] = (down_proc * num_proc_x) + x_proc;
	}

	// Left?
	int left_proc = get_proc_x(pos_x - bin_len, num_proc_x);
	if ((pos_x - bin_len > 0) && (left_proc != proc_x))
	{
		process_ids[index++] = (proc_y * num_proc_x) + left_proc;
	}

	// Right?
	int right_proc = get_proc_x(pos_x + bin_len, num_proc_x);
	if ((pos_x + bin_len < get_size()) && (right_proc != proc_x))
	{
		process_ids[index++] = (proc_y * num_proc_x) + right_proc;
	}

	// Top Right?
	if ((right_proc != proc_x) && (up_proc != proc_y))
	{
		process_ids[index++] = up_proc + 1;
	}

	// Top Left?
	if ((left_proc != proc_x) && (up_proc != proc_y))
	{
		process_ids[index++] = up_proc - 1;
	}

	// Bottom Left?
	if ((left_proc != proc_x) && (down_proc != proc_y))
	{
		process_ids[index++] = down_proc - 1;
	}

	// Top Right?
	if ((right_proc != proc_x) && (down_proc != proc_y))
	{
		process_ids[index++] = down_proc + 1;
	}



	// Debugging
	// std::cout<<"My pos x is "<< pos_x<<" and pos y is "<<pos_y<<std::endl;
	// std::cout<<"The processors I belong to are :::: \n\n\n";
	// for(int i = 0; i < 9; i++)
	// {
	// 	if (process_ids[i] != -1)
	// 		std::cout<<" "<<process_ids[i];
	// }
	// std::cout<<"\n\n\n";
	return process_ids;





}

// int* get_bin_size(int num_proc_x, int num_proc_y, int rank, int bin_len){
// 	int idx_col = rank/num_proc_x;
// 	int idx_row = rank%num_proc_y;
// 	int *num_bin = new int[2];
// 	num_bin[0] = (bin_len + num_proc_x - 1) /num_proc_x;
// 	num_bin[1] = (bin_len + num_proc_y - 1)/num_proc_y;
// 	if (idx_col == num_proc_y - 1)
// 		num_bin[1] = bin_len - idx_col*num_bin[1] + 1；
// 	else if (idx_col == 0)
// 		num_bin[1] += 1；
// 	else
// 		num_bin[1] += 2;

// 	if (idx_row == num_proc_x - 1)
// 		num_bin[0] = bin_len - idx_row*num_bin[0] + 1；
// 	else if (idx_row == 0)
// 		num_bin[0] += 1；
// 	else
// 		num_bin[0] += 2;
// 	return num_bin;

// }

//
// initialize bins locally
// 1. the global index of bins
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

// bin range: 
// - - - - -
//
// init_local_bins(bin_t* local_bins, particle_t* local_particles, int *local_bin_size ,int rank, int bin_len){
// 	for (int idx = 0; idx < sizeof(local_particles)/sizeof(particle_t); idx++){
// 		(int)floor(local_particle[idx].x/cutoff)
// 		local_particle[idx].y
// 	}
// }

