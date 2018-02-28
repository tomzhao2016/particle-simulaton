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

int* get_bin_size(int num_proc_x, int num_proc_y, int rank, int bin_len){
	int idx_col = rank/num_proc_x;
	int idx_row = rank%num_proc_x;
	int *num_bin = new int[2];
	num_bin[0] = (bin_len + num_proc_x - 1)/num_proc_x;
	num_bin[1] = (bin_len + num_proc_y - 1)/num_proc_y;
	if (idx_col == num_proc_y - 1)
		num_bin[1] = bin_len - idx_col*num_bin[1] + 1;
	else if (idx_col == 0)
		num_bin[1] += 1;
	else
		num_bin[1] += 2;

	if (idx_row == num_proc_x - 1)
		num_bin[0] = bin_len - idx_row*num_bin[0] + 1;
	else if (idx_row == 0)
		num_bin[0] += 1;
	else
		num_bin[0] += 2;
	return num_bin;

}


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
//
// This method map a global bin index into a local index inside the proc
//
int glob2loc_row(int global_row, int idx_row, int num_proc_x, int num_bin_row){
	int local_row;
	local_row = global_row - idx_row*num_bin_row;
	if (idx_row > 0)
		local_row++;
	return local_row;
}

int glob2loc_col(int global_col, int idx_col, int num_proc_y, int num_bin_col){
	int local_col;
	local_col = global_col - idx_col*num_bin_col;
	if (idx_col > 0)
		local_col++;
	return local_col;
}

void find_local_neighbors(bin_t *bins, int cur_bin, int len_row, int len_col)
{
    int init_x, init_y, end_x, end_y;
    int bin_x = cur_bin/len_row;
    int bin_y = cur_bin%len_row;  
    if (bin_x == 0) {
        init_x = 0;
        end_x = 2;
    }
    else if(bin_x == len_row - 1) {
        init_x = -1;
        end_x = 1;
    }
    else {
        init_x = -1;
        end_x = 2;
    }
    if (bin_y == 0) {
        init_y = 0;
        end_y = 2;
    }
    else if(bin_y == len_col - 1) {
        end_y = 1;
        init_y = -1;
    }
    else{
        init_y = -1;
        end_y = 2;
    }

    
    for (int i = init_x; i < end_x; i++)
        for (int j = init_y; j < end_y; j++)
            bins[cur_bin].neighbor_idx.insert((bin_x + i)*len_row + bin_y + j);
    
}
//
// initialize bins locally
//
// edge cases: when 2 processors

void init_local_bins(bin_t* local_bins, particle_t* local_particles,int local_size, 
 int *local_bin_size, int num_proc_x, int num_proc_y, int rank, int bin_len){
	
	//
	// col and row index of processor
	//
	int idx_col = rank/num_proc_x;
	int idx_row = rank%num_proc_x;
	
	//
	// num of bins in each proc which has not3 added neighbors yet
	//
	int *num_bin = new int[2];
	num_bin[0] = (bin_len + num_proc_x - 1) /num_proc_x;
	num_bin[1] = (bin_len + num_proc_y - 1)/num_proc_y;
	double cutoff = get_cut_off();
	//
	// assign each particle to bins
	//
	for (int idx = 0; idx < local_size; idx++){
	//for (int idx = 0; idx < 10; idx++){
		
		//
		// global bin index in row and col
		//
		int global_row = (int)floor(local_particles[idx].x/cutoff);
		int global_col = (int)floor(local_particles[idx].y/cutoff);
		//
		// Edge case: if particle is in the left/down most edges then it belongs to the last bin
		//
		if (global_row == bin_len)
			global_row--;
		if (global_col == bin_len)
			global_col--;

		//
		// map into local bin index in row and col
		//
		int local_row = glob2loc_row(global_row, idx_row, num_proc_x, num_bin[0]);
		int local_col = glob2loc_col(global_col, idx_col,  num_proc_y, num_bin[1]);
		//std::cout<<"I am processor "<<rank<<" "<<" I am particle "<<offsets[rank] + idx<<" with local_row and local_col"<<local_row<<" "<<local_col<<std::endl;
		
		//
		// bin idx in 1D array
		//
		int cur_bin = local_row * local_bin_size[0] + local_col;

		// 
		// insert particle into bins
		//
		local_bins[cur_bin].native_particle.insert({idx ,local_particles[idx]});
		std::cout<<"I am processor "<<rank<<" "<<" My native particle is "<<local_bins[cur_bin].native_particle.size()<<std::endl;

	}
	
	int local_col_size = local_bin_size[1];
	int local_row_size = local_bin_size[0];
	for (int i = 0; i<local_col_size*local_row_size; i++){
		local_bins[i].flag = 0;
		find_local_neighbors(local_bins, i, local_row_size, local_col_size);
	}
	if (idx_col == num_proc_y - 1){
		for (int i= 0 ; i< local_row_size;i++){
			local_bins[i].flag = 2;
			local_bins[i + local_row_size].flag = 1;
		}
	}
	else if (idx_col == 0)
		for (int i= 0 ; i < local_row_size;i++){
			local_bins[i + (local_col_size-1)*local_row_size].flag = 2;
			local_bins[i + (local_col_size-2)*local_row_size].flag = 1;
		}
	else
		for (int i= 0 ; i < local_row_size;i++){
			local_bins[i + (local_col_size-1)*local_row_size].flag = 2;
			local_bins[i].flag = 2;
			local_bins[i + local_row_size].flag = 1;
			local_bins[i + (local_col_size-2)*local_row_size].flag = 1;
		}
	


	if (idx_row == num_proc_x - 1)
		for (int i= 0 ; i< local_col_size;i++){
			local_bins[local_row_size*i].flag = 2;
			local_bins[1+local_row_size*i].flag = 1;
		}
	else if (idx_row == 0)
		for (int i= 0 ; i< local_col_size;i++){
			local_bins[local_row_size*i+local_row_size-1].flag = 2;
			local_bins[local_row_size*i+local_row_size-2].flag = 1;
		}
	else
		for (int i = 0 ; i< local_col_size;i++){
			local_bins[local_row_size*i+local_row_size-1].flag = 2;
			local_bins[local_row_size*i].flag = 2;
			local_bins[local_row_size*i+local_row_size-2].flag = 1;
			local_bins[1+local_row_size*i].flag = 1;
		}
	
}


