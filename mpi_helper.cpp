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
	//
	// col and row for each processor
	//

	int idx_col = rank/num_proc_x;
	int idx_row = rank%num_proc_x;
	int *num_bin = new int[2];

	//
	// num_bin[0] is the numbers of bins in each processor x axis
	// num_bin[1] is the numbers of bins in each processor y axis
	//
	num_bin[0] = bin_len/num_proc_x;
	num_bin[1] = bin_len/num_proc_y;

	//
	// Note that all processor has same size of bins
	// add neighbor bins, three cases 
	// top: +1
	// bottom: +1
	// middle: +2
	//
	if (idx_col == num_proc_y - 1 && idx_col != 0)
		num_bin[1] += 1;
	else if (idx_col == 0 && idx_col != num_proc_y - 1)
		num_bin[1] += 1;
	else if (idx_col > 0 && idx_col < num_proc_y - 1)
		num_bin[1] += 2;
	

	if (idx_row == num_proc_x - 1 && idx_row != 0)
		num_bin[0] += 1;
	else if (idx_row == 0 && idx_row != num_proc_x - 1)
		num_bin[0] += 1;
	else if(idx_row > 0 && idx_row < num_proc_x - 1)
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

	// std::cout<<"In get procs, pos_x is "<<pos_x<<std::endl;
	// std::cout<<"In get procs, pos_y is "<<pos_y<<std::endl;
	// std::cout<<"In get procs, x_proc is "<<x_proc<<std::endl;
	// std::cout<<"In get procs, y_proc is "<<y_proc<<std::endl;

	// if (pos_x == 0.471405)
	// {
		// std::cout<<"Hello testing";
	// }



	process_ids[index++] = native_proc;

	double bin_len = get_size()/bin_length(num_proc_x, num_proc_y);

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
		process_ids[index++] = (up_proc * num_proc_x) + x_proc + 1;
	}

	// Top Left?
	if ((left_proc != proc_x) && (up_proc != proc_y))
	{
		process_ids[index++] = (up_proc * num_proc_x) + x_proc - 1;
	}

	// Bottom Left?
	if ((left_proc != proc_x) && (down_proc != proc_y))
	{
		process_ids[index++] = (down_proc * num_proc_x) + x_proc - 1;
	}

	// Bottom Right?
	if ((right_proc != proc_x) && (down_proc != proc_y))
	{
		process_ids[index++] = (down_proc * num_proc_x) + x_proc + 1;
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

int glob2loc_row(int global_row, int idx_row, int num_proc_x, int num_bin_row){
	//
	// This method map a global bin index into a local index inside the proc
	// global_row: global row index for this bin before distributed
	// idx_row: row index of the processor
	// num_proc_x: row total number of processors
	// num_bin_row: this should be the native&edge row bin number
	// 

	int local_row;
	local_row = global_row - idx_row*num_bin_row;
	if (idx_row > 0)
		local_row++;
	return local_row;
}

int glob2loc_col(int global_col, int idx_col, int num_proc_y, int num_bin_col){
	//
	// This method map a global bin index into a local index inside the proc
	// global_col: global col index for this bin before distributed
	// idx_col: col index of the processor
	// num_proc_y: col total number of processors
	// num_bin_col: this should be the native&edge col bin number
	//
	int local_col;
	local_col = global_col - idx_col*num_bin_col;
	if (idx_col > 0)
		local_col++;
	return local_col;
}



void find_local_neighbors(bin_t *bins, int cur_bin, int len_row, int len_col)
{
	//
	// this method has same function as find_neighbors in serial code
	// bins: first bins address
	// cur_bin: the target bin, we want to insert neighbor idx to this bin
	// len_row: number of bins in row
	// len_col: number of bins in col
	//
    int init_x, init_y, end_x, end_y;
    int bin_x = cur_bin/len_row;
    int bin_y = cur_bin%len_row;  
    
    //
    // set offset for three cases top/bottom/middle
    // for both x and y 
    // (given cur_bin we just need to know the offsets from this cur_bin)
    // offset -1:2
    // neighbor_idx = cur_bin + offset
    //
    // x most left only 0:2
    if (bin_x == 0) {
        init_x = 0;
        end_x = 2;
    }
    // x most right only -1:1
    else if(bin_x == len_row - 1) {
        init_x = -1;
        end_x = 1;
    }
    // x most right only -1:2
    else {
        init_x = -1;
        end_x = 2;
    }
    // y top only 0:2
    if (bin_y == 0) {
        init_y = 0;
        end_y = 2;
    }
    // y bottom only -1:1
    else if(bin_y == len_col - 1) {
        end_y = 1;
        init_y = -1;
    }
    // x most right only -1:2
    else{
        init_y = -1;
        end_y = 2;
    }

    //
    // find neighbors in x and y direction
    //
    for (int i = init_x; i < end_x; i++)
        for (int j = init_y; j < end_y; j++)
            bins[cur_bin].neighbor_idx.insert((bin_x + i)*len_row + bin_y + j);
    
}
//
// initialize bins locally
//
// edge cases: when 2 processors
//
void init_local_bins(bin_t* local_bins, particle_t* local_particles,int local_size, 
 int *local_bin_size, int num_proc_x, int num_proc_y, int rank, int bin_len){

	//
	// This function intializes bins in each processor.
	// local_bins: is array of bins needed to be initialized
	// local_particles: are array of particles in current processor
	// local_size: is number of particles
	// local_bin_size: is array[2], which is the row and col num of local bin numbers
	// num_proc_x and num_proc_y: are x and y numbers of processors
	// rank: is the id of current processor
	// bin_len is the total number of bins before scattering particles
	//

	//
	// col and row index of processor
	//
	int idx_col = rank/num_proc_x;
	int idx_row = rank%num_proc_x;
	
	//
	// num of bins in each proc which has not added neighbors yet
	// num_bin :the number of native bins + edge bins
	// bin_width: width of each bin
	//
	int *num_bin = new int[2];
	num_bin[0] = bin_len /num_proc_x;
	num_bin[1] = bin_len /num_proc_y;
	double bin_width = (double)get_size()/bin_len;

	//
	// assign each particle to bins
	//
	for (int idx = 0; idx < local_size; idx++){
		
		//
		// global bin index in row and col
		//
		int global_row = (int)floor(local_particles[idx].x/bin_width);
		int global_col = (int)floor(local_particles[idx].y/bin_width);
		//
		// Edge case: if particle is in the left/downside edges then it belongs to the last bin
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
		
		//
		// calculate idx of bin in current processor
		//
		int cur_bin = local_col * local_bin_size[0] + local_row;

		//
		// insert particle into bins
		//
		// hard code since in 2 some particles are wrong

		//if (cur_bin > 0  && cur_bin < local_bin_size[0]*local_bin_size[1]){
			
		//}
		local_bins[cur_bin].native_particle.insert({ local_particles[idx].id,local_particles[idx]});
		// Debug
		// if( rank == 0){
		// 	std::cout<<"This particle x is "<<local_particles[idx].x<<std::endl;
		// 	std::cout<<"This particle y is "<<local_particles[idx].y<<std::endl;
		// 	std::cout<<"This particle local_row is "<<local_row<<std::endl;
		// 	std::cout<<"This particle local_col is "<<local_col<<std::endl;
		// 	std::cout<<"This cur_bin is "<<cur_bin<<std::endl;
		// 	std::cout<<"This idx_row  is "<<idx_row<<std::endl;
		// 	std::cout<<"This idx_col is "<<idx_col<<std::endl;
		// 	std::cout<<"This global_row is "<<global_row<<std::endl;
		// 	std::cout<<"This global_col is "<<global_col<<std::endl;
		// 	std::cout<<"This bin_width is "<<bin_width<<std::endl;
		// 	std::cout<<"This bin_len is "<<bin_len<<std::endl;
		// 	std::cout<<"The size is "<<get_size()<<std::endl;
		// }

	}



	//
	// local col and row bin num,
	//
	int local_col_size = local_bin_size[1];
	int local_row_size = local_bin_size[0];

	//DEBUG - if insert correctly
	// if (rank == 0){
	// 	for (int i = 0; i<local_col_size*local_row_size; i++){
	// 		std::cout<<"I am "<<rank<<" particle size for bins is "<<local_bins[i].native_particle.size()<<std::endl;
	// 	}
		
	// }
	//
	// find neighbors and reset neighbor_index set 
	// and initialize all flags as 0
	//
	for (int i = 0; i<local_col_size*local_row_size; i++){
		local_bins[i].flag = 0;
		find_local_neighbors(local_bins, i, local_row_size, local_col_size);
	}


	//
	// set edge and neighbor flag
	//
	// most down, only up set as 2 and 1

	if(num_proc_x == 1 || num_proc_y == 1){
		// bottom
		if (num_proc_x == 1 && num_proc_y != 1){
			if(idx_col>0 && idx_col<num_proc_y-1){
				for(int i = 0; i<local_row_size;i++){
					local_bins[i].flag = 2;
					local_bins[i + local_row_size].flag = 1;
					local_bins[(local_col_size-1)*local_row_size + i].flag = 2;
					local_bins[(local_col_size-2)*local_row_size + i].flag = 1;
				}
			}
			else if(idx_col == num_proc_y - 1)
				for(int i = 0; i<local_row_size;i++){
					local_bins[i].flag = 2;
					local_bins[i + local_row_size].flag = 1;
				}
			else if(idx_col == 0)
				for(int i = 0; i<local_row_size;i++){
					local_bins[(local_col_size-1)*local_row_size + i].flag = 2;
					local_bins[(local_col_size-2)*local_row_size + i].flag = 1;
				}

		}
	}
	else{
	if (idx_col == num_proc_y - 1){
		for (int i= 0 ; i< local_row_size;i++){
			local_bins[i].flag = 2;
			local_bins[i + local_row_size].flag = 1;
		}
		// most right, from 1 to end
		if (idx_row == num_proc_x - 1)
			for (int i= 1 ; i< local_col_size;i++){
				local_bins[local_row_size*i].flag = 2;
				local_bins[1+local_row_size*i].flag = 1;
			}
		// most left, from 1 to end
		else if (idx_row == 0)
			for (int i= 1 ; i< local_col_size;i++){
				local_bins[local_row_size*i+local_row_size-1].flag = 2;
				local_bins[local_row_size*i+local_row_size-2].flag = 1;
			}
		// if in the middle from 1 to end 
		else
			for (int i = 1 ; i< local_col_size;i++){
				local_bins[local_row_size*i+local_row_size-1].flag = 2;
				local_bins[local_row_size*i].flag = 2;
				local_bins[local_row_size*i+local_row_size-2].flag = 1;
				local_bins[1+local_row_size*i].flag = 1;
			}
	}
	// most up, only bottom
	else if (idx_col == 0){
		for (int i= 0 ; i < local_row_size;i++){
			local_bins[i + (local_col_size-1)*local_row_size].flag = 2;
			local_bins[i + (local_col_size-2)*local_row_size].flag = 1;
		}
		// most right, from 0 to end-1
		if (idx_row == num_proc_x - 1)
			for (int i= 0 ; i< local_col_size-1;i++){
				local_bins[local_row_size*i].flag = 2;
				local_bins[1+local_row_size*i].flag = 1;
			}
		// most left, from 0 to end-1
		else if (idx_row == 0)
			for (int i= 0 ; i< local_col_size-1;i++){
				local_bins[local_row_size*i+local_row_size-1].flag = 2;
				local_bins[local_row_size*i+local_row_size-2].flag = 1;
			}
		// if in the middle from 0 to end -1
		else
			for (int i = 0 ; i< local_col_size-1;i++){
				local_bins[local_row_size*i+local_row_size-1].flag = 2;
				local_bins[local_row_size*i].flag = 2;
				local_bins[local_row_size*i+local_row_size-2].flag = 1;
				local_bins[1+local_row_size*i].flag = 1;
			}

	}
	// if in the middle it should set all the surroundings
	else{
		for (int i= 0 ; i < local_row_size;i++){
			local_bins[i + (local_col_size-1)*local_row_size].flag = 2;
			local_bins[i].flag = 2;

			local_bins[i + local_row_size].flag = 1;
			local_bins[i + (local_col_size-2)*local_row_size].flag = 1;
			
		}
		// most right, from 1 to end-1
		if (idx_row == num_proc_x - 1)
			for (int i= 1 ; i< local_col_size-1;i++){
				local_bins[local_row_size*i].flag = 2;

				local_bins[1+local_row_size*i].flag = 1;
			}
		// most left, from 1 to end-1
		else if (idx_row == 0)
			for (int i= 1 ; i< local_col_size-1;i++){
				local_bins[local_row_size*i+local_row_size-1].flag = 2;
				if(i < local_col_size-1)
					local_bins[local_row_size*i+local_row_size-2].flag = 1;
			}
		// if in the middle from 1 to end -1
		else
			for (int i = 1 ; i< local_col_size-1;i++){
				local_bins[local_row_size*i+local_row_size-1].flag = 2;
				local_bins[local_row_size*i].flag = 2;
				if(i > 0&&i < local_col_size-1){
					local_bins[local_row_size*i+local_row_size-2].flag = 1;
					local_bins[1+local_row_size*i].flag = 1;
				}
			}
		}
	}
	for (int i = 0; i<local_col_size*local_row_size; i++){
		std::cout<<"flag "<<local_bins[i].flag<<std::endl;
	}

	
}

void clean_local_bins(bin_t *local_bins, int local_bin_size){
	//
	// This method cleans the particles in all local_bins
	// local_bin: is array of bins in each processor
	// local_bin_size: is the size of this local_bin
	//
	for (int idx = 0; idx<local_bin_size; idx++){
		local_bins[idx].native_particle.clear();
	}
}

void update_local_bins(bin_t *local_bins, std::map<double,particle_t>local_particles_native_map,
	int *local_bin_size, int num_proc_x, int num_proc_y, int rank, int bin_len){

	//
	// This method assign each particle into bins in this processor
	// local_bins: is empty array of bins needed to be updated
	// local_particles_native_map: are array of new native particles(map) in current processor
	// local_bin_size: is array[2], which is the row and col num of local bin numbers
	// num_proc_x and num_proc_y: are x and y numbers of processors
	// rank: is the id of current processor
	// bin_len is the total number of bins before scattering particles
	//
	// index row of this processor
	//
	int idx_row = rank%num_proc_x;
	int idx_col = rank/num_proc_x;
	
	int local_col_size = local_bin_size[1];
	int local_row_size = local_bin_size[0];
	//
	// nuber of native bins in row and col
	//
	int num_bin_row = bin_len/num_proc_x;
	int num_bin_col = bin_len/num_proc_y;
	//
	// width of each bin
	//
	double bin_width = get_size()/bin_len;
	for (std::map<double, particle_t>::iterator it_p = local_particles_native_map.begin() ;it_p != local_particles_native_map.end(); ++it_p){
		//
		// global index
		//
		int global_row = (int)floor(it_p->second.x/bin_width);
		int global_col = (int)floor(it_p->second.y/bin_width);

		//
		// Edge case
		//
		if (global_row == bin_len){
			global_row--;
		}
		if (global_col == bin_len){
			global_col--;
		}

		//
		// convert to local index
		//
		int local_row = glob2loc_row(global_row, idx_row, num_proc_x, num_bin_row);
		int local_col = glob2loc_col(global_col, idx_col, num_proc_y, num_bin_col);
		//
		// find cur_bin index
		//
		int cur_bin = local_col * local_row_size + local_row;
		//if( rank == 4 ){
			// std::cout<<"This particle is "<<it_p->first<<std::endl;
			// std::cout<<"This particle x is "<<it_p->second.x<<std::endl;
			// std::cout<<"This particle y is "<<it_p->second.y<<std::endl;
			// std::cout<<"This particle local_row is "<<local_row<<std::endl;
			// std::cout<<"This particle local_col is "<<local_col<<std::endl;
			// std::cout<<"This cur_bin is "<<cur_bin<<std::endl;
			// std::cout<<"This idx_row  is "<<idx_row<<std::endl;
			// std::cout<<"This idx_col is "<<idx_col<<std::endl;
			// std::cout<<"This global_row is "<<global_row<<std::endl;
			// std::cout<<"This global_col is "<<global_col<<std::endl;
			// std::cout<<"This bin_width is "<<bin_width<<std::endl;
			// std::cout<<"This bin_len is "<<bin_len<<std::endl;
			// std::cout<<"This local_bin_size is "<<local_row_size*local_col_size<<std::endl;
			// std::cout<<"The size is "<<get_size()<<std::endl;
			std::cout<<"Ima do it using "<<cur_bin<<std::endl;
			local_bins[cur_bin];
			std::cout<<"I did it"<<std::endl;

			local_bins[cur_bin].native_particle.insert({it_p->second.id, it_p->second});
		//}
		
		
	}

}


