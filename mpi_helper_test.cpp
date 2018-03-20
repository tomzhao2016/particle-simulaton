#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common2.h"
#include "mpi_helper.h"
#include <iostream>

void findLocalNeighborsSimpleTest(){

	bin_t *bins = new bin_t[1];

	int cur_bin;
	int len_row = 1;
	int len_col = 1;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.clear();
	nidxSet.insert(0);

	find_local_neighbors(bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);
	std::cout<<" find_local_neighbors pass simple tests "<<std::endl;
	
}

void findLocalNeighborsEdgeTest(){

	bin_t *bins = new bin_t[2];

	int cur_bin;
	int len_row = 2;
	int len_col = 1;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	find_local_neighbors(bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 1;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	find_local_neighbors(bins, 1, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);
	std::cout<<" find_local_neighbors other edge tests "<<std::endl;
	
}

void findLocalNeighborsEdge2Test(){

	bin_t *bins = new bin_t[4];

	int cur_bin;
	int len_row = 1;
	int len_col = 4;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	find_local_neighbors(bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 1;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(2);
	find_local_neighbors(bins, 1, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 2;
	nidxSet.clear();
	nidxSet.insert(1);
	nidxSet.insert(2);
	nidxSet.insert(3);
	find_local_neighbors(bins, 2, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 3;
	nidxSet.clear();
	nidxSet.insert(2);
	nidxSet.insert(3);
	find_local_neighbors(bins, 3, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);
	std::cout<<" find_local_neighbors other edge 2 tests "<<std::endl;
	
}

void findLocalNeighborsTest(){

	bin_t *bins = new bin_t[6];

	int cur_bin;
	int len_row = 3;
	int len_col = 2;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(3);
	nidxSet.insert(4);

	find_local_neighbors(bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 1;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(2);
	nidxSet.insert(3);
	nidxSet.insert(4);
	nidxSet.insert(5);
	find_local_neighbors(bins, 1, len_row, len_col);
	//std::cout << " I am here Line 41" << std::endl;
	// for (std::set<int>::iterator it = bins[cur_bin].neighbor_idx.begin(); it != bins[cur_bin].neighbor_idx.end(); ++it ){
	// 	std::cout<< *it <<std::endl;
	// }

	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 2;
	nidxSet.clear();
	nidxSet.insert(1);
	nidxSet.insert(2);
	nidxSet.insert(4);
	nidxSet.insert(5);
	find_local_neighbors(bins, 2, len_row, len_col);
	//std::cout << " I am here Line 51" << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 3;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(3);
	nidxSet.insert(4);
	find_local_neighbors(bins, 3, len_row, len_col);
	//std::cout << " I am here Line 61" << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	std::cout<<" find_local_neighbors pass tests "<<std::endl;

}
void findLocalNeighborsGeneralTest(){

	bin_t *bins = new bin_t[24];

	int cur_bin;
	int len_row = 4;
	int len_col = 6;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(4);
	nidxSet.insert(5);
	//std::cout << " I am here Line 105 " << std::endl;
	find_local_neighbors(bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 7;
	nidxSet.clear();
	nidxSet.insert(2);
	nidxSet.insert(3);
	nidxSet.insert(6);
	nidxSet.insert(7);
	nidxSet.insert(10);
	nidxSet.insert(11);
	find_local_neighbors(bins, 7, len_row, len_col);
	//std::cout << " I am here Line 118 " << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 10;
	nidxSet.clear();
	nidxSet.insert(5);
	nidxSet.insert(6);
	nidxSet.insert(7);
	nidxSet.insert(9);
	nidxSet.insert(10);
	nidxSet.insert(11);
	nidxSet.insert(13);
	nidxSet.insert(14);
	nidxSet.insert(15);
	find_local_neighbors(bins, 10, len_row, len_col);
	//std::cout << " I am here Line 137 " << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 12;
	nidxSet.clear();
	nidxSet.insert(8);
	nidxSet.insert(9);
	nidxSet.insert(12);
	nidxSet.insert(13);
	nidxSet.insert(16);
	nidxSet.insert(17);
	find_local_neighbors(bins, 12, len_row, len_col);
	//std::cout << " I am here Line 152 " << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 23;
	nidxSet.clear();
	nidxSet.insert(18);
	nidxSet.insert(19);
	nidxSet.insert(22);
	nidxSet.insert(23);
	find_local_neighbors(bins, 23, len_row, len_col);
	//std::cout << " I am here Line162" << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	std::cout<<" find_local_neighbors pass general tests "<<std::endl;

}

void getBinSizeTest(){

	int n = 500;
	set_size(n);
	int number_of_processors = 3; 

	int* num_bin;


    int num_proc_x = (int) floor(sqrt(number_of_processors)); // The number of processors along the x-axis.
    int num_proc_y = (int) floor(number_of_processors) / num_proc_x; // The number of processors along the y-axis.

	int rank =0;
	int bin_len = bin_length(num_proc_x, num_proc_y);
	// std::cout<<"bin length is: "<<bin_len<<std::endl;
	num_bin = get_bin_size(num_proc_x, num_proc_y, rank, bin_len);
	//std::cout<<"number of bins in y axis is: "<<num_bin[0]<<std::endl;
	//std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(num_bin[0] == 48);
	// std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(num_bin[1] == 17);


	rank =1;
	// std::cout<<"bin length is: "<<bin_len<<std::endl;
	num_bin = get_bin_size(num_proc_x, num_proc_y, rank, bin_len);
	//std::cout<<"number of bins in y axis is: "<<num_bin[0]<<std::endl;
	//std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(num_bin[0] == 48);
	// std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(num_bin[1] == 18);

	rank =2;
	// std::cout<<"bin length is: "<<bin_len<<std::endl;
	num_bin = get_bin_size(num_proc_x, num_proc_y, rank, bin_len);
	//std::cout<<"number of bins in y axis is: "<<num_bin[0]<<std::endl;
	//std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(num_bin[0] == 48);
	// std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(num_bin[1] == 17);

}


void glob2locTest(){

	int global_row = 35;
	int global_col = 23;
	int idx_row = 0;
	int idx_col = 1;
	int num_proc_x = 1;
	int num_proc_y = 3;
	int num_bin_row = 24;
	int num_bin_col = 16;
	int loc_row = glob2loc_row(global_row, idx_row, num_proc_x, num_bin_row);
	int loc_col = glob2loc_col(global_col, idx_col, num_proc_y, num_bin_col);
	//std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(loc_row == 35);
	// std::cout<<"number of bins in y axis is: "<<num_bin[1]<<std::endl;
	assert(loc_col == 8);
	std::cout<<" Pass glob2loc Test! "<<std::endl;
}

void initLocalBins23Test(){
	// rank = 3 proc_x = 2 proc_y = 3
	bin_t* local_bins = new bin_t[25*18];
	particle_t* local_particles = new particle_t[1];
	int local_size = 0;
	int *local_bin_size = new int[2];
	local_bin_size[0] = 25;
	local_bin_size[1] = 18;
	int num_proc_x = 2;
	int num_proc_y = 3;
	int rank = 3;
	int bin_len = 48;

	init_local_bins(local_bins, local_particles, local_size, 
 local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
	// std::cout<<" I am at 308! "<<std::endl;
	for (int x = 0; x< 25; x++){
		for (int y = 0; y < 18; y++){
			if(x == 0 || y == 0 || y == local_bin_size[1] - 1){
				assert(local_bins[y*25+x].flag==2);
			}
			else if(x == 1){
				assert(local_bins[y*25+x].flag==1);
			}
			else if(y == 1){
				assert(local_bins[y*25+x].flag==1);
			}
			else if(y == 16){
				assert(local_bins[y*25+x].flag==1);
			}
			else{
				assert(local_bins[y*25+x].flag==0);
			}
		}
	}

	local_bins = new bin_t[25*17];
	local_size = 0;
	local_bin_size[0] = 25;
	local_bin_size[1] = 17;
	num_proc_x = 2;
	num_proc_y = 3;
	rank = 0;
	bin_len = 48;

	init_local_bins(local_bins, local_particles, local_size, 
 local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
	//std::cout<<" I am at 308! "<<std::endl;
	for (int x = 0; x< 25; x++){
		for (int y = 0; y < 17; y++){
			if(x == 24 || y == 16){
				assert(local_bins[y*25+x].flag==2);
			}
			else if(x == 23){
				assert(local_bins[y*25+x].flag==1);
			}
			else if(y == 15){
				assert(local_bins[y*25+x].flag==1);
			}
			else{
				assert(local_bins[y*25+x].flag==0);
			}
		}
	}



	std::cout<<" Pass init_local_bins flag 23test "<<std::endl;
}

void initLocalBins22Test(){
	// rank = 3 proc_x = 2 proc_y = 3
	bin_t* local_bins = new bin_t[26*26];
	particle_t* local_particles = new particle_t[1];
	int local_size = 0;
	int *local_bin_size = new int[2];
	local_bin_size[0] = 26;
	local_bin_size[1] = 26;
	int num_proc_x = 2;
	int num_proc_y = 2;
	int rank = 3;
	int bin_len = 50;

	init_local_bins(local_bins, local_particles, local_size, 
 local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
	// std::cout<<" I am at 308! "<<std::endl;
	for (int x = 0; x< 26; x++){
		for (int y = 0; y < 26; y++){
			if(x == 0 || y == 0){
				assert(local_bins[y*26+x].flag==2);
			}
			else if(x == 1){
				assert(local_bins[y*26+x].flag==1);
			}
			else if(y == 1){
				assert(local_bins[y*26+x].flag==1);
			}
			else{
				assert(local_bins[y*26+x].flag==0);
			}
		}
	}

	local_bins = new bin_t[26*26];
	local_size = 0;
	local_bin_size[0] = 26;
	local_bin_size[1] = 26;
	rank = 0;

	init_local_bins(local_bins, local_particles, local_size, 
 local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
	//std::cout<<" I am at 308! "<<std::endl;
	for (int x = 0; x< 26; x++){
		for (int y = 0; y < 26; y++){
			if(x == 25 || y == 25){
				assert(local_bins[y*26+x].flag==2);
			}
			else if(x == 24){
				assert(local_bins[y*26+x].flag==1);
			}
			else if(y == 24){
				assert(local_bins[y*26+x].flag==1);
			}
			else{
				assert(local_bins[y*26+x].flag==0);
			}
		}
	}



	std::cout<<" Pass init_local_bins flag 22test "<<std::endl;
}

void initLocalBins12Test(){
	// rank = 3 proc_x = 2 proc_y = 3
	bin_t* local_bins = new bin_t[50*26];
	particle_t* local_particles = new particle_t[1];
	int local_size = 0;
	int *local_bin_size = new int[2];
	local_bin_size[0] = 50;
	local_bin_size[1] = 26;
	int num_proc_x = 1;
	int num_proc_y = 2;
	int rank = 0;
	int bin_len = 50;

	init_local_bins(local_bins, local_particles, local_size, 
 local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
	// std::cout<<" I am at 308! "<<std::endl;
	for (int x = 0; x< 50; x++){
		for (int y = 0; y < 26; y++){
			if(y == 25){
				assert(local_bins[y*50+x].flag==2);
			}
			else if(y == 24){
				assert(local_bins[y*50+x].flag==1);
			}
			else{
				assert(local_bins[y*50+x].flag==0);
			}
		}
	}

	local_bins = new bin_t[50*26];
	local_size = 0;
	local_bin_size[0] = 50;
	local_bin_size[1] = 26;
	rank = 1;

	init_local_bins(local_bins, local_particles, local_size, 
 local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
	//std::cout<<" I am at 308! "<<std::endl;
	for (int x = 0; x< 50; x++){
		for (int y = 0; y < 26; y++){
			if(y == 0){
				assert(local_bins[y*50+x].flag==2);
			}
			else if(y == 1){
				assert(local_bins[y*50+x].flag==1);
			}
			else{
				assert(local_bins[y*50+x].flag==0);
			}
		}
	}



	std::cout<<" Pass init_local_bins flag 12test "<<std::endl;
}


void findProcNeighborsTest(){

	int rank = 0;
	int num_proc_x = 2;
	int num_proc_y = 3;
	int init_x, init_y, end_x, end_y;
	int num_neighbors = find_proc_neighbors(rank, num_proc_x, num_proc_y, &init_x, &init_y, &end_x, &end_y);
	assert(num_neighbors == 3);
	assert(init_x == 0);
	assert(init_y == 0);
	assert(end_x == 1);
	assert(end_y == 1);

	rank = 1;
	num_neighbors = find_proc_neighbors(rank, num_proc_x, num_proc_y, &init_x, &init_y, &end_x, &end_y);
	assert(num_neighbors == 3);	
	assert(init_x == -1);
	assert(init_y == 0);
	assert(end_x == 0);
	assert(end_y == 1);

	rank = 2;
	num_neighbors = find_proc_neighbors(rank, num_proc_x, num_proc_y, &init_x, &init_y, &end_x, &end_y);
	assert(num_neighbors == 5);
	assert(init_x == 0);
	assert(init_y == -1);
	assert(end_x == 1);
	assert(end_y == 1);

	num_proc_x = 1;
	num_proc_y = 3;
	rank = 0;
	num_neighbors = find_proc_neighbors(rank, num_proc_x, num_proc_y, &init_x, &init_y, &end_x, &end_y);
	assert(num_neighbors == 1);
	assert(init_x == 0);
	assert(init_y == 0);
	assert(end_x == 0);
	assert(end_y == 1);

	rank = 1;
	num_neighbors = find_proc_neighbors(rank, num_proc_x, num_proc_y, &init_x, &init_y, &end_x, &end_y);
	assert(num_neighbors == 2);
	assert(init_x == 0);
	assert(init_y == -1);
	assert(end_x == 0);
	assert(end_y == 1);

	rank = 2;
	num_neighbors = find_proc_neighbors(rank, num_proc_x, num_proc_y, &init_x, &init_y, &end_x, &end_y);
	assert(num_neighbors == 1);
	assert(init_x == 0);
	assert(init_y == -1);
	assert(end_x == 0);
	assert(end_y == 0);
	std::cout<<" Pass find_proc_neighbors test "<<std::endl;

}


void findIdxTest(){
 	// n = 500;
	// int num_proc_x = 2;
	// int num_proc_y = 3;
	int offset_x = -1;
	int offset_y = 0;
	int *local_bin_size = new int[2];
	local_bin_size[0] = 24;
	local_bin_size[1] = 16;
	std::set<int> send_idx = find_idx(offset_x, offset_y, local_bin_size);
	std::set<int> true_set;
	// if proc_x_current = 1, proc_y_current = 0
	for (int i  = 0; i<local_bin_size[1]; i++){
		for(int j = 0; j<local_bin_size[0]; j++){
			if(j == 1&&i<local_bin_size[1]-1)
				true_set.insert(i*local_bin_size[0]+j);
		}

	}
	//assert(send_idx == true_set);

	true_set.clear();
	send_idx = find_idx(offset_x, offset_y, local_bin_size);
	// if proc_x_current = 1, proc_y_current = 1
	for (int i  = 0; i<local_bin_size[1]; i++){
		for(int j = 0; j<local_bin_size[0]; j++){
			if(j == 1&&i<local_bin_size[1]-1&&i>0)
				true_set.insert(i*local_bin_size[0]+j);
		}
	}
	assert(send_idx == true_set);


	true_set.clear();
	send_idx = find_idx(offset_x, offset_y, local_bin_size);
	// if proc_x_current = 1, proc_y_current = 2
	for (int i  = 0; i<local_bin_size[1]; i++){
		for(int j = 0; j<local_bin_size[0]; j++){
			if(j == 1&&i > 0)
				true_set.insert(i*local_bin_size[0]+j);
		}
	}
	assert(send_idx == true_set);
	std::cout<<" Pass find_idx tests "<<std::endl;

}
void main(){

	findLocalNeighborsTest();
	findLocalNeighborsGeneralTest();
	findLocalNeighborsSimpleTest();
	findLocalNeighborsEdgeTest();
	findLocalNeighborsEdge2Test();
	getBinSizeTest();
	glob2locTest();
	initLocalBins23Test();
	initLocalBins22Test();
	initLocalBins12Test();
	findProcNeighborsTest();
	findIdxTest();

}