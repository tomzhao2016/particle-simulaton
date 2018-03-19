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

void initLocalBinsTest(){
	// rank = 3 proc_x = 2 proc_y = 3
	bin_t* local_bins = new bin_t[25*18];
	particle_t* local_particles = new particle_t[1];
	int local_size = 0;
	int *local_bin_size;
	local_bin_size[0] = 25;
	local_bin_size[1] = 18;
	int num_proc_x = 2;
	int num_proc_y = 3;
	int rank = 3;
	int bin_len = 48;

	init_local_bins(local_bins, local_particles, local_size, 
 local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
	std::cout<<" I am at 308! "<<std::endl;
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
}

void main(){

	findLocalNeighborsTest();
	findLocalNeighborsGeneralTest();
	findLocalNeighborsSimpleTest();
	findLocalNeighborsEdgeTest();
	findLocalNeighborsEdge2Test();
	getBinSizeTest();
	glob2locTest();
	initLocalBinsTest();

}