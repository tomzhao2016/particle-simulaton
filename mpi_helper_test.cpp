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

	bin_t bins = new bins[1];

	int cur_bin;
	int len_row = 1;
	int len_col = 1;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.clear();
	nidxSet.insert(0);

	find_local_neighbors(bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidxSet);
	std::cout<<" find_local_neighbors pass tests "<<std::endl;
	
}

void findLocalNeighborsTest(){

	bin_t *bins = new bin_t[6];

	int cur_bin;
	int len_row = 2;
	int len_col = 3;
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
	std::cout << " I am here Line 41" << std::endl;
	for (std::set<int>::iterator it = bins[cur_bin].neighbor_idx.begin(); it != bins[cur_bin].neighbor_idx.end(); ++it ){
		std::cout<< *it <<std::endl;
	}

	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 2;
	nidxSet.clear();
	nidxSet.insert(1);
	nidxSet.insert(2);
	nidxSet.insert(4);
	nidxSet.insert(5);
	find_local_neighbors(bins, 2, len_row, len_col);
	std::cout << " I am here Line 51" << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 3;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(3);
	nidxSet.insert(4);
	find_local_neighbors(bins, 3, len_row, len_col);
	std::cout << " I am here Line 61" << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	std::cout<<" find_local_neighbors pass general tests "<<std::endl;

}
void findLocalNeighborsGeneralTest(){

	bin_t *bins = new bin_t[24];

	int cur_bin;
	int len_row = 6;
	int len_col = 4;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.clear();
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(4);
	nidxSet.insert(5);
	std::cout << " I am here Line 105 " << std::endl;
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
	std::cout << " I am here Line 118 " << std::endl;
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
	std::cout << " I am here Line 137 " << std::endl;
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
	std::cout << " I am here Line 152 " << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	cur_bin = 23;
	nidxSet.clear();
	nidxSet.insert(18);
	nidxSet.insert(19);
	nidxSet.insert(22);
	nidxSet.insert(23);
	find_local_neighbors(bins, 23, len_row, len_col);
	std::cout << " I am here Line162" << std::endl;
	assert(bins[cur_bin].neighbor_idx  == nidxSet);

	std::cout<<" find_local_neighbors pass tests "<<std::endl;

}

void main(){
	findLocalNeighborsTest();
	findLocalNeighborsGeneralTest();
}