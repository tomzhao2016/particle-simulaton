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

void findLocalNeighborsTest(){

	bin_t *bins = (bin_t*)malloc(6 * sizeof(bin_t));
	int cur_bin;
	int len_row = 2;
	int len_col = 3;
	std::set<int> nidxSet; 

	cur_bin = 0;
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(3);
	nidxSet.insert(4);
	find_local_neighbors(bin_t *bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidx);

	cur_bin = 1;
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(2);
	nidxSet.insert(3);
	nidxSet.insert(4);
	nidxSet.insert(5);
	find_local_neighbors(bin_t *bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidx);

	cur_bin = 2;
	nidxSet.insert(1);
	nidxSet.insert(2);
	nidxSet.insert(4);
	nidxSet.insert(5);
	find_local_neighbors(bin_t *bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidx);

	cur_bin = 3;
	nidxSet.insert(0);
	nidxSet.insert(1);
	nidxSet.insert(3);
	nidxSet.insert(4);
	find_local_neighbors(bin_t *bins, 0, len_row, len_col);
	assert(bins[cur_bin].neighbor_idx  == nidx);

	std::cout<<" find_local_neighbors pass tests "<<std::endl;

}
void main(){
	findLocalNeighborsTest();
}