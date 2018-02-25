#include <set>
#include <map>
#include "common.h"

#ifndef __CS267_MPI_HELPER_H__
#define __CS267_MPI_HELPER_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

typedef struct{
	std::map<int, particle_t> native_particle;
	std::set<int> neighbor_idx;
} mbin_t;

typedef struct{
	std::map<int, mbin_t> native_bins;
	std::map<int, mbin_t> neighbor_bins;
	std::map<int, mbin_t> edge_bins;
} processor_t;

void get_dest_bin(int curpos, int newpos);

void init_mbins(mbin_t mbins, int n, particle_t* particles);



