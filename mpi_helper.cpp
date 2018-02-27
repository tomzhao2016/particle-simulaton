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