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
	return (int) floor(pos_x / num_proc_x);
}

int get_proc_y(double pos_y, int num_proc_y)
{
	// Returns the position of the particle processor along the y-direction
	return (int) floor(pos_y / num_proc_y);
}