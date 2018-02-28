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
#include <algorithm>


double size;
//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

double get_size ()
{
    return size;
}

double get_cut_off ()
{
    return cutoff;
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    printf("I am in line 1 of init particles.... \n");
    srand48( time( NULL ) );
    printf("I am trying to init particles.... \n");
    fflush(stdout);
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //

        // printf("I am before, %d\n", i);
        // fflush(stdout);
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);


        // std::cout<<"I am after the first reference to p\n"<<std::endl;
        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    // printf("I am before free shuffle \n");
    // fflush(stdout);
    free( shuffle );
    // printf("Now I have freed shuffle \n");
    // fflush(stdout);
    // TODO: FREE THIS LATER. 

}


int bin_length(int num_proc_x, int num_proc_y)
{
    //
    // smallest gcd
    //
    int s_gcd = std::__gcd(num_proc_x, num_proc_y);
    //
    // original bin size, which may not divided by num_pro_x and num_proc_y
    //
     int bin_len = (int)ceil(size/cutoff);
     bin_len = (int)floor(bin_len/s_gcd) * s_gcd;
     return bin_len;
}

//
// Initialize the bins and assign each particle into bins
//
// void init_bins(bin_t *bins, int n, particle_t *p )
// {
//     // cutoff should dividable by 1
//     int len_bin = (int)ceil(size/cutoff);
//     long int num_bin = len_bin * len_bin;
//     int sum = 0;

//     //bin_t* bins = (bin_t*) malloc( num_bin * sizeof(bin_t) );

    
//     for (int i = 0; i < n; i++ )
//     {
//         int bin_x = (int)floor(p[i].x/cutoff);
//         int bin_y = (int)floor(p[i].y/cutoff);

//         bins[bin_x * len_bin + bin_y].particle_idx.insert(i);
//         // 
//         // assign each particle into a bin
//         //
//         p[i].cur_bin = bin_x * len_bin + bin_y;
//     }
       


    

// }

//
// find neighbours of current bin 
// return the bin idx as a set, including itself
//
void find_neighbors(bin_t *bins, int cur_bin, int len_bin)
{
    int init_x, init_y, end_x, end_y;
    int bin_x = cur_bin/len_bin;
    int bin_y = cur_bin%len_bin;  
    if (bin_x == 0) {
        init_x = 0;
        end_x = 2;
    }
    else if(bin_x == len_bin - 1) {
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
    else if(bin_y == len_bin - 1) {
        end_y = 1;
        init_y = -1;
    }
    else{
        init_y = -1;
        end_y = 2;
    }

    
    for (int i = init_x; i < end_x; i++)
        for (int j = init_y; j < end_y; j++)
            bins[bin_x * len_bin + bin_y].neighbor_idx.insert((bin_x + i)*len_bin + bin_y + j);

    
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
	if (r2 != 0)
        {
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	      *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
	
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//
//
// void update_bin(particle_t &p, bin_t *bins, int p_idx)
// {
//     int len_bin = (int)ceil(size/cutoff);
//     int bin_x = (int)floor(p.x/cutoff);
//     int bin_y = (int)floor(p.y/cutoff);

//     bins[p.cur_bin].particle_idx.erase(p_idx);
//     bins[bin_x * len_bin + bin_y].particle_idx.insert(p_idx);

//     p.cur_bin = bin_x * len_bin + bin_y;
// }
//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
