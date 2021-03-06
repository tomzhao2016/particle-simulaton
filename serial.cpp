#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <set>
#include <iostream>

#define cutoff  0.01

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    double mysize = set_size( n );
    init_particles( n, particles );

    /*
        assign all the particles to one of the n_row*n_col bins, my_bins is a pointer to the 
        n_row * n_col array of hash sets, each set contains the indices of particle_t objects in that bin
    */
    int n_row = (int) floor(mysize/cutoff);
    int n_col = n_row;
    std::set<int>  *my_bins = new std::set<int>[n_row*n_col];
    double xloc, yloc;
    int row, col;
    int bin_index;
    int particle_count = 0;

    /*
        The following 2 things need to be done: each particle needs to be given a bin number,
        and each bin needs to have particles added to it
    */
    for( int i = 0; i < n; i++ )
        {
            xloc = particles[i].x;
            yloc = particles[i].y;
            row = floor(yloc/mysize*n_row);
            col = floor(xloc/mysize*n_col);
            bin_index = row + col*n_row;
            my_bins[bin_index].insert(i);
            particles[i].bin_number = bin_index;
        }

    /*
        create an array of sets, each set contains the indices of the neighboring bins of each bin.
        Since this computation does not change at each time step, we can just do it once. 
        
    */
    std::set<int>  *bin_neighbors = new std::set<int>[n_row*n_col];
    for( int j = 0; j < n_col; j++ ){
        for( int i = 0; i < n_row; i++ ){
                 // upper left
                if (i - 1 >= 0 && j - 1 >= 0){
                    bin_neighbors[i+j*n_row].insert((i-1)+(j-1)*n_row);
                }
                // left
                if (j - 1 >= 0){
                    bin_neighbors[i+j*n_row].insert(i+(j-1)*n_row);
                }
                // lower left
                if (j - 1 >= 0 && i+1 < n_row){
                    bin_neighbors[i+j*n_row].insert(i+1+(j-1)*n_row);
                }
                // up
                if (i - 1 >= 0){
                    bin_neighbors[i+j*n_row].insert(i-1+j*n_row);
                }
                // down
                if (i + 1 < n_row){
                    bin_neighbors[i+j*n_row].insert(i+1+j*n_row);
                }
                // upper right
                if (i - 1 >= 0 && j + 1 < n_col){
                    bin_neighbors[i+j*n_row].insert((i-1)+(j+1)*n_row);
                }
                // right
                if (j + 1 < n_col){
                    bin_neighbors[i+j*n_row].insert(i+(j+1)*n_row);
                }
                // lower right
                if (i + 1 < n_row && j + 1 < n_col){
                    bin_neighbors[i+j*n_row].insert((i+1)+(j+1)*n_row);
                }

            }
        }  


    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    
    std::set<int>::iterator it2;
    std::set<int>::iterator it3;
    std::set<int>::iterator it4;

    for( int step = 0; step < NSTEPS; step++ )
    {   
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        
        /*
            an implementation to compute forces, O(N) implementation
        */
        for (int i = 0; i < n; i++)
        {
            particles[i].ax = particles[i].ay = 0;
            int current_bin_number = particles[i].bin_number;

            // first deal with particles in the same bin
            if(my_bins[current_bin_number].size() > 1)
            {
                for (it2 = my_bins[current_bin_number].begin(); it2 != my_bins[current_bin_number].end(); ++it2)
                {   
                    if(&particles[i] != &particles[*it2])
                    {   // Do not compute forces on the particle itself.
                        apply_force( particles[i], particles[*it2], &dmin, &davg, &navg ); 
                    }
                }
            }

            // next deal with particels in the neighboring bins
            for (it3 = bin_neighbors[current_bin_number].begin(); it3 != bin_neighbors[current_bin_number].end(); ++it3)
            { 
                if(my_bins[*it3].size() > 0)
                {
                    for (it4 = my_bins[*it3].begin(); it4 != my_bins[*it3].end(); ++it4)
                    {   
                        apply_force( particles[i], particles[*it4], &dmin, &davg, &navg);
                    }
                } 
            }

        }

        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );       

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
        
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    /*
        clear the bins and reassign particles to bins
    */
    // first clear all the bins
    for (int i = 0; i < n_row*n_col; i++)
    {
        my_bins[i].clear();
    }
    // then reassign the particles
    for( int i = 0; i < n; i++ )
        {
            xloc = particles[i].x;
            yloc = particles[i].y;
            row = floor(yloc/mysize*n_row);
            col = floor(xloc/mysize*n_col);
            bin_index = row + col*n_row;
            my_bins[bin_index].insert(i);
            particles[i].bin_number = bin_index;
        }
    }

    simulation_time = read_timer( ) - simulation_time;
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    
    return 0;
}
