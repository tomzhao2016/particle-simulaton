#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common2.h"
#include <set>
#include <iostream>



#define cutoff  0.01

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
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
    set_size( n );
    double mysize = get_size();
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
    // std::cout<<"the size of the board is "<<mysize<<std::endl;
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
        create an array of sets, each set contains the indices of the neighboring bins of each bin
        
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

    // for( int step = 0; step < NSTEPS; step++ )
    for( int step = 0; step < NSTEPS; step++ )
    {   
        /*
            an implementation to compute forces, O(N) implementation
        */
        for (int i = 0; i < n; i++){
            particles[i].ax = particles[i].ay = 0;
            int current_bin_number = particles[i].bin_number;
            // first deal with particles in the same bin
            if(my_bins[current_bin_number].size() > 1){
                for (it2 = my_bins[current_bin_number].begin(); it2 != my_bins[current_bin_number].end(); ++it2){   
                    if(&particles[i] != &particles[*it2]){ // not the same particle
                        apply_force( particles[i], particles[*it2]); 
                    }
                }
            }
            // next deal with particels in the neighboring bins
            for (it3 = bin_neighbors[current_bin_number].begin(); it3 != bin_neighbors[current_bin_number].end(); ++it3){ 
                if(my_bins[*it3].size() > 0){
                    for (it4 = my_bins[*it3].begin(); it4 != my_bins[*it3].end(); ++it4){   
                        apply_force( particles[i], particles[*it4]);
                    }
                } 
            }

        }
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );       

          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
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
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time);
    //
    // Clearing space
    // 
    if(particles) 
        free( particles );
    //if(my_bins)
    //   free(my_bins);
    if( fsave )
        fclose( fsave );
    
    
    return 0;
}
