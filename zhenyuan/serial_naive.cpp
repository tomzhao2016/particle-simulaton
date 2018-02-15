#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <set>
#include <iostream>

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;
    //

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
    int n_row = 5;
    int n_col = 5;
    std::set<int>  *my_bins = new std::set<int>[n_row*n_col];
    double xloc, yloc;
    int row, col;
    int bin_index;
    std::cout<<"the size of the board is "<<mysize<<std::endl;
    for( int i = 0; i < n; i++ )
        {
            xloc = particles[i].x;
            yloc = particles[i].y;
            row = floor(yloc/mysize*n_row);
            col = floor(xloc/mysize*n_col);
            bin_index = row + col*n_row;
            //std::cout<<typeid(&particles[i]).name()<<std::endl;
            my_bins[bin_index].insert(i);
        }

    // test the size of each bin
    // for( int i = 0; i < n_row*n_col; i++ )
    //     {
    //         std::cout<<(my_bins[i]).size()<<std::endl;
    //     }

    // how to use an iterator
    // std::set<int>::iterator it;
    // std::set<int>::iterator it2;
    // for (it = my_bins[0].begin(); it != my_bins[0].end(); ++it)
    // {   
    //     // for (it2 = my_bins[0].begin(); it2 != my_bins[0].end(); ++it2)
    //     // {
    //     //     std::cout<<particles[it].x - (*it2)->x<<std::endl;
    //     // }
    //     std::cout<<(*it)<<std::endl;
    //     std::cout<<(particles[*it]).x<<std::endl;
    // }
    
    // test using set.clear
    // for (int i = 0; i < n_row*n_col; i++)
    // {
    //     my_bins[0].clear();
    // }
    // std::cout<<"after removing"<<std::endl;
    // std::cout<<(my_bins[0]).size()<<std::endl;
    

    /*
        create an array of sets, each set contains the indices of the neighboring bins of each bin
        
    */
    std::set<int>  *bin_neighbors = new std::set<int>[n_row*n_col];
    for( int i = 0; i < n_row; i++ )
        {
            for( int j = 0; j < n_col; j++ )
            {
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
    //  check bin_neighbors is correct
    // std::set<int>::iterator it;
    // for (int i = 0; i < n_row; i++){
    //     for(int j = 0; j < n_col; j++){
    //         std::cout<<"row "<<i<<" col: "<<j<<std::endl;
    //         for (it = bin_neighbors[i+j*n_row].begin(); it != bin_neighbors[i+j*n_row].end(); ++it)
    //         {   
    //             std::cout<<(*it)<<std::endl;
    //          }
    //     }
    // }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    std::set<int>::iterator it;
    std::set<int>::iterator it2;
    std::set<int>::iterator it3;
    std::set<int>::iterator it4;

    for( int step = 0; step < NSTEPS; step++ )
    {
	    navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute forces, O(N^2) implementation
        //
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }

        /*
            an implementation to compute forces, O(N) implementation

        */
        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                for (it = my_bins[i+j*n_row].begin(); it != my_bins[i+j*n_row].end(); ++it){   
                    // first deal with particles in the same bin
                    for (it2 = my_bins[i+j*n_row].begin(); it2 != my_bins[i+j*n_row].end(); ++it2){   
                        if(&particles[*it] != &particles[*it2]){ // not the same particle
                            apply_force( particles[*it], particles[*it2],&dmin,&davg,&navg);
                        }
                    }
                    // next deal with particels in the neighboring bins
                    for (it3 = bin_neighbors[i+j*n_row].begin(); it3 != bin_neighbors[i+j*n_row].end(); ++it3){   
                        for (it4 = my_bins[*it3].begin(); it4 != my_bins[*it3].end(); ++it4){   
                            apply_force( particles[*it], particles[*it4],&dmin,&davg,&navg);
                        } 
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
    
    /*
        clear the bins and reassign particles to bins

    */
    // first clear all the bins
    for (int i = 0; i < n_row*n_col; i++)
    {
        my_bins[0].clear();
    }
    // then reassign the particles
    for( int i = 0; i < n; i++ )
        {
            xloc = particles[i].x;
            yloc = particles[i].y;
            row = floor(yloc/mysize*n_row);
            col = floor(xloc/mysize*n_col);
            bin_index = row + col*n_row;
            //std::cout<<typeid(&particles[i]).name()<<std::endl;
            my_bins[bin_index].insert(i);
        }
    
    return 0;
}







