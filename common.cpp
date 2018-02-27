#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0, numthreads;
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
    
    set_size( n );
    int bin_len = bin_length(n);
    bin_t *bins = new bin_t[bin_len*bin_len]; 

    init_particles( n, particles );
    init_bins(bins, n, particles);
    
    for (int i = 0; i<bin_len ;i++ )
       for (int j = 0;j<bin_len;j++ )
           find_neighbours(bins, i*bin_len+j, bin_len);

    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        
        #pragma omp parallel private(dmin)
        {
            numthreads = omp_get_num_threads();
            dmin = 1.0;
            #pragma omp for reduction(+:navg) reduction(+:davg)
            for( int i = 0; i < n; i++ )
            {
                particles[i].ax = particles[i].ay = 0;
                int cur_bin = particles[i].cur_bin;

                std::set<int> neighbour_idx = bins[cur_bin].neighbour_idx;
                // iterate in 9 neighbourhoods (INCLUDING ITSELF)
                for (std::set<int>::iterator j = neighbour_idx.begin();j != neighbour_idx.end(); ++j){
                    std::set<int> particle_idx = bins[*j].particle_idx;
                    //  iterate over particles inside each bin
                    for(std::set<int>::iterator k = particle_idx.begin();k != particle_idx.end(); ++k){
                        apply_force( particles[i], particles[*k],&dmin,&davg,&navg);
                    }
                }      
            }
            //
            //  move particles
            //
            #pragma omp for
            for( int i = 0; i < n; i++ ) {
                move( particles[i] );
            }  
            #pragma omp master
            for( int i = 0; i < n; i++ ) {
                // update bins after move particles
                update_bin( particles[i], bins, i );
            } 
            if( find_option( argc, argv, "-no" ) == -1 )
            {
              //
              // Computing statistical data
              //
                #pragma omp master
                if (navg) {
                absavg +=  davg/navg;
                nabsavg++;
                }
                #pragma omp critical
                if (dmin < absmin) absmin = dmin;

                //
                //  save if necessary
                //
                #pragma omp master
                if( fsave && (step%SAVEFREQ) == 0 )
                  save( fsave, n, particles );
            }
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
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);
 
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