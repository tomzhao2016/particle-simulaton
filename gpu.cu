#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common2.h"
#include <set>
#include <iostream>


#define NUM_THREADS 256

extern double size;
//
//  benchmarking program
//

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__global__ void compute_forces_gpu(particle_t * particles, int * number_neighbors, int * neighbors_indices, int max_neighbors, int n)
{
  // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;
    // 
    particles[tid].ax = particles[tid].ay = 0;
    int temp1 = number_neighbors[tid];

    for(int j = 0 ; j <  temp1; j++)
        apply_force_gpu(particles[tid], particles[number_neighbors[tid*max_neighbors + j]]);
}

__global__ void move_gpu (particle_t * particles, int n, double size)
{

  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t * p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

}



int main( int argc, char **argv )
{    
    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));
    
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
            // std::cout << "within init particles, i = "<<i << " line 131"<< std::endl;
            xloc = particles[i].x;
            yloc = particles[i].y;
            row = (int)floor(yloc/mysize*n_row);
            col = (int) floor(xloc/mysize*n_col);
            bin_index = row + col*n_row;
            // std::cout << "within init particles, i = "<<i << " xloc "<<xloc<< std::endl;
            // std::cout << "within init particles, i = "<<i << " yloc "<<yloc<< std::endl;
            // std::cout << "within init particles, i = "<<i << " row "<<row<< std::endl;
            // std::cout << "within init particles, i = "<<i << " col "<<col<< std::endl;
            // std::cout << "within init particles, i = "<<i << " mysize "<<mysize<< std::endl;
            // std::cout << "within init particles, i = "<<i << " nrow "<<n_row<< std::endl;
            // std::cout << "within init particles, i = "<<i << " bin_index"<<bin_index<< std::endl;
            my_bins[bin_index].insert(i);
            // std::cout << "within init particles, i = "<<i << " line 139"<< std::endl;
            particles[i].bin_number = bin_index;
            // std::cout << "within init particles, i = "<<i << " line 141"<< std::endl;
        }
    /*
        create an array of sets, each set contains the indices of the neighboring bins of each bin,
        this doesn't need to be updated through the simulation
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


    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );
    
    std::set<int>::iterator it2;
    std::set<int>::iterator it3;
    std::set<int>::iterator it4;

    // density is constant, use a factor of 3
    int max_neighbors = (int) 3 * n/n_row/n_col;

    // a 1-D array of size n that contains the number of effective neighbors
    int *number_neighbors = (int *) malloc(  n * sizeof(int) );

    // a 2-D array of size max_neighbors * n that contains the number of effective neighbors of n particles
    int length_neighbors_indices = max_neighbors * n;
    int *neighbors_indices = (int *) malloc(  length_neighbors_indices * sizeof(int) );


    for( int step = 0; step < NSTEPS; step++ )
    {
        /*
        store the indices of neighboring nodes in a 2-D array
        */
        // std::cout << "step "<< step<< "before convert bins to arrays" << std::endl;
        for(int i = 0; i < n; i++){
            int index_temp = 0;
            int current_bin_number = particles[i].bin_number;
            // first deal with particles in the same bin
            for(it2 = my_bins[current_bin_number].begin(); it2 != my_bins[current_bin_number].end(); ++it2){   
                if(particles[i].id != particles[*it2].id){ // not the same particle
                    neighbors_indices[i * max_neighbors + index_temp] = *it2;
                    index_temp++;
                }
            }
            // next deal with particels in the neighboring bins
            for (it3 = bin_neighbors[current_bin_number].begin(); it3 != bin_neighbors[current_bin_number].end(); ++it3){ 
                for (it4 = my_bins[*it3].begin(); it4 != my_bins[*it3].end(); ++it4){   
                        neighbors_indices[i * max_neighbors + index_temp] = *it4;
                        index_temp++;
                    }
            }
            number_neighbors[i] = index_temp;
        }
        if(step == 0){
            std::cout << number_neighbors << std::endl;
        }

        // copy neighbors_indices, number_neighbors to GPU for fast access
        int *d_number_neighbors;
        cudaMalloc((void **) &d_number_neighbors, n * sizeof(int));
        // 
        int *d_neighbors_indices;
        cudaMalloc((void **) &d_neighbors_indices, length_neighbors_indices * sizeof(int));
        // 
        cudaMemcpy(d_number_neighbors, number_neighbors, n * sizeof(int), cudaMemcpyHostToDevice);
        cudaThreadSynchronize();
        //
        cudaMemcpy(d_neighbors_indices, neighbors_indices, length_neighbors_indices * sizeof(int), cudaMemcpyHostToDevice);
        cudaThreadSynchronize();


        //
        //  compute forces
        //

	    int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
	    compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, d_number_neighbors, d_neighbors_indices, max_neighbors, n);
        
        //
        //  move particles
        //
	    move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);
        
        // copy back to memory in each step, so as to update my_bins
        cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
        cudaThreadSynchronize();
        /*
        clear the bins and reassign particles to bins
        */
        // first clear all the bins
        for (int i = 0; i < n_row*n_col; i++)
        {
            my_bins[i].clear();
        }
        // then reassign the particles
        for( int i = 0; i < n; i++ ){
                xloc = particles[i].x;
                yloc = particles[i].y;
                row = floor(yloc/mysize*n_row);
                col = floor(xloc/mysize*n_col);
                bin_index = row + col*n_row;
                my_bins[bin_index].insert(i);
                particles[i].bin_number = bin_index;
            }

        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
        // Copy the particles back to the CPU
            save( fsave, n, particles);
        }
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
