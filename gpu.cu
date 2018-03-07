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

__global__ void compute_forces_gpu(particle_t * particles, int n, int max_particles_per_bin, int * my_bins, int * my_bins_count, int * bin_neighbors, int * bin_neighbors_count)
{
  // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;
    // 
    particles[tid].ax = particles[tid].ay = 0;
    int bin_number = particles[tid].bin_number;
    // particles from the same bin
    int num_par_same_bin = my_bins_count[bin_number];
    int start_index = max_particles_per_bin * bin_number;
    for(int i = 0; i < num_par_same_bin; i++){
        apply_force_gpu(particles[tid], particles[my_bins[start_index + i]]);
    }
    // particles from the neighboring bins
    int num_nb = bin_neighbors_count[bin_number];
    for(int j = 0 ; j <  num_nb; j++){
        int nb_index = bin_neighbors[8 * bin_number + j];
        int num_par_same_bin = my_bins_count[nb_index];
        int start_index = max_particles_per_bin * nb_index;
        for (int k = 0; k < num_par_same_bin; k++){
            apply_force_gpu(particles[tid], particles[my_bins[start_index + k]]);
        }
        
    }
}

__global__ void set_zero_array(int * my_bins_count, int num_bins)
{
  // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= num_bins) return;
    // 
    my_bins_count[tid] = 0;
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

__global__ void update_bin_number (particle_t * particles, int n, double mysize, int n_row, int n_col, int * my_bins, int * my_bins_count, int max_particles_per_bin){
    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;
    double xloc = particles[tid].x;
    double yloc = particles[tid].y;
    int row = (int) floor(yloc/mysize*n_row);
    int col = (int) floor(xloc/mysize*n_col);
    int bin_index = row + col*n_row;
    particles[tid].bin_number = bin_index;
    int old_count = atomicAdd(&my_bins_count[bin_index], 1);
    my_bins[max_particles_per_bin * bin_index + old_count] = particles[tid].id;
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
        assign all the particles to one of the n_row*n_col bins, my_bins is a 2d array with fixed number of rows
        , and n_row*n_col columns, each column stores the id of particles in that bin
    */
    int n_row = (int) floor(mysize/cutoff);
    int n_col = n_row;
    int num_bins = n_col * n_row;

    // mamimum number of particles in a bin
    int max_particles_per_bin = 5 * (int) floor(n/num_bins); 

    // my_bins is a 2d array with size max_particles_per_bin by num_bins
    int  *my_bins = (int *) malloc( max_particles_per_bin * num_bins * sizeof(int) );

    // my_bins_count is a 1-d array with size 1 by num_bins
    int  *my_bins_count = (int *) malloc(num_bins * sizeof(int) );
    for(int i = 0; i < num_bins; i++){
        my_bins_count[i] = 0;
    }

    double xloc, yloc;
    int row, col;
    int bin_index;
    //std::cout << " line 167"<< std::endl;
    // assign each particle to bins, represented by my_bins and my_bins_count
    for( int i = 0; i < n; i++ )
        {
            // std::cout << "within init particles, i = "<<i << " line 131"<< std::endl;
            xloc = particles[i].x;
            yloc = particles[i].y;
            row = (int)floor(yloc/mysize*n_row);
            col = (int) floor(xloc/mysize*n_col);
            bin_index = row + col * n_row;
            particles[i].bin_number = bin_index;
            my_bins[max_particles_per_bin * bin_index + my_bins_count[bin_index]] = i;
            my_bins_count[bin_index]++;
        }
    //std::cout << " line 181"<< std::endl;
    /*
        use a 2d array of size 8 by num_bins to store the indices of neighbor bins of each bin
    */

    // bin_neighbors is a 2d array with size 8 by num_bins
    int  *bin_neighbors = (int *) malloc( 8 * num_bins * sizeof(int) );

    // bin_neighbors_count is a 1-d array with size 1 by num_bins
    int  *bin_neighbors_count = (int *) malloc(num_bins * sizeof(int) );
    for(int i = 0; i < num_bins; i++){
        bin_neighbors_count[i] = 0;
    }

    for( int j = 0; j < n_col; j++ ){
        for( int i = 0; i < n_row; i++ ){
                 // upper left
                if (i - 1 >= 0 && j - 1 >= 0){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = ((i-1)+(j-1)*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }
                // left
                if (j - 1 >= 0){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = (i+(j-1)*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }
                // lower left
                if (j - 1 >= 0 && i+1 < n_row){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = (i+1+(j-1)*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }
                // up
                if (i - 1 >= 0){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = (i-1+j*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }
                // down
                if (i + 1 < n_row){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = (i+1+j*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }
                // upper right
                if (i - 1 >= 0 && j + 1 < n_col){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = ((i-1)+(j+1)*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }
                // right
                if (j + 1 < n_col){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = (i+(j+1)*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }
                // lower right
                if (i + 1 < n_row && j + 1 < n_col){
                    bin_neighbors[(i+j*n_row) * 8 + bin_neighbors_count[i+j*n_row]] = ((i+1)+(j+1)*n_row);
                    bin_neighbors_count[i+j*n_row]++;
                }

            }
        } 
    //    std::cout << " line 240"<< std::endl;

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    // copy my_bins to GPU
    int *d_my_bins;
    cudaMalloc((void **) &d_my_bins, max_particles_per_bin * num_bins * sizeof(int));

    // copy my_bins_count to GPU
    int *d_my_bins_count;
    cudaMalloc((void **) &d_my_bins_count, num_bins * sizeof(int));

    // copy bin_neighbors to GPU, bin_neighbors need not update
    int *d_bin_neighbors;
    cudaMalloc((void **) &d_bin_neighbors, 8 * num_bins * sizeof(int));

    // copy bin_neighbors_count to GPU, bin_neighbors_count need not update
    int *d_bin_neighbors_count;
    cudaMalloc((void **) &d_bin_neighbors_count, num_bins * sizeof(int));
    //
    cudaMemcpy(d_my_bins, my_bins, max_particles_per_bin * num_bins * sizeof(int), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    //
    cudaMemcpy(d_my_bins_count, my_bins_count, num_bins * sizeof(int), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    //
    cudaMemcpy(d_bin_neighbors, bin_neighbors, 8 * num_bins * sizeof(int), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    //
    cudaMemcpy(d_bin_neighbors_count, bin_neighbors_count, num_bins * sizeof(int), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    //std::cout << " line 278"<< std::endl;
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
        //double  before_compute_forces = read_timer();
        //std::cout << "step "<< step<< " line 290"<< std::endl;
	    int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
	    compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n, max_particles_per_bin, d_my_bins, d_my_bins_count, d_bin_neighbors, d_bin_neighbors_count);
        //double  after_compute_forces = read_timer();
        //std::cout << "compute_forces time " << after_compute_forces - before_compute_forces << std::endl;
        //std::cout << "step "<< step<< " line 295"<< std::endl;
        //
        //  move particles
        //
	    move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);
        //std::cout << "step "<< step<< " line 300"<< std::endl;
        // reset d_my_bins_count
        set_zero_array<<< blks, NUM_THREADS >>> (d_my_bins_count, num_bins);
        //std::cout << "step "<< step<< " line 305"<< std::endl;
        // update particles' bin_number 
        update_bin_number <<< blks, NUM_THREADS >>> (d_particles, n, mysize, n_row, n_col, d_my_bins, d_my_bins_count, max_particles_per_bin);
        //std::cout << "step "<< step<< " line 308"<< std::endl;
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
        // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
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
