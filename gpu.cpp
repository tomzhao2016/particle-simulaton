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
 //    if(tid == 0){
 //    printf("before apply forces \n");
 //    printf("particles[tid].x %f\n",  particles[tid].x);
 //    printf("particles[tid].y %f\n",  particles[tid].y);
 //    printf("particles[tid].ax %f\n",  particles[tid].ax);
 //    printf("particles[tid].ay %f\n",  particles[tid].ay);
 //    printf("particles[tid].vx %f\n",  particles[tid].vx);
 //    printf("particles[tid].vy %f\n",  particles[tid].vy);
 // }
    particles[tid].ax = particles[tid].ay = 0;
    int bin_number = particles[tid].bin_number;
    // particles from the same bin
    int num_par_same_bin = my_bins_count[bin_number];
    int start_index = max_particles_per_bin * bin_number;
    for(int i = 0; i < num_par_same_bin; i++){
        if(particles[tid].id != particles[my_bins[start_index + i]].id){
            apply_force_gpu(particles[tid], particles[my_bins[start_index + i]]);
        }
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
    // for(int j = 0 ; j < n ; j++)
    // apply_force_gpu(particles[tid], particles[j]);


 //    if(tid == 0){
 //    printf("after apply forces \n");
 //    printf("particles[tid].x %f\n",  particles[tid].x);
 //    printf("particles[tid].y %f\n",  particles[tid].y);
 //    printf("particles[tid].ax %f\n",  particles[tid].ax);
 //    printf("particles[tid].ay %f\n",  particles[tid].ay);
 //    printf("particles[tid].vx %f\n",  particles[tid].vx);
 //    printf("particles[tid].vy %f\n",  particles[tid].vy);
 // }
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

 //    if(tid == 0){
 //    printf("before move \n");
 //    printf("particles[tid].x %f\n",  particles[tid].x);
 //    printf("particles[tid].y %f\n",  particles[tid].y);
 //    printf("particles[tid].ax %f\n",  particles[tid].ax);
 //    printf("particles[tid].ay %f\n",  particles[tid].ay);
 //    printf("particles[tid].vx %f\n",  particles[tid].vx);
 //    printf("particles[tid].vy %f\n",  particles[tid].vy);
 // }

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
 //    if(tid == 0){
 //    printf("after move \n");
 //    printf("particles[tid].x %f\n",  particles[tid].x);
 //    printf("particles[tid].y %f\n",  particles[tid].y);
 //    printf("particles[tid].ax %f\n",  particles[tid].ax);
 //    printf("particles[tid].ay %f\n",  particles[tid].ay);
 //    printf("particles[tid].vx %f\n",  particles[tid].vx);
 //    printf("particles[tid].vy %f\n",  particles[tid].vy);
 // }

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
    // atomic add ensures that only one particle can be added to my_bins at a time
    int old_count = atomicAdd(&my_bins_count[bin_index], 1);
    my_bins[max_particles_per_bin * bin_index + old_count] = particles[tid].id;
}

__global__ void test_bin_update (int n, int * my_bins, int * my_bins_count, int num_bins){
    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= 1) return;
    int tmp = 0;
    for (int i = 0; i < num_bins; i++){
        tmp += my_bins_count[i];
    }
    printf("total number of particles in bins %d\n", tmp);
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

    // printf("after init \n");
    // printf("particles[0].x %f\n",  particles[0].x);
    // printf("particles[0].y %f\n",  particles[0].y);
    // printf("particles[0].ax %f\n",  particles[0].ax);
    // printf("particles[0].ay %f\n",  particles[0].ay);
    // printf("particles[0].vx %f\n",  particles[0].vx);
    // printf("particles[0].vy %f\n",  particles[0].vy);

    // printf("after init \n");
    // printf("particles[1].x %f\n",  particles[1].x);
    // printf("particles[1].y %f\n",  particles[1].y);
    // printf("particles[1].ax %f\n",  particles[1].ax);
    // printf("particles[1].ay %f\n",  particles[1].ay);
    // printf("particles[1].vx %f\n",  particles[1].vx);
    // printf("particles[1].vy %f\n",  particles[1].vy);



    /*
        assign all the particles to one of the n_row*n_col bins, my_bins is a 2d array with fixed number of rows
        , and n_row*n_col columns, each column stores the id of particles in that bin
    */
    int n_row = (int) floor(mysize/cutoff);
    int n_col = n_row;
    int num_bins = n_col * n_row;

    // mamimum number of particles in a bin
    int max_particles_per_bin = 20; 

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
    // for(int i = 0; i < n; i++){
    //     std::cout << " id of particle  "<< i << "is"<< particles[i].bin_number<<std::endl;
    // }

    // test my_bins_count add up to n
    // int totol_particle = 0;
    // for(int i = 0; i < num_bins; i++){
    //     totol_particle += my_bins_count[i];
    // }
    // std::cout << "total number of particles is "<<totol_particle<<std::endl;

    // for(int i = 0; i < num_bins; i++){
    //     //std::cout << " # of particle in bin "<< i << "is"<< my_bins_count[i]<<std::endl;
    //     for(int j = 0; j < my_bins_count[i]; j++){
    //         std::cout << " id of particle in bin "<< i << "is"<< my_bins[i*max_particles_per_bin + j]<<std::endl;
    //     }
    // }

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
    // for(int i = 0; i < num_bins; i++){
    //         std::cout << " # of bin_neighbors in bin "<< i << "is"<< bin_neighbors_count[i]<<std::endl;
    // }
    // std::cout << " bin_neighbors in bin 100"<< "is "<< bin_neighbors[100 * 8 + 0]<<std::endl;
    // std::cout << " bin_neighbors in bin 100"<< "is "<< bin_neighbors[100 * 8 + 1]<<std::endl;
    // std::cout << " bin_neighbors in bin 100"<< "is "<< bin_neighbors[100 * 8 + 2]<<std::endl;
    // std::cout << " bin_neighbors in bin 100"<< "is "<< bin_neighbors[100 * 8 + 3]<<std::endl;
    // std::cout << " bin_neighbors in bin 100"<< "is "<< bin_neighbors[100 * 8 + 4]<<std::endl;

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
    //
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //std::cout << " step " <<step << std::endl;

        //
        //  compute forces
        //
	    int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
	    compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n, max_particles_per_bin, d_my_bins, d_my_bins_count, d_bin_neighbors, d_bin_neighbors_count);

        //
        //  move particles
        //
	    move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);

        // reset d_my_bins_count
        int blks2 = (num_bins + NUM_THREADS - 1) / NUM_THREADS;
        set_zero_array<<< blks2, NUM_THREADS >>> (d_my_bins_count, num_bins);

        // update particles' bin_number 
        update_bin_number <<< blks, NUM_THREADS >>> (d_particles, n, mysize, n_row, n_col, d_my_bins, d_my_bins_count, max_particles_per_bin);
        
        // test_bin_update  <<< 1, 1>>> (n, d_my_bins, d_my_bins_count, num_bins);
        
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
