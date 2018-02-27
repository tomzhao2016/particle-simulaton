#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common2.h"
#include <math.h>
#include "mpi_helper.h"
#include <iostream>

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //
    //  process command line parameters
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


    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // The total number of processes available to us are n_proc.
    int num_proc_x = (int) floor(sqrt(n_proc)); // The number of processors along the x-axis.
    int num_proc_y = (int) ceil(n_proc / sqrt(n_proc)); // The number of processors along the y-axis. 
    printf("Number of total processes %d", n_proc);
    fflush(stdout);
    printf("Number of processes_x %d\n", num_proc_x);
    fflush(stdout);
    printf("Number of processes_y %d\n", num_proc_y);
    fflush(stdout);


    // Create MPI Datatype for particle. 
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    set_size( n );
    if( rank == 0 ){
        printf("Trying to init the particles in rank 0 : \n");
        fflush(stdout);
        init_particles( n, particles );
        printf("Finished init-ing particles in rank 0: \n");
        fflush(stdout);
    }

    std::cout<<"Afer initialization"<<std::endl;

    // sizes array stores the number of particles that each processor is to be sent in the beginning.
    // NOTE: This sizes array will change.  
    int *partition_sizes = (int*) malloc (n_proc * sizeof(int));
    if (rank == 0)
    {
        // Initialize the sizes array to be empty. 
        for (int i = 0; i < n_proc; i++)
            partition_sizes[i] = 0;


        for (int i = 0; i < n; i++)
        {
            // Iterate through all the particles. 
            int x_proc = get_proc_x(particles[i].x, num_proc_x);
            int y_proc = get_proc_y(particles[i].y, num_proc_y);
            int proc_for_p = (y_proc * num_proc_x) + x_proc;
            // Populate the sizes array.
            partition_sizes[proc_for_p] += 1;
        }
    }

    // Send an array of sizes (array of ints) to each processor first. 
    int *local_size; // This is where we will recieve the size. 
    MPI_Scatter(partition_sizes, 1, MPI_INT, local_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Debugging
    if (rank == 1)
    {
        std::cout<<"I am processor 1 \n";
        std::cout<<(*local_size)<<std::endl;
    }
    else if (rank == 2)
    {
        std::cout<<"I am processor 2 \n";
        std::cout<<(*local_size)<<std::endl;
    }

    // Each processor allocates space.
    particle_t *local_particles = (particle_t*) malloc( *local_size * sizeof(particle_t) );


    int *partition_offsets = (int*) malloc( n_proc * sizeof(int) );

    partition_offsets[0] = partition_sizes[0];

    for (int i = 1; i < n_proc; i ++)
    {
        partition_offsets[i] = partition_offsets[i-1] + partition_sizes[i];
    }



    // Recieve the particles for this processor into local_partices
    int *offsets_copy = (int*) malloc (n_proc * sizeof(int));
    for (int i = 0; i < n_proc; i++)
    {
        offsets_copy[i] = partition_offsets[i];
    }

    particle_t *particles_to_scatter = (particle_t*) malloc (n * sizeof(particle_t));

    if (rank == 0)
    {
        for (int i = 0; i < n; i++)
        {
            int x_proc = get_proc_x(particles[i].x, num_proc_x);
            int y_proc = get_proc_y(particles[i].y, num_proc_y);
            int proc_for_p = (y_proc * num_proc_x) + x_proc;

            particles_to_scatter[offsets_copy[proc_for_p] - 1] = particles[i]; 
            offsets_copy[proc_for_p]--;
        }
    }

    

    MPI_Scatterv( particles_to_scatter, partition_sizes, partition_offsets, PARTICLE,
             local_particles, *local_size, PARTICLE, 0, MPI_COMM_WORLD );



    




    int bin_per_proc;
    // std::cout << "TO DO HERE" << std::endl;

    //
    //  allocate storage for local partition
    //

    // TO DO HERE
    // std::cout << "TO DO HERE" << std::endl;

    // int nlocal = partition_sizes[rank];
    // processor_t *local = new mbin_t[nlocal];
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    // TO DO HERE
    // std::cout << "TO DO HERE" << std::endl;

    // fill in init_mbins() in mpi_helper.cpp
    // fill in init_mblocks() in mpi_helper.cpp
    // Or we may use MPI_Comm_split???

    // set_size( n );
    // if( rank == 0 ){
        // init_particles( n, particles );
        // init_mbins(mbins, n, particles); 
    // }
    // MPI_Scatterv( mbins, partition_sizes, partition_offsets, MBIN, local, nlocal, MBIN, 0, MPI_COMM_WORLD );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  1.Update forces
        //
        // for b in native_bins & egde_bins:
        //   for p1 in b:
        //      for p2 in b.neighbors:
        //          apply_force(p1, p2);
        // 
        



        // NOT SURE how to change avg and min
        if( find_option( argc, argv, "-no" ) == -1 )
        {
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); 
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

 
          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        //
        // 2.move particles
        // for b in native_bin & edge_bin:
        //   for p in b:
        //      mpi_move(p, M); 
        // M: map processor_id to particle_t
        // std::map<int, particle_t> M;

        //
        // 3.1 send particle to other processor
        // for processor_id,particle_t in M:
        //   MPI_send(particle_t to native/edge)
        //

        //
        // 3.2receive from other processor
        // for n in neighbor_processor_id:
        //   MPI_receive(particle_t into native/edge)
        //

        //
        // 4.update_bins
        // 

        // barrier

        //
        // 5.1 send edge bins to neighbor processor
        // 

        //
        // 5.2 receive from edge processor
        //

    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    if (partition_offsets)
        free( partition_offsets );
    if (partition_sizes)
        free( partition_sizes );
    // // free( local );
    if (particles)
        free( particles );
    if (local_particles)
        free(local_particles);
    // if( fsave )
        // fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
