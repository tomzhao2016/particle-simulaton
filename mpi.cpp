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
    int num_proc_y = (int) floor(n_proc) / num_proc_x; // The number of processors along the y-axis. 
    printf("Number of total processes %d", n_proc);
    fflush(stdout);
    printf("Number of processes_x %d\n", num_proc_x);
    fflush(stdout);
    printf("Number of processes_y %d\n", num_proc_y);
    fflush(stdout);


    // Create MPI Datatype for particle. 
    MPI_Datatype PARTICLE;
    //MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    set_size( n );
    if( rank == 0 ){
        printf("Trying to init the particles in rank 0 : \n");
        fflush(stdout);
        init_particles( n, particles );
        printf("Finished init-ing particles in rank 0: \n");
        fflush(stdout);
    }

    //std::cout<<"Afer initialization"<<std::endl;

    // sizes array stores the number of particles that each processor is to be sent in the beginning.
    // NOTE: This sizes array will change.  
    int *partition_sizes = (int*) malloc (n_proc * sizeof(int));
    if (rank == 0)
    {
        // Initialize the sizes array to be empty. 
        for (int i = 0; i < n_proc; i++)
            partition_sizes[i] = 0;

        //std::cout<<"Reached line 89 in rank 0"<<std::endl;
        for (int i = 0; i < n; i++) // for each particle
        {
            int *process_ids = (int *) malloc(9 * sizeof(int));
            process_ids = get_procs(particles[i].x, particles[i].y, num_proc_x, num_proc_y);

            for (int j = 0; j < 9; j++)
            {
                if (process_ids[j] != -1)
                {
                    partition_sizes[process_ids[j]] += 1;
                    //std::cout<<"\n\nParticle "<<i<<" goes to processor "<<" "<<process_ids[j];
                }

            }
            free(process_ids);
        }
        //std::cout<<"Reached line 103 in rank 0";
    }



    // Debugging
    if (rank == 0)
    {
        for (int i = 0; i < n_proc; i++)
        {
            //std::cout<<"Partition size "<<i<<" is equal to "<<partition_sizes[i]<<std::endl;
        }
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    // Send an array of sizes (array of ints) to each processor first. 
    int *local_size = (int *)malloc(sizeof(int)); // This is where we will recieve the size. 
    //std::cout<<"Reached line 104 in rank "<<rank<<std::endl;
    MPI_Scatter(partition_sizes, 1, MPI_INT, local_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //std::cout<<"Reached line 111 in rank "<<rank<<std::endl;

    // Debugging
    if (rank == 1)
    {
        //std::cout<<"I am processor 1 \n";
        //std::cout<<(*local_size)<<std::endl;
    }
    else if (rank == 2)
    {
        //std::cout<<"I am processor 2 \n";
        //std::cout<<(*local_size)<<std::endl;
    } 



    // Each processor allocates space.
    particle_t *local_particles = (particle_t*) malloc( *local_size * sizeof(particle_t) );
    //std::cout<<"Rank "<<rank<<" allocated space for local particles"<<std::endl;

    int *partition_offsets = (int*) malloc( n_proc * sizeof(int) );
    int *offsets_copy = (int*) malloc (n_proc * sizeof(int));



    particle_t *particles_to_scatter;


    if (rank == 0)
    {
        //std::cout<<"Reached line 149 in rank 0 "<<std::endl;
        partition_offsets[0] = 0;

        for (int i = 1; i < n_proc; i ++)
        {
            partition_offsets[i] = partition_offsets[i-1] + partition_sizes[i-1];
        }
    

        particles_to_scatter = (particle_t*) malloc ((partition_offsets[n_proc-1] + partition_sizes[n_proc-1] + 1) * sizeof(particle_t));

        // Recieve the particles for this processor into local_partices
        
        for (int i = 0; i < n_proc; i++)
        {
            offsets_copy[i] = partition_offsets[i];
        }
        //std::cout<<"Reached line 166 in rank 0"<<std::endl;
        
        for (int i = 0; i < n; i++)
        {
            int *process_ids = (int *) malloc(9 * sizeof(int));
            process_ids = get_procs(particles[i].x, particles[i].y, num_proc_x, num_proc_y);
            // std::cout<<"Working on line 173."


            for (int j = 0; j < 9; j++)
            {
                if (process_ids[j] != -1)
                {
                    // partition_sizes[process_ids[j]] += 1;
                    //std::cout<<"\n\nParticle "<<i<<" goes to processor "<<" "<<process_ids[j];

                    particles_to_scatter[offsets_copy[process_ids[j]]++] = particles[i];

                }

            }
            free(process_ids);

            // for (int j = 0; j < 9; i++) 
            // {
            //     if (process_ids[j] != -1)
            //     {
            //         std::cout<<"In line 178, i and j are "<<i<<" "<<j<<std::endl;
            //         std::cout<<"In like 179, process_ids[j] and offsets_copy[process_ids[j]] are "<<process_ids[j]<<" "<<offsets_copy[process_ids[j]]<<std::endl;
            //         particles_to_scatter[offsets_copy[process_ids[j]]++] = particles[i];
            //     }
            // }
            // free(process_ids);
        }
        


        //std::cout<<"Reached line 180 in rank 0 "<<std::endl;
    }

    
    MPI_Barrier(MPI_COMM_WORLD); //TODO: Possibly remove this.
    MPI_Scatterv( particles_to_scatter, partition_sizes, partition_offsets, PARTICLE,
             local_particles, *local_size, PARTICLE, 0, MPI_COMM_WORLD );
    
    //
    // bin_len is the total number of bins before divide into local processors
    //
    int bin_len = bin_length(num_proc_x, num_proc_y);
    //
    // local_bin_size is the an array, first elements being row bin number in local processor
    // second elements being col bin number in local processor
    //
    int *local_bin_size = get_bin_size(num_proc_x, num_proc_y, rank, bin_len);
    int local_bin_row = local_bin_size[0];
    int local_bin_col = local_bin_size[1];
    
    //
    // local bin is the bins allocate in each processor
    //
    bin_t *local_bins = new bin_t[local_bin_row*local_bin_col];

    //
    // each bins include particles on left and up edges 
    // and the right most particles belongs to the right most bins
    //
    init_local_bins(local_bins, local_particles, *local_size,
        local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
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
        //  iterate over each native bins
        //

        
        for(int idx = 0; idx < local_bin_row*local_bin_col ; idx++){
            //
            // if flag !=2 it is a native/edge bin
            //
            if (local_bins[idx].flag != 2){
                //  
                // store map of particles in this bin
                //
                std::map<double,particle_t> p1_map = local_bins[idx].native_particle;
                //
                // iterate over all the particles in this map
                //
                // if (rank == 0 && step < 3){
                //    std::cout<<"this bin's particle num is : "<<p1_map.size()<<std::endl;
                // }
                for(std::map<double,particle_t>::iterator p1 = p1_map.begin(); p1!=p1_map.end(); ++p1){
                    //  
                    // store set of neighbor index of this bin
                    //
                    p1->second.ax = p1->second.ay = 0;
                    std::set<int> neighbor_idx = local_bins[idx].neighbor_idx;
                    
                    //
                    // iterate over all the neighbor bins
                    //
                    for(std::set<int>::iterator j = neighbor_idx.begin();j != neighbor_idx.end(); ++j){
                        //  
                        // store map of particles in this neighbor bin
                        //
                        std::map<double,particle_t> p2_map = local_bins[*j].native_particle;

                        // DEBUG - print the neighbor index
                        // if (step == 0 && rank ==0){
                        //     std::cout<<"I am "<<rank<<" bin "<<idx<<std::endl;
                        //     std::cout<<"I am "<<rank<<" My neighbors are "<<*j<<std::endl;
                        // }
                        
                        //
                        // iterate over particles in this bin
                        //
                        for(std::map<double,particle_t>::iterator p2 = p2_map.begin();p2 != p2_map.end(); ++p2){
                            if (p1->first != p2->first){

                                // DEBUG - print acceleration before applying force
                                // if (step == 0){
                                //    std::cout<<"I am "<<rank<<" before p1 acceleration x "<<p1->second.ax<<std::endl;
                                //     std::cout<<"I am "<<rank<<" before p1 acceleration y "<<p1->second.ay<<std::endl;
                                //     std::cout<<"I am "<<rank<<" My x is"<<p1->second.x<<std::endl;
                                //     std::cout<<"I am "<<rank<<" My y is"<<p1->second.y<<std::endl;
                                //     std::cout<<"I am "<<rank<<" My neighbour x is "<<p2->second.x<<std::endl;
                                //     std::cout<<"I am "<<rank<<" My neighbour y is "<<p2->second.y<<std::endl;
                                // }
                                apply_force( p1->second, p2->second,&dmin,&davg,&navg);
                                // DEBUG - print acceleration after applying force
                                // if (step == 0 ){
                                //      std::cout<<"I am "<<rank<<" after p1 acceleration x "<<p1->second.ax<<std::endl;
                                //      std::cout<<"I am "<<rank<<" after p1 acceleration y "<<p1->second.ay<<std::endl;
                                // }
                            }
                        }
                    }
                }   
            }     
        }

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
        // 2.move particles and save all the native particles to a map
        // for b in native_bin & edge_bin:
        //   for p in b:
        //      move(p); 
        //
        std::map<double, particle_t> local_particles_native;
        int local_size_native = 0;
        for(int idx = 0; idx < local_bin_row*local_bin_col ; idx++){
            //
            // if flag !=2 it is a native/edge bin
            //
            if (local_bins[idx].flag != 2){
                //  
                // store map of particles in this bin
                //
                std::map<double,particle_t> p1_map = local_bins[idx].native_particle;
                //
                // iterate over all the particles in this map
                //
                // if (rank == 0 && step < 3){
                //    std::cout<<"this bin's particle num is : "<<p1_map.size()<<std::endl;
                // }
                for(std::map<double,particle_t>::iterator p1 = p1_map.begin(); p1!=p1_map.end(); ++p1){
                    //  
                    // move particles
                    //
                    move(p1->second);
                    local_particles_native.insert({p1.id, p1});
                    local_size_native++;
                }   
            }     
        }

        //
        // 3.1 send and receove particles to/from other processor
        // for processor_id,particle_t in M:
        //   MPI_send(particle_t to native/edge)
        //


        //
        // 4.update_bins
        // local_particles_native and local_size_native contains the new particles and new size 
        // of all particles
        //

        //
        // clean native particles in bins
        //
        clean_local_bins(local_bins, local_bin_row*local_bin_col);
        //
        // reassign each particles in each bins
        //
        update_local_bins(local_bins, local_particles_native, local_size_native);

        // barrier
        MPI_Barrier(MPI_COMM_WORLD);
        
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
