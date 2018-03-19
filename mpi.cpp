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


    /*********************************************
     * Set up MPI
     ********************************************/
    int number_of_processors, rank;
    MPI_Request send_request0,send_request1,
    send_request2, send_request3, send_request4, 
    send_request5, send_request6, send_request7,
    recv_request0, recv_request1, recv_request2, recv_request3,
    recv_request4, recv_request5, recv_request6, recv_request7;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &number_of_processors );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    /******************************************
     *  Allocate Generic Resources
     ********************************************/
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    /************************************************************************
     * Calculate total number of processors available, and how to split them into a grid.
     ************************************************************************************/
    // The total number of processes available to us are n_proc.
    int num_proc_x = (int) floor(sqrt(number_of_processors)); // The number of processors along the x-axis.
    int num_proc_y = (int) floor(number_of_processors) / num_proc_x; // The number of processors along the y-axis.


    /*****************************************
     * Create MPI Datatype
     ********************************************/
    MPI_Datatype PARTICLE;
    //MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    MPI::COMM_WORLD.Set_errhandler ( MPI::ERRORS_THROW_EXCEPTIONS );

    set_size( n );
    if( rank == 0 )
    {
        init_particles( n, particles );
    }

    /**********************************************************************************************
     * Initialization ---- sending the particles from root node to the oher node */
    /****************************************************************************************/
    // sizes array stores the number of particles that each processor is to be sent in the beginning.
    int *partition_sizes = (int*) malloc (number_of_processors * sizeof(int));
    if (rank == 0)
    {
        // Initialize the sizes array to be empty. 
        for (int i = 0; i < number_of_processors; i++)
            partition_sizes[i] = 0;

        for (int i = 0; i < n; i++) // for each particle
        {
            int *process_ids = (int *) malloc(9 * sizeof(int));
            process_ids = get_procs(particles[i].x, particles[i].y, num_proc_x, num_proc_y);

            for (int j = 0; j < 9; j++)
            {
                if (process_ids[j] != -1)
                {
                    partition_sizes[process_ids[j]] += 1;
                }

            }
            free(process_ids);
        }
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    // Send an array of sizes (array of ints) to each processor first. 
    int *local_size = (int *)malloc(sizeof(int)); // This is where we will recieve the size. 
    MPI_Scatter(partition_sizes, 1, MPI_INT, local_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Each processor allocates space.
    particle_t *local_particles = (particle_t*) malloc( *local_size * sizeof(particle_t) );
    //std::cout<<"Rank "<<rank<<" allocated space for local particles"<<std::endl;

    int *partition_offsets = (int*) malloc( number_of_processors * sizeof(int) );
    int *offsets_copy = (int*) malloc (number_of_processors * sizeof(int));

    // Is initialized later.
    particle_t *particles_to_scatter;

    if (rank == 0)
    {
        partition_offsets[0] = 0;
        for (int i = 1; i < number_of_processors; i ++)
        {
            partition_offsets[i] = partition_offsets[i-1] + partition_sizes[i-1];
        }

        particles_to_scatter = (particle_t*) malloc ((partition_offsets[number_of_processors-1] + partition_sizes[number_of_processors-1] + 1) * sizeof(particle_t));
        // Recieve the particles for this processor into local_partices
        
        for (int i = 0; i < number_of_processors; i++)
        {
            offsets_copy[i] = partition_offsets[i];
        }
        
        for (int i = 0; i < n; i++)
        {
            int *process_ids = (int *) malloc(9 * sizeof(int));
            process_ids = get_procs(particles[i].x, particles[i].y, num_proc_x, num_proc_y);

            for (int j = 0; j < 9; j++)
            {
                if (process_ids[j] != -1)
                {
                    particles_to_scatter[offsets_copy[process_ids[j]]++] = particles[i];
                }

            }
            free(process_ids);
        }
    }

    
    MPI_Barrier(MPI_COMM_WORLD); //TODO: Possibly remove this.
    MPI_Scatterv( particles_to_scatter, partition_sizes, partition_offsets, PARTICLE,
             local_particles, *local_size, PARTICLE, 0, MPI_COMM_WORLD );
    

    /******************************************
     * Get current indices of the processors
     ******************************************/

    // The current indices of the processors.
    int proc_y_current = rank/num_proc_x;
    int proc_x_current = rank%num_proc_x;


    // The new processor of the particle.
    int proc_x_new, proc_y_new;

    /*************************************************************
     * Variable initialization
     *********************************************************/



    //
    // bin_len is the total number of bins
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

    /**********************************************************************
     * Simulation begins
     * ******************************************************************
     */
    double simulation_time = read_timer();
    for( int step = 0; step < NSTEPS; step++ ) {


        MPI_Barrier(MPI_COMM_WORLD);
        navg = 0;
        dmin = 1.0;
        davg = 0.0;


        // std::cout<<"I am beginning of step: "<<step<<std::endl;

        //std::cout << "I am beginning of step, in rank " << step << " " << rank << std::endl;

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if (find_option(argc, argv, "-no") == -1)
            if (fsave && (step % SAVEFREQ) == 0)
                save(fsave, n, particles);

        //
        //  1.Update forces 
        //  iterate over each native bins
        //
        /*********************************************************************
         * Update forces
         ****************************************************************/
        for (int idx = 0; idx < local_bin_row * local_bin_col; idx++) {
            // if flag !=2 it is a native/edge bin
            if (local_bins[idx].flag != 2) {
                //  
                // store map of particles in this bin
                //
                std::map<double, particle_t> p1_map = local_bins[idx].native_particle;
                //
                // iterate over all the particles in this map

                for (std::map<double, particle_t>::iterator p1 = p1_map.begin(); p1 != p1_map.end(); ++p1) {
                    //  
                    // store set of neighbor index of this bin
                    //
                    p1->second.ax = p1->second.ay = 0;
                    std::set<int> neighbor_idx = local_bins[idx].neighbor_idx;

                    //
                    // iterate over all the neighbor bins
                    //
                    for (std::set<int>::iterator j = neighbor_idx.begin(); j != neighbor_idx.end(); ++j) {
                        //  
                        // store map of particles in this neighbor bin
                        //
                        std::map<double, particle_t> p2_map = local_bins[*j].native_particle;

                        //
                        // iterate over particles in this bin
                        //
                        for (std::map<double, particle_t>::iterator p2 = p2_map.begin(); p2 != p2_map.end(); ++p2) {
                            if (p1->first != p2->first) {
                                apply_force(p1->second, p2->second, &dmin, &davg, &navg);
                            }
                        }
                    }
                }
            }
        }

        /****************************
         * Statistical data
         *************************/
        // NOT SURE how to change avg and min
        if (find_option(argc, argv, "-no") == -1) {

            MPI_Reduce(&davg, &rdavg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&navg, &rnavg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&dmin, &rdmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);


            if (rank == 0) {
                //
                // Computing statistical data
                //
                if (rnavg) {
                    absavg += rdavg / rnavg;
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

        /******************************
         * Move particles **********
         ******************************/

        std::map<double, particle_t> local_particles_native_map;
        // Move the particles.
        for (int idx = 0; idx < local_bin_row * local_bin_col; idx++) {

            //
            // if flag !=2 it is a native/edge bin
            //
            if (local_bins[idx].flag != 2) {

                // store map of particles in this bin
                //
                std::map<double, particle_t> p1_map = local_bins[idx].native_particle;
                //
                // iterate over all the particles in this map
                //
                for (std::map<double, particle_t>::iterator p1 = p1_map.begin(); p1 != p1_map.end(); ++p1) {
                    // move particles
                    // std::cout<<"Old Info, Step, ID, X, Y "<<step<<" "<<p1->second.id<<" "<<p1->second.x<<" "<<p1->second.y<<std::endl;
                    move(p1->second);
                    local_bins[idx].native_particle.erase(p1->first);

                    int proc_x_next = get_proc_x(p1->second.x, num_proc_x);
                    int proc_y_next = get_proc_y(p1->second.y, num_proc_y);

                    if(proc_x_next != proc_x_current || proc_y_next != proc_y_current){

                        int target = num_proc_x * proc_y_next + proc_x_next;
                        // assume all go to neighbors
                        if(abs(proc_x_next - proc_x_current)<=1 &&abs(proc_y_next-proc_y_current)<=1){
                            MPI_Request request;
                            MPI_Isend(&p1->second ,1 , PARTICLE, target, target, MPI_COMM_WORLD, &request);
                        }
                        continue;             
                    }
                    addto_local_bins(local_bins, p1->second, local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
                }
            }
        }

        // Tell neighbors finish sending
        particle_t end;
        end.id = -1;
        int init_x, init_y;
        int end_x, end_y;
        int num_neighbors = find_proc_neighbors(rank, num_proc_x, num_proc_y, &init_x, &init_y, &end_x, &end_y);
        MPI_Request request;
        for (int offset_x = init_x; offset_x <= end_x; offset_x++){
            for(int offset_y = init_y; offset_y <= end_y; offset_y++){
                int send_to_idx = (proc_y_current + offset_y) * num_proc_x + proc_x_current + offset_x;
                if (offset_x == 0 && offset_y == 0) continue;
                else MPI_Isend(&end, 1, PARTICLE, send_to_idx, rank, MPI_COMM_WORLD, &request);
            }
        }

        // receive new particles from neighbors
        int rec_cnt = 0;
        particle_t new_particle;
        while(rec_cnt < num_neighbors){
            MPI_Status stat;
            MPI_Recv(&new_particle, 1, PARTICLE,MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            if(new_particle.id == -1){
                rec_cnt++;
                continue;
            }
            addto_local_bins(local_bins, new_particle, local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        // send neighbors 
        for (int offset_x = init_x; offset_x <= end_x; offset_x++){
            for(int offset_y = init_y; offset_y <= end_y; offset_y++){
                int send_to_idx = (proc_y_current + offset_y) * num_proc_x + proc_x_current + offset_x;
                if (offset_x == 0 && offset_y == 0) continue;
                else {
                    std::set<int> neighbor_idx = find_idx(offset_x, offset_y, local_bin_size);
                    for (std::set<int> iterator:: p= neighbor_idx.begin(); p != neighbor_idx.end(); ++p){
                        std::map<double, particle_t> p1_map = local_bins[*p].native_particle;
                        for (std::map<double, particle_t>::iterator p1 = p1_map.begin(); p1 != p1_map.end(); ++p1)
                            MPI_Isend(&end, 1, PARTICLE, send_to_idx, rank, MPI_COMM_WORLD, &request);
                        }
                    } 
                }
            }
        }

        // end 
        for (int offset_x = init_x; offset_x <= end_x; offset_x++){
            for(int offset_y = init_y; offset_y <= end_y; offset_y++){
                int send_to_idx = (proc_y_current + offset_y) * num_proc_x + proc_x_current + offset_x;
                if (offset_x == 0 && offset_y == 0) continue;
                else MPI_Isend(&end, 1, PARTICLE, send_to_idx, rank, MPI_COMM_WORLD, &request);
            }
        }

        rec_cnt = 0;
        while(rec_cnt < num_neighbors){
            MPI_Status stat;
            MPI_Recv(&new_particle, 1, PARTICLE,MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            if(new_particle.id == -1){
                rec_cnt++;
                continue;
            }
            addto_local_bins(local_bins, new_particle, local_bin_size, num_proc_x, num_proc_y, rank, bin_len);
        }

        




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
        fprintf(fsum,"%d %d %g\n",n,number_of_processors,simulation_time);
    }
  

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
