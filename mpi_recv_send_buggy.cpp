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
    int num_proc_y = (int) floor(n_proc / num_proc_x); // The number of processors along the y-axis. 
    // printf("Number of total processes %d", n_proc);
    // fflush(stdout);
    // printf("Number of processes_x %d\n", num_proc_x);
    // fflush(stdout);
    // printf("Number of processes_y %d\n", num_proc_y);
    // fflush(stdout);


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

        std::cout<<"Reached line 89 in rank 0"<<std::endl;
        for (int i = 0; i < n; i++) // for each particle
        {
            int *process_ids = (int *) malloc(9 * sizeof(int));
            process_ids = get_procs(particles[i].x, particles[i].y, num_proc_x, num_proc_y);

            for (int j = 0; j < 9; j++)
            {
                if (process_ids[j] != -1)
                {
                    partition_sizes[process_ids[j]] += 1;
                    // std::cout<<"\n\nParticle "<<i<<" goes to processor "<<" "<<process_ids[j];
                }

            }
            free(process_ids);
        }
        // std::cout<<"Reached line 103 in rank 0";
    }



    // Debugging
    if (rank == 0)
    {
        for (int i = 0; i < n_proc; i++)
        {
            std::cout<<"Partition size "<<i<<" is equal to "<<partition_sizes[i]<<std::endl;
        }
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    // Send an array of sizes (array of ints) to each processor first. 
    int *local_size = (int *)malloc(sizeof(int)); // This is where we will recieve the size. 
    std::cout<<"Reached line 104 in rank "<<rank<<std::endl;
    MPI_Scatter(partition_sizes, 1, MPI_INT, local_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::cout<<"Reached line 111 in rank "<<rank<<std::endl;

    // Debugging
    // if (rank == 0)
    // {
    //     std::cout<<"I am processor 0 \n";
    //     std::cout<<(*local_size)<<std::endl;
    // }
    // else if (rank == 1)
    // {
    //     std::cout<<"I am processor 1 \n";
    //     std::cout<<(*local_size)<<std::endl;
    // } 



    // Each processor allocates space.
    particle_t *local_particles = (particle_t*) malloc( *local_size * sizeof(particle_t) );
    std::cout<<"Rank "<<rank<<" allocated space for local particles"<<std::endl;

    int *partition_offsets = (int*) malloc( n_proc * sizeof(int) );
    int *offsets_copy = (int*) malloc (n_proc * sizeof(int));



    particle_t *particles_to_scatter;


    if (rank == 0)
    {
        std::cout<<"Reached line 149 in rank 0 "<<std::endl;
        partition_offsets[0] = 0;

        for (int i = 1; i < n_proc; i ++)
        {
            partition_offsets[i] = partition_offsets[i-1] + partition_sizes[i-1];
        }
    

        particles_to_scatter = (particle_t*) malloc ((partition_offsets[n_proc-1] + partition_sizes[n_proc-1]) * sizeof(particle_t));

        // Recieve the particles for this processor into local_partices
        
        for (int i = 0; i < n_proc; i++)
        {
            offsets_copy[i] = partition_offsets[i];
        }
        std::cout<<"Reached line 166 in rank 0"<<std::endl;
        
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
                    // std::cout<<"\n\nParticle "<<i<<" goes to processor "<<" "<<process_ids[j];

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
        std::cout<<"Reached line 180 in rank 0 "<<std::endl;
    }



    MPI_Barrier(MPI_COMM_WORLD); //TODO: Possibly remove this.
    MPI_Scatterv( particles_to_scatter, partition_sizes, partition_offsets, PARTICLE,
             local_particles, *local_size, PARTICLE, 0, MPI_COMM_WORLD );



    /*
        ******************************************************************
        code to send/receive particles between processors
        ******************************************************************
    */

    // manually change the location of some particles to other processors for experiment
    if (rank == 0){
        for (int i = 0; i < 10; i++){
            local_particles[i].x = get_size() - 0.00001;
        }
    }
    if (rank == 0){
        for (int i = 0; i < 10; i++){
            local_particles[i].x = 0 + 0.00001 ;
        }
    }
    // assuming that I know the native particles, the number of native particles, and their new x and y in each processor
    int proc_x_new, proc_y_new;
    int proc_y_current = rank/num_proc_x;
    int proc_x_current = rank - proc_y_current * num_proc_x;

    // size of array of particles to be sent to 8 neighboring processors
    int *send_size_up = (int *)malloc(sizeof(int));
    *send_size_up = 0;
    int *send_size_upperleft = (int *)malloc(sizeof(int));
    *send_size_upperleft = 0;
    int *send_size_left = (int *)malloc(sizeof(int));
    *send_size_left = 0;
    int *send_size_lowerleft = (int *)malloc(sizeof(int));
    *send_size_lowerleft = 0;
    int *send_size_down = (int *)malloc(sizeof(int));
    *send_size_down = 0;
    int *send_size_lowerright = (int *)malloc(sizeof(int));
    *send_size_lowerright = 0;
    int *send_size_right = (int *)malloc(sizeof(int));
    *send_size_right = 0;
    int *send_size_upperright = (int *)malloc(sizeof(int));
    *send_size_upperright = 0;
    // store the indices of particles to be sent, save time
    std::set<int>  index_send; 

    /*
      calculate the send_size of array of particles to be sent to 8 neighboring processors
    */
    for (int i = 0; i < *local_size; i++){
            proc_x_new =  get_proc_x(local_particles[i].x, num_proc_x);
            proc_y_new =  get_proc_y(local_particles[i].y, num_proc_y);
            if ((proc_y_new * num_proc_x) + proc_x_new != rank){ // if the native particles moves to another processor
                index_send.insert(i);
                // up
                if(proc_x_new == proc_x_current && proc_y_new == proc_y_current + 1){
                    *send_size_up += 1;
                }
                // upper left
                if(proc_x_new == proc_x_current - 1 && proc_y_new == proc_y_current + 1){
                    *send_size_upperleft += 1;
                }
                // left
                if(proc_x_new == proc_x_current - 1 && proc_y_new == proc_y_current){
                    *send_size_left += 1;
                }
                // lower left
                if(proc_x_new == proc_x_current - 1 && proc_y_new == proc_y_current - 1){
                    *send_size_lowerleft += 1;
                }
                // down
                if(proc_x_new == proc_x_current && proc_y_new == proc_y_current - 1){
                    *send_size_down += 1;
                }
                // lower right
                if(proc_x_new == proc_x_current + 1 && proc_y_new == proc_y_current - 1){
                    *send_size_lowerright += 1;
                }
                // right
                if(proc_x_new == proc_x_current + 1 && proc_y_new == proc_y_current){
                    *send_size_right += 1;
                }
                // upper right
                if(proc_x_new == proc_x_current + 1 && proc_y_new == proc_y_current + 1){
                    *send_size_upperright += 1;
                }

            }
        }
    /*
      assign memory for 8 arrays of particles to be sent
    */
    particle_t *particles_send_up = (particle_t*) malloc (*send_size_up * sizeof(particle_t));
    particle_t *particles_send_upperleft = (particle_t*) malloc (*send_size_upperleft * sizeof(particle_t));
    particle_t *particles_send_left = (particle_t*) malloc (*send_size_left * sizeof(particle_t));
    particle_t *particles_send_lowerleft = (particle_t*) malloc (*send_size_lowerleft * sizeof(particle_t));
    particle_t *particles_send_down = (particle_t*) malloc (*send_size_down * sizeof(particle_t));
    particle_t *particles_send_lowerright = (particle_t*) malloc (*send_size_lowerright * sizeof(particle_t));
    particle_t *particles_send_right = (particle_t*) malloc (*send_size_right * sizeof(particle_t));
    particle_t *particles_send_upperright = (particle_t*) malloc (*send_size_upperright * sizeof(particle_t));
   
    /*
      populate these 8 arrays of particles to be sent
    */
    std::set<int>::iterator it2;
    int index_up = 0;
    int index_upperleft = 0;
    int index_left = 0;
    int index_lowerleft = 0;
    int index_down = 0;
    int index_lowerright = 0;
    int index_right = 0;
    int index_upperright = 0;

    for (it2 = index_send.begin(); it2 != index_send.end(); ++it2){
        std::cout<<"test of it2 "<<*it2<<std::endl;

        proc_x_new =  get_proc_x(local_particles[*it2].x, num_proc_x);
        proc_y_new =  get_proc_y(local_particles[*it2].y, num_proc_y);
        // up
        if(proc_x_new == proc_x_current && proc_y_new == proc_y_current + 1){
            particles_send_up[index_up++] = local_particles[*it2];
        }
        // upper left
        if(proc_x_new == proc_x_current - 1 && proc_y_new == proc_y_current + 1){
            particles_send_upperleft[index_upperleft++] = local_particles[*it2];
        }
        // left
        if(proc_x_new == proc_x_current - 1 && proc_y_new == proc_y_current){
            particles_send_left[index_left++] = local_particles[*it2];
        }
        // lower left
        if(proc_x_new == proc_x_current - 1 && proc_y_new == proc_y_current - 1){
            particles_send_lowerleft[index_lowerleft++] = local_particles[*it2];
        }
        // down
        if(proc_x_new == proc_x_current && proc_y_new == proc_y_current - 1){
            particles_send_down[index_down++] = local_particles[*it2];
        }
        // lower right
        if(proc_x_new == proc_x_current + 1 && proc_y_new == proc_y_current - 1){
            particles_send_lowerright[index_lowerright++] = local_particles[*it2];
        }
        // right
        if(proc_x_new == proc_x_current + 1 && proc_y_new == proc_y_current){
            particles_send_right[index_right++] = local_particles[*it2];
        }
        // upper right
        if(proc_x_new == proc_x_current + 1 && proc_y_new == proc_y_current + 1){
            particles_send_upperright[index_upperright++] = local_particles[*it2];
        }
    }
    
    /*
      first send the 8 integer of number of particles, MPI can send empty messages, so always send
    */
    // up
    if(rank + num_proc_x >=0 && rank + num_proc_x < num_proc_x*num_proc_y){
        MPI_Send(send_size_up, 1, MPI_INT, rank + num_proc_x, 0, MPI_COMM_WORLD);
    }
    // upperleft
    if(rank + num_proc_x - 1 >=0 && rank + num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Send(send_size_upperleft, 1, MPI_INT, rank + num_proc_x - 1, 0, MPI_COMM_WORLD);
    }
    // left
    if(rank - 1 >=0 && rank - 1 < num_proc_x*num_proc_y){
        MPI_Send(send_size_left, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
    }
    // lowerleft
    if(rank - num_proc_x - 1 >=0 && rank - num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Send(send_size_lowerleft, 1, MPI_INT, rank - num_proc_x - 1, 0, MPI_COMM_WORLD);
    }
    // down
    if(rank - num_proc_x >=0 && rank - num_proc_x < num_proc_x*num_proc_y){
        MPI_Send(send_size_down, 1, MPI_INT, rank - num_proc_x, 0, MPI_COMM_WORLD);
    }
    // lowerright
    if(rank - num_proc_x + 1 >=0 && rank - num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Send(send_size_lowerright, 1, MPI_INT, rank - num_proc_x + 1, 0, MPI_COMM_WORLD);
    }
    // right
    if(rank + 1 >=0 && rank + 1 < num_proc_x*num_proc_y){
        MPI_Send(send_size_right, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
    // upperright
    if(rank + num_proc_x + 1 >=0 && rank + num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Send(send_size_upperright, 1, MPI_INT, rank + num_proc_x + 1, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD); //

    /*
      first receive the 8 integer of number of particles
    */
    int *receive_size_up = (int *)malloc(sizeof(int));
    *receive_size_up = 0;
    int *receive_size_upperleft = (int *)malloc(sizeof(int));
    *receive_size_upperleft = 0;
    int *receive_size_left = (int *)malloc(sizeof(int));
    *receive_size_left = 0;
    int *receive_size_lowerleft = (int *)malloc(sizeof(int));
    *receive_size_lowerleft = 0;
    int *receive_size_down = (int *)malloc(sizeof(int));
    *receive_size_down = 0;
    int *receive_size_lowerright = (int *)malloc(sizeof(int));
    *receive_size_lowerright = 0;
    int *receive_size_right = (int *)malloc(sizeof(int));
    *receive_size_right = 0;
    int *receive_size_upperright = (int *)malloc(sizeof(int));
    *receive_size_upperright = 0;
    // up
    if(rank + num_proc_x >=0 && rank + num_proc_x < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_up, 1, MPI_INT, rank + num_proc_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // upperleft
    if(rank + num_proc_x - 1 >=0 && rank + num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_upperleft, 1, MPI_INT, rank + num_proc_x - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // left
    if(rank - 1 >=0 && rank - 1 < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_left, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // lowerleft
    if(rank - num_proc_x - 1 >=0 && rank - num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_lowerleft, 1, MPI_INT, rank - num_proc_x - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // down
    if(rank - num_proc_x >=0 && rank - num_proc_x < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_down, 1, MPI_INT, rank - num_proc_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // lowerright
    if(rank - num_proc_x + 1 >=0 && rank - num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_lowerright, 1, MPI_INT, rank - num_proc_x + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // right
    if(rank + 1 >=0 && rank + 1 < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_right, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // upperright
    if(rank + num_proc_x + 1 >=0 && rank + num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Recv(receive_size_upperright, 1, MPI_INT, rank + num_proc_x + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    /*
      then send the 8 arrays of particles, MPI can send empty messages, so always send
    */  
    if(rank + num_proc_x >=0 && rank + num_proc_x < num_proc_x*num_proc_y){
        MPI_Send(particles_send_up, *send_size_up, PARTICLE, rank + num_proc_x, 0, MPI_COMM_WORLD);
    }
    // upperleft
    if(rank + num_proc_x - 1 >=0 && rank + num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Send(particles_send_upperleft, *send_size_upperleft, PARTICLE, rank + num_proc_x - 1, 0, MPI_COMM_WORLD);
    }
    // left
    if(rank - 1 >=0 && rank - 1 < num_proc_x*num_proc_y){
        MPI_Send(particles_send_left, *send_size_left, PARTICLE, rank - 1, 0, MPI_COMM_WORLD);
    }
    // lowerleft
    if(rank - num_proc_x - 1 >=0 && rank - num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Send(particles_send_lowerleft, *send_size_lowerleft, PARTICLE, rank - num_proc_x - 1, 0, MPI_COMM_WORLD);
    }
    // down
    if(rank - num_proc_x >=0 && rank - num_proc_x < num_proc_x*num_proc_y){
        MPI_Send(particles_send_down, *send_size_down, PARTICLE, rank - num_proc_x, 0, MPI_COMM_WORLD);
    }
    // lowerright
    if(rank - num_proc_x + 1 >=0 && rank - num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Send(particles_send_lowerright, *send_size_lowerright, PARTICLE, rank - num_proc_x + 1, 0, MPI_COMM_WORLD);
    }
    // right
    if(rank + 1 >=0 && rank + 1 < num_proc_x*num_proc_y){
        MPI_Send(particles_send_right, *send_size_right, PARTICLE, rank + 1, 0, MPI_COMM_WORLD);
    }
    // upperright
    if(rank + num_proc_x + 1 >=0 && rank + num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Send(particles_send_upperright, *send_size_upperright, PARTICLE, rank + num_proc_x + 1, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD); //

    /*
      assign memory for 8 arrays of particles to be received
    */ 
    particle_t *particles_receive_up = (particle_t*) malloc (*receive_size_up * sizeof(particle_t));
    particle_t *particles_receive_upperleft = (particle_t*) malloc (*receive_size_upperleft * sizeof(particle_t));
    particle_t *particles_receive_left = (particle_t*) malloc (*receive_size_left * sizeof(particle_t));
    particle_t *particles_receive_lowerleft = (particle_t*) malloc (*receive_size_lowerleft * sizeof(particle_t));
    particle_t *particles_receive_down = (particle_t*) malloc (*receive_size_down * sizeof(particle_t));
    particle_t *particles_receive_lowerright = (particle_t*) malloc (*receive_size_lowerright * sizeof(particle_t));
    particle_t *particles_receive_right = (particle_t*) malloc (*receive_size_right * sizeof(particle_t));
    particle_t *particles_receive_upperright = (particle_t*) malloc (*receive_size_upperright * sizeof(particle_t));
    
    /*
      receive the 8 arrays of particles
    */
    // up
    if(rank + num_proc_x >=0 && rank + num_proc_x < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_up, *receive_size_up, PARTICLE, rank + num_proc_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // upperleft
    if(rank + num_proc_x - 1 >=0 && rank + num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_upperleft, *receive_size_upperleft, PARTICLE, rank + num_proc_x - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // left
    if(rank - 1 >=0 && rank - 1 < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_left, *receive_size_left, PARTICLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // lowerleft
    if(rank - num_proc_x - 1 >=0 && rank - num_proc_x - 1 < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_lowerleft, *receive_size_lowerleft, PARTICLE, rank - num_proc_x - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // down
    if(rank - num_proc_x >=0 && rank - num_proc_x < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_down, *receive_size_down, PARTICLE, rank - num_proc_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // lowerright
    if(rank - num_proc_x + 1 >=0 && rank - num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_lowerright, *receive_size_lowerright, PARTICLE, rank - num_proc_x + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // right
    if(rank + 1 >=0 && rank + 1 < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_right, *receive_size_right, PARTICLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // upperright
    if(rank + num_proc_x + 1 >=0 && rank + num_proc_x + 1 < num_proc_x*num_proc_y){
        MPI_Recv(particles_receive_upperright, *receive_size_upperright, PARTICLE, rank + num_proc_x + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD); //

    /*
        finally do some test
    */
    // should receive 10 points from the right
    if(rank == 0){
        std::cout<<"rank 0**************************************************"<<std::endl;
        std::cout<<"number of particles received from up"<<receive_size_up<<std::endl;
        std::cout<<"number of particles received from upperleft"<<receive_size_upperleft<<std::endl;
        std::cout<<"number of particles received from left"<<receive_size_left<<std::endl;
        std::cout<<"number of particles received from lowerleft"<<receive_size_lowerleft<<std::endl;
        std::cout<<"number of particles received from down"<<receive_size_down<<std::endl;
        std::cout<<"number of particles received from lowerright"<<receive_size_lowerright<<std::endl;
        std::cout<<"number of particles received from right"<<receive_size_right<<std::endl;
        std::cout<<"number of particles received from upperright"<<receive_size_upperright<<std::endl;
    }
    // should receive 10 points from the left
    if(rank == 1){
        std::cout<<"rank 1**************************************************"<<std::endl;
        std::cout<<"number of particles received from up"<<receive_size_up<<std::endl;
        std::cout<<"number of particles received from upperleft"<<receive_size_upperleft<<std::endl;
        std::cout<<"number of particles received from left"<<receive_size_left<<std::endl;
        std::cout<<"number of particles received from lowerleft"<<receive_size_lowerleft<<std::endl;
        std::cout<<"number of particles received from down"<<receive_size_down<<std::endl;
        std::cout<<"number of particles received from lowerright"<<receive_size_lowerright<<std::endl;
        std::cout<<"number of particles received from right"<<receive_size_right<<std::endl;
        std::cout<<"number of particles received from upperright"<<receive_size_upperright<<std::endl;
    }

    /*
        ******************************************************************
        code to send/receive particles between processors
        ******************************************************************
    */

    // Debugging
    // if (rank == 0)
    // {
    //     for (int i = 0; i < *local_size; i++)
    //     {
    //         std::cout<<"\n\n i, and the x pos is "<<i<<" "<<local_particles[i].x;
    //     }
    // }
    
    // int bin_len = bin_length(n);
    // int *local_bin_size = get_bin_size(num_proc_x, num_proc_y, rank, bin_len);
    // int local_bin_row = local_bin_size[0];
    // int local_bin_col = local_bin_size[1];
    // bin_t *local_bins = new bin_t[local_bin_row*local_bin_col];

    // // each bins include particles on left and up edges 
    // // and the right most particles belongs to the right most bins
    // init_local_bins(local_bins, local_particles, *local_size, local_bin_size ,rank, bin_len);

    // int bin_per_proc;
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
