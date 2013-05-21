#include <mpi.h>

#include <time.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include "physics.h"
#include "definitions.h"

#define ROOT 0

int destination(const int area_size, pcord_t c) {
    return (int)c.x / area_size;
}

int main(int argc, char *argv[]) {
    int taskid = 0, ntasks = 2;
    struct timespec stime, etime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    MPI_Status status;
    MPI_Datatype pcoord_mpi_type;

    MPI_Type_contiguous( 4, MPI_FLOAT, &pcoord_mpi_type );
    MPI_Type_commit(&pcoord_mpi_type);

    const int AREA_SIZE = BOX_HORIZ_SIZE / ntasks;
    std::vector<pcord_t>* travellers = new std::vector<pcord_t>[ntasks];

    double total_momentum, local_total_momentum = 0.0;

    std::vector<particle_t> particles;
    cord_t wall;
    wall.x0 = 0;
    wall.y0 = 0;
    wall.x1 = BOX_HORIZ_SIZE;
    wall.y1 = BOX_VERT_SIZE;

    srand(std::time(0));

    for (int i = 0; i < INIT_NO_PARTICLES; ++i) {
        pcord_t c;
        particle_t p;

        float r = std::rand() % MAX_INITIAL_VELOCITY;
        float angle = std::rand() * M_PI * 2;

        const float LOW = AREA_SIZE * taskid;
        const float HIGH = AREA_SIZE * (taskid + 1);
        c.x = LOW + (float)std::rand()/((float)RAND_MAX/(HIGH - LOW));
        c.y = std::rand() % (int)BOX_VERT_SIZE;
        c.vx = r * cos(angle);
        c.vy = r * sin(angle);

        p.pcord = c;
        p.ptype = 0;

        particles.push_back(p);
    }

    int new_destination;
    float collission;
    pcord_t* recieve_buffer;
    recieve_buffer = new pcord_t[PARTICLE_BUFFER_SIZE];

    clock_gettime(CLOCK_REALTIME, &stime);


    for (int t = 0; t < MAX_TIME; ++t) {
        int i = 0;
        int particles_size = particles.size();
        // For each particle
        while (i < particles_size) {
            collission = -1;

            // Check if the particle collides with any other particle
            for (int j = i + 1; (j < particles.size()); ++j) {
                collission = collide(&particles[i].pcord,
                                     &particles[j].pcord);

                if (collission >= 0) {
                    interact(&particles[i].pcord,
                             &particles[j].pcord, collission);
                    // Collission with particle and then collission with wall
                    local_total_momentum += wall_collide(&particles[j].pcord, wall);
                    continue;
                }
            }

            if (collission < 0) {
                // The particle did not collide, move it
                feuler(&particles[i].pcord, 1);

                // Check if particle collides with a wall and add that to total momentum
                local_total_momentum += wall_collide(&particles[i].pcord, wall);
            }

            new_destination = destination(AREA_SIZE, particles[i].pcord);
            if(new_destination != taskid) {
                // Mark a new destination for the particle
                travellers[new_destination].push_back(particles[i].pcord);

                particles.erase(particles.begin() + i);
                particles_size = particles.size();
            } else {
                i++;
            }
        } // end particle-loop


        pcord_t cord;
        particle_t particle;
        MPI_Request req;
        int sent_particles = 0;

        // Send data asynchronously to all the other nodes
        for (int i = 0; i < ntasks; ++i) {
            if (i != taskid) {
                // Send data to node i
                MPI_Isend(&travellers[i][0], travellers[i].size(), pcoord_mpi_type, i, 0, MPI_COMM_WORLD, &req);
                sent_particles += travellers[i].size();
            }
        }

        printf("Task %d sent %d particles, t=%d\n", taskid, sent_particles, t);

        // Recieve data from all the other nodes
        for (int j = 0; j < ntasks - 1; ++j) {
            int recieved_length, flag = 0;

            while(!flag) {
                MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
            }
            MPI_Get_count(&status, pcoord_mpi_type, &recieved_length);
            MPI_Recv(recieve_buffer, recieved_length, pcoord_mpi_type, status.MPI_SOURCE, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);

            for (int k = 0; k < recieved_length; ++k) {
                pcord_t coord;
                coord.x = recieve_buffer[k].x;
                coord.y = recieve_buffer[k].y;
                coord.vx = recieve_buffer[k].vx;
                coord.vy = recieve_buffer[k].vy;

                particle.pcord = coord;
                particle.ptype = 0;
                particles.push_back(particle);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < ntasks; ++i) {
            travellers[i].clear();
        }
    } // end timestep-loop

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&local_total_momentum, &total_momentum, ntasks, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

    if (taskid == ROOT) {
        printf("Total pressure is %e\n", total_momentum / (MAX_TIME * WALL_LENGTH));

        clock_gettime(CLOCK_REALTIME, &etime);
        printf("Everything took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
               1e-9*(etime.tv_nsec  - stime.tv_nsec));
    }

    MPI_Finalize();

    return 0;
}
