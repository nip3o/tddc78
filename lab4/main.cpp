#include <vector>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

// #define _MPI

#ifdef _MPI
#include <mpi.h>
#endif

#include "physics.h"
#include "definitions.h"

#define ROOT 0

void print_pcoord(pcord_t c) {
    printf("x %.2f, y %.2f, vx %.2f, vy %.2f\n", c.x, c.y, c.vx, c.vy);
}

int destination(int area_size, pcord_t c) {
    return (int)c.x % (int)area_size;
}

int main(int argc, char const *argv[]) {
    int taskid = 0, ntasks = 4;

#ifdef _MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#endif

    const int AREA_SIZE = BOX_HORIZ_SIZE / ntasks;
    std::vector<particle_t>* travellers = new std::vector<particle_t>[ntasks];

    double total_momentum = 0.0;

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

        c.x = std::rand() % (int)BOX_HORIZ_SIZE;
        c.y = std::rand() % (int)BOX_VERT_SIZE;
        c.vx = r * cos(angle);
        c.vy = r * sin(angle);

//        print_pcoord(c);

        p.pcord = c;
        p.ptype = 0;

        particles.push_back(p);
    }

    printf("Created %d particles\n", (int)particles.size());
    int new_destination;
    float collission;

    for (int t = 0; t < MAX_TIME; ++t) {
        // For each particle
        for (int i = 0; i < particles.size(); ++i) {
            collission = -1;

            // Check if the particle collides with any other particle
            for (int j = i + 1; (j < particles.size()); ++j) {
                collission = collide(&particles[i].pcord,
                                           &particles[j].pcord);

                if (collission >= 0) {
                    interact(&particles[i].pcord,
                             &particles[j].pcord, collission);

                    // Collission with particle and then collission with wall
                    total_momentum += wall_collide(&particles[j].pcord, wall);
                    continue;
                }
            }

            if (collission < 0) {
                // The particle did not collide, move it
                feuler(&particles[i].pcord, 1);

                // Check if particle collides with a wall and add that to total momentum
                total_momentum += wall_collide(&particles[i].pcord, wall);
            }

            new_destination = destination(AREA_SIZE, particles[i].pcord);
            if(new_destination != taskid) {
                // Mark a new destination for the particle
                travellers[new_destination].push_back(particles[i]);
                // Add it to the list of particles to be deleted from this node
                travellers[taskid].push_back(particles[i]);
            }
        }
    }
    printf("Total pressure is %.2f\n", total_momentum / (MAX_TIME * WALL_LENGTH));

#ifdef _MPI
    MPI_Finalize();
#endif

    return 0;
}
