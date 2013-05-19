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


int main(int argc, char const *argv[])
{
    int taskid, ntasks;

#ifdef _MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#endif

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


    for (int t = 0; t < MAX_TIME; ++t) {
        // For each particle
        for (int i = 0; i < particles.size(); ++i) {
            // Check if the particle collides with any other particle
            for (int j = i + 1; (j < particles.size()); ++j) {
                float collission = collide(&particles[i].pcord,
                                           &particles[j].pcord);

                if (collission >= 0) {
                    interact(&particles[i].pcord,
                             &particles[j].pcord, collission);

                    // Collission with particle and then collission with wall.
                    // This should be very unlikely to happen
                    total_momentum += wall_collide(&particles[j].pcord, wall);
                } else {
                    feuler(&particles[i].pcord, 1);
                }
            }
            // Check if particle collides with a wall and and that to total momentum
            total_momentum += wall_collide(&particles[i].pcord, wall);
        }
    }
    printf("Total pressure is %.2f\n", total_momentum / (MAX_TIME * WALL_LENGTH));

#ifdef _MPI
    MPI_Finalize();
#endif

    return 0;
}
