#include <vector>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include "physics.h"
#include "definitions.h"


int main(int argc, char const *argv[])
{
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
        float angle = std::rand() * 2;

        printf("angle %.2f\n", angle);

        c.x = std::rand() % (int)BOX_HORIZ_SIZE;
        c.y = std::rand() % (int)BOX_VERT_SIZE;
        c.vx = r * cos(angle);
        c.vy = r * sin(angle);

        printf("x %.2f, y %.2f, vx %.2f, vy %.2f\n", c.x, c.y, c.vx, c.vy);

        p.pcord = c;
        p.ptype = 0;

        particles.push_back(p);
    }

    printf("Created %d particles\n", particles.size());


    for (int t = 0; t < MAX_TIME; ++t) {

        for (int i = 0; i < particles.size(); ++i) {
            for (int j = i + 1; (j < particles.size()); ++j) {
                float collission = collide(&particles[i].pcord,
                                           &particles[j].pcord);

                if (collission >= 0) {
                    printf("Collission! \n");
                    interact(&particles[i].pcord,
                             &particles[j].pcord, collission);
                } else {
                    feuler(&particles[i].pcord, t);
                }

                wall_collide(&particles[i].pcord, wall);
                wall_collide(&particles[j].pcord, wall);
            }
        }
    }

    return 0;
}
