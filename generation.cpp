#include<generation.hpp>

int main() {

    Domain simDomain (Point3d(0,0,0),Point3d(1,1,1));
    randomParticleGenerator(particles, simDomain);
    deleteParticles(particles);
    printParticlesList(particles);

    return 0;
}
