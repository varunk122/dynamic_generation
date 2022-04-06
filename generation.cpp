#include<generation.hpp>

int main() {

    Domain simDomain (Point3d(0,0,0),Point3d(1,1,1));
    randomParticleGenerator(particles, simDomain);
    std::cout <<"Packing fraction after inserting a lot of particles: "<< calculatePackingFraction(particles,simDomain) << std::endl;
    deleteParticles(particles);
    std::cout << "Packing fraction after deleting overlap particles: " << calculatePackingFraction(particles,simDomain) << std::endl;
    randomDeleting(particles, simDomain, 0.6);
    std::cout << "Packing fraction after deleting particles randomly to achieve desired packing fraction: " << calculatePackingFraction(particles,simDomain) << std::endl;
    // printParticlesList(particles);
    saveToCSV(particles);

    return 0;
}
