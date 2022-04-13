#include<generator.hpp>

std::vector<Particle> particles; //this list will store all the particles


int main() {

    RandomGenerator generator;

    Domain simDomain (Point3d(0,0,0),Point3d(1,1,1));
    generator.randomParticleGenerator(particles, simDomain);
    std::cout <<"Packing fraction after inserting a lot of particles: "<< generator.calculatePackingFraction(particles,simDomain) << std::endl;
    // deleteParticles(particles);
    generator.addParticlesToGrid(particles, simDomain);
    // printGrid();
    generator.deleteInGrid();
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction after deleting overlap particles: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.randomDeleting(particles, simDomain, 0.6);
    std::cout << "Packing fraction after deleting particles randomly to achieve desired packing fraction: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    // printParticlesList(particles);
    generator.saveToCSV(particles);

    return 0;
}
