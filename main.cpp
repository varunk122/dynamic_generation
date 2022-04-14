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
    generator.newDeleteOverlappingParticles(100);
    generator.newDeleteOverlappingParticles(60);
    generator.newDeleteOverlappingParticles(50);
    generator.newDeleteOverlappingParticles(40);
    generator.newDeleteOverlappingParticles(30);
    generator.newDeleteOverlappingParticles(25);
    generator.newDeleteOverlappingParticles(20);
    generator.newDeleteOverlappingParticles(15);
    generator.newDeleteOverlappingParticles(13);
    generator.newDeleteOverlappingParticles(12);
    generator.newDeleteOverlappingParticles(10);
    generator.newDeleteOverlappingParticles(8);
    generator.newDeleteOverlappingParticles(6);
    generator.newDeleteOverlappingParticles(5);
    generator.newDeleteOverlappingParticles(4);
    generator.newDeleteOverlappingParticles(3);
    generator.newDeleteOverlappingParticles(2);

    std::cout << "Total cells: " << generator.totalCellsInGrid() << std::endl; 
    std::cout << "Total number of empty cells: " << generator.emptyCellsCount() << std::endl; 
    generator.deleteInGrid();
    std::cout << "Total number of empty cells: " << generator.emptyCellsCount() << std::endl; 

    // generator.addParticlesInEmptyCells();
    // std::cout << "Total number of empty cells: " << generator.emptyCellsCount() << std::endl; 

    generator.addParticlesToList(particles);
    std::cout << "Packing fraction after deleting overlap particles: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.randomDeleting(particles, simDomain, 0.6);
    std::cout << "Packing fraction after deleting particles randomly to achieve desired packing fraction: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    // printParticlesList(particles);
    generator.saveToCSV(particles);

    return 0;
}
