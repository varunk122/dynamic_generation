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
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 100 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(60);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 60 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(50);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 50 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(40);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 40 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.newDeleteOverlappingParticles(30);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 30 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(25);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 25 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(20);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 20 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(15);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 15 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(13);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 13 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(12);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 12 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(10);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 10 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(8);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 8 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(6);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 6 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(5);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 5 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);
    generator.newDeleteOverlappingParticles(4);
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction 4 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToGrid(particles, simDomain);

    generator.deleteInGrid(1e-4);


    // generator.newDeleteOverlappingParticles(3);
    // generator.addParticlesToList(particles);
    // std::cout << "Packing fraction 3 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    // generator.addParticlesToGrid(particles, simDomain);
    // generator.newDeleteOverlappingParticles(2);
    // generator.addParticlesToList(particles);
    // std::cout << "Packing fraction 2 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.addParticlesToList(particles);
    std::cout << "Packing fraction after deleting overlap particles: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;

    generator.overlapAnalysis(particles);

    // generator.addParticlesToGrid(particles, simDomain);

    // generator.newDeleteOverlappingParticles(1);
    // generator.addParticlesToList(particles);
    // std::cout << "Packing fraction 1 overlap: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    // generator.addParticlesToGrid(particles, simDomain);

    // std::cout << "Total cells: " << generator.totalCellsInGrid() << std::endl; 
    // std::cout << "Total number of empty cells: " << generator.emptyCellsCount() << std::endl; 
    // generator.deleteInGrid(1e-4);
    // std::cout << "Total number of empty cells: " << generator.emptyCellsCount() << std::endl; 

    // generator.addParticlesInEmptyCells();
    // std::cout << "Total number of empty cells: " << generator.emptyCellsCount() << std::endl; 
    // generator.deleteInGrid();
    // std::cout << "Total number of empty cells: " << generator.emptyCellsCount() << std::endl; 
    // generator.addParticlesToList(particles);
    std::cout << "Packing fraction after deleting overlap particles: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    generator.randomDeleting(particles, simDomain, 0.6);
    std::cout << "Packing fraction after deleting particles randomly to achieve desired packing fraction: " << generator.calculatePackingFraction(particles,simDomain) << std::endl;
    // printParticlesList(particles);
    generator.saveToCSV(particles);

    return 0;
}
