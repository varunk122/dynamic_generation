#include<generator.hpp>

std::vector<Particle> particles; //this list will store all the particles


int main(int argc, char *argv[]) {

    if (argc !=9) {
        std::cout << "wrong number of input arguments\n";
        std::cout << "total number of input arguments provided " << argc - 1 << "\n";
        return 0;
    }

    RandomGenerator generator;

    generator.radius_ = std::atof (argv[1]);
    int max_iter = std::atoi (argv[2]);

    int x0 = std::atoi(argv[3]);
    int y0 = std::atoi(argv[4]);
    int z0 = std::atoi(argv[5]);

    int x1 = std::atoi(argv[6]);
    int y1 = std::atoi(argv[7]);
    int z1 = std::atoi(argv[8]);


   
    Domain simDomain (Point3d(x0,y0,z0),Point3d(x1,y1,z1));

    generator.carefulInsertion(max_iter,simDomain);
    generator.addParticlesToList(particles);
    std::cout << "Final Packing Fraction obtained " << generator.calculatePackingFraction(particles, simDomain);
    
    generator.saveToCSV(particles);

    return 0;
}
