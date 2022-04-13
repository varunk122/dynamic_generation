#include<generator.hpp>

	
RandomGenerator::RandomGenerator() {}

double RandomGenerator::calculateOverlap(Particle& particle1, Particle& particle2) {
	Point3d point1 = particle1.centre;
	Point3d point2 = particle2.centre;
	double dist_ = sqrt(pow((point1(0) - point2(0)),2) + pow((point1(1) - point2(1)),2) + pow((point1(2) - point2(2)),2));
	return dist_ - (particle1.radius + particle2.radius);
}


void RandomGenerator::addParticlesToGrid(std::vector<Particle>& particles, Domain& dom) {
	
	double diameter = 2 * particles[0].radius;
	int nx_ = std::ceil((dom.second(0) - dom.first(0)) / diameter);
	int ny_ = std::ceil((dom.second(1) - dom.first(1)) / diameter);
	int nz_ = std::ceil((dom.second(2) - dom.first(2)) / diameter);

	grid.resize(nx_);
	for (unsigned iX = 0; iX < nx_; iX++) {
		grid[iX].resize(ny_);
		for(unsigned iY = 0; iY < ny_; iY++) {
			grid[iX][iY].resize(nz_);
		}
	}

	while (particles.size() > 0) {
		Point3d coord_ = particles[0].centre - dom.first;
		int idX = std::floor(coord_(0) / diameter); 
		int idY = std::floor(coord_(1) / diameter); 
		int idZ = std::floor(coord_(2) / diameter); 

		grid[idX][idY][idZ].push_back(particles[0]);
		std::swap(particles[0], particles[particles.size() -1]);
		particles.pop_back();
	}

}

void RandomGenerator::printParticle(Particle& particle) {
	std::cout << "id: " << particle.id <<"\n" << "particle centre: " << particle.centre(0) << " " << particle.centre(1)
	<< " " << particle.centre(2) << "\n" << "particle radius " << particle.radius << "\n"; 
} 

void RandomGenerator::printGrid() {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	for (unsigned iX = 0; iX < nx_ ; iX++) {
		for (unsigned iY = 0; iY < ny_; iY++) {
			for (unsigned iZ = 0; iZ < nz_; iZ++) {

				std::cout << iX << " " << iY << " " << iZ << "\n";
				for (unsigned p = 0; p < grid[iX][iY][iZ].size(); p++) {
					printParticle(grid[iX][iY][iZ][p]);
				}
				std::cout << "\n";
			}
		}
	}	

}

void RandomGenerator::deleteOverlappingParticles(Point3d cell1, Point3d cell2, bool same_cell)  {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	if ((cell1(0) >= nx_) || (cell1(1) >= ny_) || (cell1(2) >= nz_)) return;
	if ((cell2(0) >= nx_) || (cell2(1) >= ny_) || (cell2(2) >= nz_)) return;
	if ((cell2(0) < 0) || (cell2(1) < 0) || (cell2(2) < 0)) return;

	std::vector<Particle>& particles_cell1_ = grid[cell1(0)][cell1(1)][cell1(2)];
	std::vector<Particle>& particles_cell2_ = grid[cell2(0)][cell2(1)][cell2(2)];

	if(particles_cell1_.size() == 0) return;
	if(particles_cell2_.size() == 0) return;

	for (unsigned i = 0; i < particles_cell1_.size(); i++) {
		for (unsigned j =0; j < particles_cell2_.size(); j++) {
			if (same_cell && (j<=i)) continue;
			// printParticle(particles_cell2_[j]);

			double overlap_ = calculateOverlap(particles_cell1_[i], particles_cell2_[j]);
			if(overlap_ < 0 && fabs(overlap_) / (2*0.05)  < 1e-4 ) {
				std::swap(particles_cell2_[j], particles_cell2_[particles_cell2_.size()-1]);
				particles_cell2_.pop_back();
				j--;
			}

		}
	}
}

void RandomGenerator::deleteInGrid() {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	for (unsigned iX = 0; iX < nx_ ; iX++) {
		for (unsigned iY = 0; iY < ny_; iY++) {
			for (unsigned iZ = 0; iZ < nz_; iZ++) {

				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX, iY, iZ), true);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX, iY , iZ + 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX, iY - 1, iZ + 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX, iY - 1, iZ), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX, iY - 1, iZ - 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY - 1, iZ - 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY, iZ - 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY + 1, iZ - 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY - 1, iZ), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY, iZ), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY + 1, iZ), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY - 1, iZ + 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY, iZ + 1), false);
				deleteOverlappingParticles(Point3d(iX, iY, iZ), Point3d(iX + 1, iY + 1, iZ + 1), false);

			}
		}
	}	
}

void RandomGenerator::addParticlesToList(std::vector<Particle>& particles) {
	if(particles.size() !=0) return;

	for (unsigned iX = 0; iX < grid.size(); iX++) {
		for	(unsigned iY =0; iY < grid[iX].size(); iY++) {
			for (unsigned iZ =0; iZ < grid[iX][iY].size(); iZ++ ) {
				std::vector<Particle>& p = grid[iX][iY][iZ]; 
				while(p.size() > 0) {
					particles.push_back(p[0]);
					std::swap(p[0], p[p.size()-1]);
					p.pop_back();
				}
			}
		}
	}
}

// std::unordered_map<int,std::unordered_set<std::pair<Particle,Particle>>> distance_list;

double RandomGenerator::calculateSphereVolume(double& r) {
	return (4.0*PI*pow(r,3)) / 3.0;
}

// void generateDistanceList(std::vector<Particle>& particles) {
// 	for(int i =0; i< particles.size(); i++) {
// 		for(int j=0; j< particles.size();j++) {
// 			double overlap_ = calculateOverlap(particle)
// 		}
// 	}
// }

double RandomGenerator::calculatePackingFraction(std::vector<Particle>& particles, Domain & dom) {
	double spheres_vol_ = 0;
	for(int i =0; i< particles.size(); i++) {
		spheres_vol_ += calculateSphereVolume(particles[i].radius);
	} 

	double simulation_vol_ = (dom.second(0)- dom.first(0))*(dom.second(1)- dom.first(1))*(dom.second(2)- dom.first(2));

	return spheres_vol_ / simulation_vol_;

}

void RandomGenerator::updatePackingFraction(double& packing_fraction, Domain& dom, Particle& deleted_particle) {

	double simulation_vol_ = (dom.second(0)- dom.first(0))*(dom.second(1)- dom.first(1))*(dom.second(2)- dom.first(2));
	packing_fraction = (packing_fraction * simulation_vol_ - calculateSphereVolume(deleted_particle.radius)) / simulation_vol_;

}

void RandomGenerator::randomDeleting(std::vector<Particle>& particles, Domain& dom, double required_pf) {
	int total = particles.size();
	double current_pf = calculatePackingFraction(particles, dom);
	while(current_pf - required_pf > 0 && total) {
		int random_particle = std::rand() % total;
		updatePackingFraction(current_pf, dom, particles[random_particle]);
		// particles.erase(particles.begin() + random_particle);
		std::swap(particles[random_particle], particles.back());
		particles.pop_back();
		total--;
	}
}


void RandomGenerator::printParticlesList(std::vector<Particle>& particles) {
	
	for (int i =0; i< particles.size();i++) {
		std::cout << particles[i].centre(0) << " " << particles[i].centre(1) << " " << particles[i].centre(2) << std::endl;
	}
}

void RandomGenerator::saveToCSV(std::vector<Particle>& particles) {

	std::fstream particle_file_;
	particle_file_.open("./values.csv", std::ios::out | std::ios::app);

	for (int i =0; i< particles.size();i++) {
		particle_file_ << particles[i].centre(0) << ", " << particles[i].centre(1) << ", " << particles[i].centre(2) << ", "
						<< particles[i].radius << "\n";
	}

}

void RandomGenerator::deleteParticles(std::vector<Particle>& particles) {
	std::vector<unsigned> toDelete_;
	for (int i =0; i < particles.size(); i++) {
		for (int j =i+1; j <particles.size(); j++) {
			double overlap_ = calculateOverlap(particles[i], particles[j]);
			if( overlap_ <= 0 && fabs(overlap_) / (2*0.05)  < 1e-4 ) {
				// particles.erase(particles.begin() + j);
				// j--;
				std::swap(particles[j], particles.back());
				particles.pop_back();
				j--;
			}
		}
	}
}

void RandomGenerator::randomParticleGenerator(std::vector<Particle>& particles, Domain& dom) {
	int id_ = 0;
	double raidus_ = 0.05;
	float packingFraction_ = 0;
	double simulationVolume_ = (dom.second(0)- dom.first(0))*(dom.second(1)- dom.first(1))*(dom.second(2)- dom.first(2));
	double spheresVolume_ =0;

	// uniform real distribution

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution_x (dom.first(0), dom.second(0));
	std::uniform_real_distribution<double> distribution_y (dom.first(1), dom.second(1));
	std::uniform_real_distribution<double> distribution_z (dom.first(2), dom.second(2));

	while (packingFraction_ < 2) {
		particles.push_back(Particle(id_++, Point3d(distribution_x(generator), distribution_y(generator), distribution_z(generator)),raidus_));
		spheresVolume_ += calculateSphereVolume(raidus_);
		packingFraction_ = spheresVolume_/simulationVolume_;
	}

	std::cout << "Total Particles added into simulation domain: " << particles.size() << std::endl;
}
