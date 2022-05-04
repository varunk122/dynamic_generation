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

		if(idX < nx_ && idY < ny_ && idZ < nz_) {
			grid[idX][idY][idZ].push_back(particles[0]);
		}
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


unsigned RandomGenerator::totalCellsInGrid() {
	unsigned nx_ = grid.size();
	unsigned ny_ = grid[0].size();
	unsigned nz_ = grid[0][0].size();

	return nx_ * ny_ * nz_;

}

unsigned RandomGenerator::emptyCellsCount() {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	unsigned count = 0;

	for (unsigned iX = 0; iX < nx_ ; iX++) {
		for (unsigned iY = 0; iY < ny_; iY++) {
			for (unsigned iZ = 0; iZ < nz_; iZ++) {

				if(grid[iX][iY][iZ].size() == 0) {
					count += 1;
				}	
			}
		}
	}	

	return count;
}

void RandomGenerator::addParticlesInEmptyCells() {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	// double radius_ = 0.02;

	for (unsigned iX = 0; iX < nx_ ; iX++) {
		for (unsigned iY = 0; iY < ny_; iY++) {
			for (unsigned iZ = 0; iZ < nz_; iZ++) {
				if(grid[iX][iY][iZ].size() == 0) {
					Point3d centre_(iX, iY, iZ);
					centre_ = centre_ * 2 * radius_ + Point3d(radius_,radius_,radius_);
					grid[iX][iY][iZ].push_back(Particle(ID++, centre_, radius_));
				}
			}
		}
	}

}

int RandomGenerator::checkIfOverlappingWithParticlesInCell(Point3i cell, Particle particle) {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	int overlaps_ = 0;

	if ((cell(0) >= nx_) || (cell(1) >= ny_) || (cell(2) >= nz_)) return 0;
	if ((cell(0) < 0) || (cell(1) < 0) || (cell(2) < 0)) return 0;

	std::vector<Particle>& cell_particles = grid[cell(0)][cell(1)][cell(2)];

	for (unsigned i = 0; i < cell_particles.size(); i++) {
			if(cell_particles[i].id == particle.id) continue;
			double overlap_ = calculateOverlap(cell_particles[i], particle);
			if(overlap_ < 0 && fabs(overlap_) / (2*radius_)  > 1e-4 ) {
				overlaps_++;
			}
	}

	return overlaps_;

}

bool RandomGenerator::checkIfOverlapping(Point3i cell, Particle particle, int overlap_count) {

	int x = cell(0);
	int y = cell(1);
	int z = cell(2);

	int total_overlaps = 0;

	for (int8_t dx = -1; dx <= 1; dx++) {
		for (int8_t dy = -1; dy <= 1; dy++) {
			for (int8_t dz = -1; dz <= 1; dz++) {
				int x_ = x + dx;
				int y_ = y + dy;
				int z_ = z + dz;
				// if(dx == 0 && dy ==0 && dz == 0) continue;
				total_overlaps += checkIfOverlappingWithParticlesInCell(Point3i(x_,y_,z_), particle);
					if (total_overlaps >= overlap_count) return true;
				
			}
		}
	}

	return false;

}

void  RandomGenerator::newDeleteOverlappingParticles(int overlap_count) {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	for (unsigned iX = 0; iX < nx_ ; iX++) {
		for (unsigned iY = 0; iY < ny_; iY++) {
			for (unsigned iZ = 0; iZ < nz_; iZ++) {

				std::vector<Particle>& cell_particles = grid[iX][iY][iZ];
				for (unsigned p = 0; p < cell_particles.size(); p++) {

					if(checkIfOverlapping(Point3i(iX, iY, iZ), cell_particles[p], overlap_count)) {
						std::swap(cell_particles[p], cell_particles[cell_particles.size()-1]);	
						cell_particles.pop_back();
						p--;
					}

				}
			}
		}
	}

}
void RandomGenerator::deleteOverlappingParticles(Point3i cell1, Point3i cell2, bool same_cell, double overlap)  {

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
			if(overlap_ < 0 && fabs(overlap_) / (2.0*radius_)  >= overlap ) {
				std::swap(particles_cell2_[j], particles_cell2_[particles_cell2_.size()-1]);
				particles_cell2_.pop_back();
				j--;
			}

		}
	}
}

void RandomGenerator::deleteInGrid(double overlap) {

	int nx_ = grid.size();
	int ny_ = grid[0].size();
	int nz_ = grid[0][0].size();

	for (unsigned iX = 0; iX < nx_ ; iX++) {
		for (unsigned iY = 0; iY < ny_; iY++) {
			for (unsigned iZ = 0; iZ < nz_; iZ++) {

				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX, iY, iZ), true, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX, iY , iZ + 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX, iY - 1, iZ + 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX, iY - 1, iZ), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX, iY - 1, iZ - 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY - 1, iZ - 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY, iZ - 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY + 1, iZ - 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY - 1, iZ), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY, iZ), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY + 1, iZ), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY - 1, iZ + 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY, iZ + 1), false, overlap);
				deleteOverlappingParticles(Point3i(iX, iY, iZ), Point3i(iX + 1, iY + 1, iZ + 1), false, overlap);

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

void RandomGenerator::overlapAnalysis(std::vector<Particle>& particles) {
	double min_overlap = INT_MAX;
	double max_overlap = INT_MIN;
	double total_overlap = 0;
	int count = 0;
	for (int i =0; i< particles.size(); i++) {
		for (int j = i + 1; j< particles.size(); j++) {
			double overlap_ = calculateOverlap(particles[i], particles[j]);
			if(overlap_ < 0) {
				overlap_ = fabs(overlap_);
				if(overlap_ < 0.0001) continue;
				total_overlap += overlap_;
				count++;
				min_overlap = std::min(min_overlap,overlap_);
				max_overlap = std::max(max_overlap, overlap_);
			}
		}
	}
	
	std::cout << "max_overlap: " << max_overlap << " min_overlap: " << min_overlap << " avg_overlap: " << total_overlap / (count) << std::endl;
	std::cout << "total overlapping particles: " << count*2 << std::endl;
	std::cout << "total particles: " << particles.size() << std::endl; 
}

void RandomGenerator::randomParticleGenerator(std::vector<Particle>& particles, Domain& dom) {
	int id_ = 0;
	// double raidus_ = 0.02;
	float packingFraction_ = 0;
	double simulationVolume_ = (dom.second(0)- dom.first(0))*(dom.second(1)- dom.first(1))*(dom.second(2)- dom.first(2));
	double spheresVolume_ =0;

	// uniform real distribution

	std::random_device rd;
	
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<double> distribution_x (dom.first(0) + radius_ , dom.second(0) - radius_ );
	std::uniform_real_distribution<double> distribution_y (dom.first(1) + radius_ , dom.second(1) - radius_ );
	std::uniform_real_distribution<double> distribution_z (dom.first(2) + radius_ , dom.second(2) - radius_ );

	while (packingFraction_ < insertion_pf_) {
		particles.push_back(Particle(id_++, Point3d(distribution_x(generator), distribution_y(generator), distribution_z(generator)),radius_));
		spheresVolume_ += calculateSphereVolume(radius_);
		packingFraction_ = spheresVolume_/simulationVolume_;
	}

	ID = id_;
	std::cout << "Total Particles added into simulation domain: " << particles.size() << std::endl;
}
