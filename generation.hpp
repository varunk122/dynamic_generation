#include<bits/stdc++.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include<fstream>


typedef Eigen::Vector3d Point3d;
#define PI 3.14159265358979323846

struct Particle {
	Point3d centre;
	double radius;
	unsigned id;
	Particle (unsigned id_, Point3d centre_, double radius_){
		id = id_;
		centre = Point3d(centre_);
		radius = radius_;
	}
};

struct Domain {
	Point3d first;
	Point3d second;
	Domain (Point3d first_ , Point3d second_) {
		first = Point3d(first_);
		second = Point3d(second_);
	}
};

std::vector<Particle> particles;

std::vector<std::vector<std::vector<int>>> grid;

double calculateSphereVolume(double& r) {
	return (4.0*PI*pow(r,3)) / 3.0;
}

double calculatePackingFraction(std::vector<Particle>& particles, Domain & dom) {
	double spheres_vol_ = 0;
	for(int i =0; i< particles.size(); i++) {
		spheres_vol_ += calculateSphereVolume(particles[i].radius);
	} 

	double simulation_vol_ = (dom.second(0)- dom.first(0))*(dom.second(1)- dom.first(1))*(dom.second(2)- dom.first(2));

	return spheres_vol_ / simulation_vol_;

}

void updatePackingFraction(double& packing_fraction, Domain& dom, Particle& deleted_particle) {

	double simulation_vol_ = (dom.second(0)- dom.first(0))*(dom.second(1)- dom.first(1))*(dom.second(2)- dom.first(2));
	packing_fraction = (packing_fraction * simulation_vol_ - calculateSphereVolume(deleted_particle.radius)) / simulation_vol_;

}

void printParticlesList(std::vector<Particle>& particles) {
	
	for (int i =0; i< particles.size();i++) {
		std::cout << particles[i].centre(0) << " " << particles[i].centre(1) << " " << particles[i].centre(2) << std::endl;
	}
}

void saveToCSV(std::vector<Particle>& particles) {

	std::fstream particle_file_;
	particle_file_.open("./values.csv", std::ios::out | std::ios::app);

	for (int i =0; i< particles.size();i++) {
		particle_file_ << particles[i].centre(0) << ", " << particles[i].centre(1) << ", " << particles[i].centre(2) << ", "
						<< particles[i].radius << "\n";
	}

}

double calculateOverlap(Point3d& point1, Point3d& point2) {
	double dist_ = sqrt(pow((point1(0) - point2(0)),2) + pow((point1(1) - point2(1)),2) + pow((point1(2) - point2(2)),2));
	return dist_ - 0.4;
}

void deleteParticles(std::vector<Particle>& particles) {
	std::vector<unsigned> toDelete_;
	for (int i =0; i < particles.size(); i++) {
		for (int j =i+1; j <particles.size(); j++) {
			double overlap_ = calculateOverlap(particles[i].centre, particles[j].centre);
			if( overlap_ < 0 && fabs(overlap_)/0.4 > 0.2) {
				particles.erase(particles.begin() + j);
				j--;
			}
		}

	}
}

void randomParticleGenerator(std::vector<Particle>& particles, Domain& dom) {
	int id_ = 0;
	int raidus_ = 0.2;
	float packingFraction_ = 0;
	double simulationVolume_ = (dom.second(0)- dom.first(0))*(dom.second(1)- dom.first(1))*(dom.second(2)- dom.first(2));
	double spheresVolume_ =0;

	// uniform real distribution

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution_x (dom.first(0), dom.second(0));
	std::uniform_real_distribution<double> distribution_y (dom.first(1), dom.second(1));
	std::uniform_real_distribution<double> distribution_z (dom.first(2), dom.second(2));

	while (packingFraction_ < 10) {
		particles.push_back(Particle(id_++, Point3d(distribution_x(generator), distribution_y(generator), distribution_z(generator)),0.2));
		spheresVolume_ += 0.0334933;
		packingFraction_ = spheresVolume_/simulationVolume_;
	}
}




