#include<bits/stdc++.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>


typedef Eigen::Vector3d Point3d;

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

double calculatePackingFraction(std::vector<Particle>& particles, Domain & dom);

void printParticlesList(std::vector<Particle>& particles) {
	
	for (int i =0; i< particles.size();i++) {
		std::cout << particles[i].centre(0) << " " << particles[i].centre(1) << " " << particles[i].centre(2) << std::endl;
	}
}

double calculateOverlap(Point3d& point1, Point3d& point2) {
	double dist_ = sqrt(pow((point1(0) - point2(0)),2) + pow((point1(1) - point2(1)),2) + pow((point1(2) - point2(2)),2));
	return dist_ - (point1.radius + point2.radius);
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

	while (packingFraction_ < 3) {
		particles.push_back(Particle(id_++, Point3d(distribution_x(generator), distribution_y(generator), distribution_z(generator)),0.2));
		spheresVolume_ += 0.0334933;
		packingFraction_ = spheresVolume_/simulationVolume_;
	}
}




