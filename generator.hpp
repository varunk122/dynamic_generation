#pragma once

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

class RandomGenerator {
public:

	std::vector<std::vector<std::vector<std::vector<Particle>>>> grid; //single-hierarchy grid which will store particles for roboust deletion algorithm

	RandomGenerator();

	void randomParticleGenerator(std::vector<Particle>& particles, Domain& dom);

	void addParticlesToGrid(std::vector<Particle>& particles, Domain& dom);
	void addParticlesToList(std::vector<Particle>& particles);
	void deleteParticles(std::vector<Particle>& particles);
	void deleteOverlappingParticles(Point3d cell1, Point3d cell2, bool same_cell);
	void deleteInGrid();
	void randomDeleting(std::vector<Particle>& particles, Domain& dom, double required_pf);

	double calculateOverlap(Particle& particle1, Particle& particle2);
	double calculatePackingFraction(std::vector<Particle>& particles, Domain & dom);
	double calculateSphereVolume(double& r);
	void updatePackingFraction(double& packing_fraction, Domain& dom, Particle& deleted_particle);

	void printParticlesList(std::vector<Particle>& particles);
	void printParticle(Particle& particle);
	void printGrid();
	void saveToCSV(std::vector<Particle>& particles);
};
