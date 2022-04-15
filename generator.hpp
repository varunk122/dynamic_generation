#pragma once

#include<bits/stdc++.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include<fstream>


typedef Eigen::Vector3d Point3d;
typedef Eigen::Vector3i Point3i;

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
	unsigned ID;
	double radius_ = 0.05;

	RandomGenerator();

	void randomParticleGenerator(std::vector<Particle>& particles, Domain& dom);

	void addParticlesToGrid(std::vector<Particle>& particles, Domain& dom);
	void addParticlesToList(std::vector<Particle>& particles);
	void deleteParticles(std::vector<Particle>& particles);
	void deleteOverlappingParticles(Point3i cell1, Point3i cell2, bool same_cell, double overlap);
	void deleteInGrid(double overlap);
	void randomDeleting(std::vector<Particle>& particles, Domain& dom, double required_pf);

	double calculateOverlap(Particle& particle1, Particle& particle2);
	double calculatePackingFraction(std::vector<Particle>& particles, Domain & dom);
	double calculateSphereVolume(double& r);
	void updatePackingFraction(double& packing_fraction, Domain& dom, Particle& deleted_particle);

	void printParticlesList(std::vector<Particle>& particles);
	void printParticle(Particle& particle);
	void printGrid();
	void saveToCSV(std::vector<Particle>& particles);

	void newDeleteOverlappingParticles(int overlap_count);
	bool checkIfOverlapping(Point3i cell, Particle particle, int overlap_count);
	int checkIfOverlappingWithParticlesInCell(Point3i cell, Particle particle);
	unsigned emptyCellsCount();
	unsigned totalCellsInGrid();
	void addParticlesInEmptyCells();
	void overlapAnalysis(std::vector<Particle>& particles);

};
