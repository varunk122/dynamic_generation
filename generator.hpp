#pragma once

#include<bits/stdc++.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include<fstream>


typedef Eigen::Vector3d Point3d;
typedef Eigen::Vector3i Point3i;

#define PI 3.14159265358979323846

//custom data type for particles
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

//custom data type for simulation domain
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

	std::vector<std::vector<std::vector<std::vector<Particle>>>> grid; // single-hierarchy grid which will store particles for roboust deletion algorithm
	unsigned ID;
	double radius_ = 0.5; //default value of radius of particles
	double insertion_pf_ = 5; // default value of insetion packing fraction

	RandomGenerator(); //constructor

	void randomParticleGenerator(std::vector<Particle>& particles, Domain& dom); //insetion algorithm

	void addParticlesToGrid(std::vector<Particle>& particles, Domain& dom); //add particles from the list to single-hierarchy grid
	void addParticlesToList(std::vector<Particle>& particles); // add particles from grid to list
	void deleteParticles(std::vector<Particle>& particles); // brute-froce algorithm to delete all particles that are overlapping with a particluar particle
	void deleteOverlappingParticles(Point3i cell1, Point3i cell2, bool same_cell, double overlap); // helper function for algorithm-1
	void deleteInGrid(double overlap); //delete particles using algorithm-1
	void randomDeleting(std::vector<Particle>& particles, Domain& dom, double required_pf); // randomly delete particles till the required packing fraction is achieved

	double calculateOverlap(Particle& particle1, Particle& particle2); // calculate overlap between two particles
	double calculatePackingFraction(std::vector<Particle>& particles, Domain & dom); //calculate packing fraction 
	double calculateSphereVolume(double& r); // calculate volume of a sphere
	void updatePackingFraction(double& packing_fraction, Domain& dom, Particle& deleted_particle); //to update packing fraction after deletion of a particle

	void printParticlesList(std::vector<Particle>& particles); // print all particles present in the list
	void printParticle(Particle& particle); //print a single particle
	void printGrid(); //print all particles present in the grid
	void saveToCSV(std::vector<Particle>& particles); //add particles to a csv file

	void newDeleteOverlappingParticles(int overlap_count); // delete particles using algorithm-2
	bool checkIfOverlapping(Point3i cell, Particle particle, int overlap_count); //helper function for algorithm-2
	int checkIfOverlappingWithParticlesInCell(Point3i cell, Particle particle); //helper function for algorithm-2
	unsigned emptyCellsCount(); // count number of empty cells in a grid
	unsigned totalCellsInGrid(); // count number of total cells in a grid
	void addParticlesInEmptyCells(); // re-add particles at the centre of empty grid
	void overlapAnalysis(std::vector<Particle>& particles); // provide maximum, minimum and average overlap

};
