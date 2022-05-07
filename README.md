# dynamic_generation

## Dependencies

* **c++ eigen library** - for linear algebra  ```sudo apt-get install libeigen3-dev```
  
## Build

To build the code run the following commands -

```
mkdir build 
cd build
cmake ..
make
cd ..
chmod +x deletion.sh insertion.sh
```

# Parameters

For deletion algorithm - 

* **simulation domain** - specify coordinates of two opposite corners of simulation domain in m
* **radius of particles** - radius of monodisperse particles in m 
* **insertion packing fraction** - set packing fraction for insertion algorithm

For insetion algorithm - 

* **simulation domain** - specify coordinates of two opposite corners of simulation domain in m
* **radius of particles** - radius of monodisperse particles in m 
* **maximum_iterations** - maximum iteration in which a new particle must be inserted in the simulation domain



## Run

* Deletion - 
```
./deletion.sh radius insertion_packing_fraction coordinate_1 coordinate_2
```
* Insertion - 
```
./insetion.sh radius maximum_iterations coordinate_1 coordinate_2
```
## Example 

* Deletion - 

```./deletion.sh 0.5 5 0 0 0 10 10 10```

This command will run the deletion algorithm for cubical simulation domain of edge size 10m, monodisperse particle radius of size 0.5m with insertion packing fraction of 5.

* Insertion - 

```./insertion.sh 0.5 10000 0 0 0 10 10 10```

This command will run the insertion algorithm for cubical simulation domain of edge size 10m, monodisperse particle radius of size 0.5m with maximum iterations allowed for the insertion of particles are 10,000.

## Output

* values.csv file will be generated in both the cases which contain positions and radius of all non-overlapping particles.
* A MATLAB script **plotting.m** is provided for visualization.
