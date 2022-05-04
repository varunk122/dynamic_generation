# dynamic_generation

## Dependencies

* **c++ eigen library** - for linear algebra  ```sudo apt-get install libeigen3-dev```
  
## Build

To build the code run the following commands -

```
mkdir build 
cd build
make
cd ..
chmod +x run.sh
```

# Parameters

* **simulation domain** - specify coordinates of two opposite corners of simulation domain in m
* **radius of particles** - radius of monodisperse particles in m 
* **insertion packing fraction** - set packing fraction for insertion algorithm


## Run
```
./run.sh radius insertion_packing_fraction coordinate_1 coordinate_2
```

## Example 

```./run.sh 0.5 5 0 0 0 10 10 10```

This command will run the algorithm for cubical simulation domain of edge size 10m, monodisperse particle raidus of size 0.5m with insertion packing fraction of 5.

## Output

* values.csv file will be generated in which contains postion and radius of all non-overlapping partilces.
* A matlab script **plotting.m** is provided for visualization.
