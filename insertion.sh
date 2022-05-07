#!/bin/bash

rm -rf values.csv
cd build
make
cd ..

echo "radius of particle $1 "
echo "maximum number of iteration allowed $2 "
echo "simulation domain opposite corners coordinates ($3 $4 $5), ($6 $7 $8) "
./build/insertion $1 $2 $3 $4 $5 $6 $7 $8