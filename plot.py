import numpy as np
import pandas as pd
import csv
import os
from mayavi import mlab

# file_path = os.path

f = open("./values.csv", "r")
x = []
y = []
z = []
r = []

for p in csv.reader(f):
    x.append(float(p[0]))
    y.append(float(p[1]))
    z.append(float(p[2]))
    r.append(float(p[3]))   


mlab.points3d(x,y,z,r)
mlab.show()
