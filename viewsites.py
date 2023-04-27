#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

fname1='BVECTORS';
fname2='PLANE';
fname3='points.dat';

fid = open(fname1, 'r');
dim = np.array([int(i) for i in fid.readline().split()]);

vectors = np.array([]);
for i in range(dim[0]):
    vectors = np.append(vectors, np.array([float(j) for j in fid.readline().split()]));
vectors = np.reshape(vectors, dim);
#print (vectors);

plt.figure("bonds");

points = np.loadtxt(fname3);
for point in points:
    x1 = np.matmul(point[0:dim[0]],vectors);
    x2 = np.matmul(point[dim[0]:],vectors);
    plt.plot([x1[0],x2[0]], [x1[1],x2[1]], '-k');
    
plt.axis('scaled');

plt.figure("sites");
for point in points:
    x1 = np.matmul(point[0:dim[0]],vectors);
    plt.plot(x1[0],x1[1], '.k');
plt.axis('scaled');
plt.show();
