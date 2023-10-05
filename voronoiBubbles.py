

import numpy
import math
import scipy.interpolate
import scipy.spatial

import getFunctionsPython as gFP
import Numerical_Model as NM
import Bubble_Growth_Modelling_Function as BGMF

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from shapely.geometry import LineString
from multiprocessing import Pool, TimeoutError
import time
import os

cubeLength = 3

points = numpy.zeros((cubeLength*cubeLength*cubeLength,3))

perturbation = 0.4



for i in range(cubeLength):
    for j in range(cubeLength):
        for k in range(cubeLength):
            points[k+ cubeLength*(j) + cubeLength*cubeLength*(i)]=(i+numpy.random.uniform(-perturbation,perturbation), j+numpy.random.uniform(-perturbation,perturbation), k+numpy.random.uniform(-perturbation,perturbation))
   
x = points[:,0]

"""fig=plt.figure(3)
ax =plt.axes(projection='3d')
ax.scatter3D(points[:,0],points[:,1],points[:,2])
plt.show()"""


vor = scipy.spatial.Voronoi(points)

vorVert = vor.vertices

regionvolume=numpy.array([])
regionindices=numpy.array([])

for vorRegion in range(len(vor.regions)):
    if not (-1 in vor.regions[vorRegion]) and vor.regions[vorRegion] :

        regionvolume =  numpy.append(regionvolume, scipy.spatial.ConvexHull(vor.vertices[vor.regions[vorRegion]]).volume)
        regionindices =  numpy.append(regionindices, vorRegion)

outputs=[]

#for i in range(len(regionvolume)):
#    outputs.append(BGMF.growBubble(regionvolume[i]))


regionvolume = [1, 2, 5, 10]

if __name__ == '__main__':
    with Pool(processes=4) as pool:
        outputs = pool.map(BGMF.growBubble, regionvolume)
        try:
            print(outputs.get(timeout=1))
        except TimeoutError:
            print("We lacked patience and got a multiprocessing.TimeoutError")

        print("For the moment, the pool remains available for more work")
#pool.close()





plt.figure(1)
plt.show
