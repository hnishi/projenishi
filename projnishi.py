# projnishi.py v 1.0
#projection into PC space with a PCA eigenvector and average.
import numpy as np

structure = 1002
dim = 63

#vec_coordinates = vec_all_coordinates

f = open('coord.dat')
a = np.array(f.read().split(),dtype=float)
f.close()

f = open('aveq.dat')
a2 = np.array(f.read().split(),dtype=float)
f.close()

f = open('e1.dat')
a3 = np.array(f.read().split(),dtype=float)
f.close()

ccc = a.reshape(len(a)/63,63)
#print np.shape(ccc)

dots = []
for i in xrange(len(a)/63):
   #ddd[i] = ccc[i] - a3
   dots.append(np.dot(a3,ccc[i] - a2))
   
#print np.shape(dots)
print dots

"""
for i in xrange(strcture):
   vec_average_coordinates = vstack( vec_average_coordinates, pre_vec_average_coordinates )

a = np.array(vec_coordinates[i*dim] - vec_average_coordinates)

vec_motometai = vec_eig*(vec_coordinates - vec_average_coordinates)
"""

