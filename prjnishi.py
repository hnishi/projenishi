# prjnishi.py v 1.0
# projection into PC space with a PCA eigenvector and average.
import sys
import numpy as np
from optparse import OptionParser

############ argument ##############
parser = OptionParser

def get_options():
   p = OptionParser()
   p.add_option('--i-trj', dest='fn_trj',
                help="file name for trajectory of structure components")
   p.add_option('--i-ave', dest='fn_ave',
                help="file name for average of structure components")
   p.add_option('--i-egv', dest='fn_egv',
                help="file name for eigen vector of VCV by PCA")
   p.add_option('--o-pcc', dest='fn_pcc',
                help="file name for output coordinates in PC space")
   opts, args = p.parse_args()
   print "----------------------------"
   p.print_help()
   print "----------------------------"
   return opts, args

############ main ##############
print "projnishi.py"
print "Projection from principle component analysis (PCA)"
opts, args = get_options()

fn_trj = 'coord.dat'
if opts.fn_trj:
   fn_trj = opts.fn_trj
f = open(fn_trj)
a = np.array(f.read().split(),dtype=float)
f.close()

fn_ave = 'aveq.dat'
if opts.fn_ave:
   fn_ave = opts.fn_ave
f = open(fn_ave)
a2 = np.array(f.read().split(),dtype=float)
f.close()

fn_egv = 'e1.dat'
if opts.fn_egv:
   fn_egv = opts.fn_egv
f = open(fn_egv)
a3 = np.array(f.read().split(),dtype=float)
f.close()

if len(a2) != len (a3):
   print "ERROR: the num of date in fn_ave and fn_egv are not the same"
   sys.exit(1)
dim = len(a3)
structure = len(a)/dim
print "dimension: ",dim
print "the number of structures: ", structure

##### CALCULATE INNER PRODUCT  
ccc = a.reshape(len(a)/dim,dim)
#print np.shape(ccc)

dots = [] #dots:PC coordinates
for i in xrange(len(a)/dim):
   #print ccc[i] - a2
   #print np.dot(a3,ccc[i] - a2)
   dots.append(np.dot(a3,ccc[i] - a2)) 
   
#print np.shape(dots)
#print dots 

##### OUTPUT PC COORDINATES OF ALL STRUCTURES  
fn_pcc = 'c1.dat'
if opts.fn_pcc:
   fn_pcc = opts.fn_pcc
output = open(fn_pcc,'w')
for i in dots:
   print >> output, i
output.close()

