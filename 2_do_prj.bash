#!/bin/bash

#python prjnishi.py \
#   --i-trj coord.dat \
#   --i-ave aveq.dat \
#   --i-egv e1.dat \
#   --o-pcc c1.dat

for i in {1..2}
do
   python prjnishi.py \
      --i-trj coord.dat \
      --i-ave aveq.dat \
      --i-egv e${i}.dat \
      --o-pcc c${i}.dat
done

paste c1.dat c2.dat > out.dat
