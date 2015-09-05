#!/bin/bash

#g++ extr_coord.cpp quatnishi.cpp pdbnishi.cpp tranishi.cpp inpnishi.cpp math_nishi.cpp -o extr_coord.exe -I ~/Dropbox/software/eigen-eigen-1306d75b4a21
g++ extr_coord.cpp quatnishi.cpp pdbnishi.cpp tranishi.cpp inpnishi.cpp math_nishi.cpp -o extr_coord.exe -I ~/nishitool/eigen-eigen-1306d75b4a21

time ./extr_coord.exe prm_extr.inp > res_extr.res 

cat res_extr.res 

numst=$(wc coord.dat | awk '{print $1}')
numall=$(wc coord.dat | awk '{print $2}')
echo The number of structures = $numst

numco=$(python -c "print  ${numall}/${numst}")
echo The number of components = $numco



