#!/bin/bash

g++ extr_dista.cpp quatnishi.cpp pdbnishi.cpp tranishi.cpp inpnishi.cpp math_nishi.cpp -o extr_dista.exe -I ~/Dropbox/software/eigen-eigen-1306d75b4a21

./extr_dista.exe prm_dista.inp > res_dista.res 

cat res_dista.res


