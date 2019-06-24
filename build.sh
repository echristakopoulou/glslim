#!/bin/bash

cd bcls-0.1/ 
./configure 
make 
cd .. 
mkdir BCLSlib 
cd BCLSlib 
cp ../bcls-0.1/include/* . 
cp ../bcls-0.1/src/.libs/*a . 
cd .. 
mkdir build 
cd build 
cmake .. 
make
