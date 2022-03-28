#! /bin/bash

gfortran msd_polymer_droplet.f90 -o msd_droplet
./msd_droplet
#gfortran msd_polymer_free.f90 -o msd_free_poly
#./msd_free_poly
rm -rf temp.dat 
