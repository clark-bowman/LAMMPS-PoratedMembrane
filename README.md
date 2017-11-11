# LAMMPS-PoratedMembrane

## membrane_finite.cpp
C++ script which generates a LAMMPS script `membrane_finite.in` and data file `membrane_finite.dat`. Parameters are read in from `params_membrane_finite.txt`.

## params_membrane_finite.txt
Variable parameters for simulation. Each parameter is the only entry on a line. For the order and description of parameters, see the commented parameter input section of `membrane_finite.cpp`.

## membrane_finite.in
LAMMPS script modeling a bilipid membrane in a periodic simulation box filled with DPD fluid. `membrane_finite.dat` must be in the same directory to run. This script runs with the version of LAMMPS compiled Aug. 10, 2015.

_Outputs_:
* `membrane_finite.lammpstrj` - LAMMPS trajectory file for visualization. Loadable in, e.g., [VMD](http://www.ks.uiuc.edu/Research/vmd/). Does not include fluid

## membrane_finite.dat
[Data file](http://lammps.sandia.gov/doc/read_data.html) which stores the initial positions of the lipids in a bilayer membrane. Includes positions, bonds, and angles.
