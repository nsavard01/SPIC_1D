This is a 1D electrostatic Particle-in-Cell (PIC) code. The S stands for Savard, because I am creative. This is a re-written version of Fortran code written in C++ due to its inherent ability to use Object-oriented-program (Fortran does have this, but it is not a foundation of the language), and its continued popularity, which I hope will somewhat future-proof the code. It includes explicit formulations (momentum-conserving and energy-conserving), as well as implicit energy-and-charge conserving to the lowest degree of interpolation (nearest-grid-point field to particle, linear particle to charge density). The following describes a rough outline of what is needed to run the code in its current iteration. Now that I no longer live in the academic bubble, it may be a while until this code is touched, but the hope is that it will be an appropriate base for graduate students in the future.

You'll need the Intel oneAPI HPC Toolkit installed on a Linux platform. I have found that everything discussed here also works well with WSL on windows (Linux subsystem). The toolkit includes the icpx and mpiicpx compilers. It also includes mkl and openmp libraries used for some mathematical functions, which can be conveniently linked with -qmkl and -qopenmp. The mpiicpx compiler (the default) allows the program to use the mpi libraries. You will also need CMake for compilation of the source files.

One can use other C++ compilers and presumably change the CMakeLists.txt, but then you'll need to have appropriate linking to all libraries. I am not well versed in this area and therefore choose the convenient method that allows me to remain ignorant of the inner workings of compilation. 


-------------
COMPILATION
-------------

Apparently a build folder is a good idea to separate all your compilation files from the source code, so I will follow this standard naively. First, download the SPIC_1D files at "https://github.com/nsavard01/SPIC_1D" locally to your computer. Now enter the SPIC_1D directory:

cd ./(intermediate-paths-to-SPIC_1D)/SPIC_1D

Now create a build folder within this directory, and enter it:

mkdir build
cd build

Now you let CMake do its thing:

cmake ..

You should get notification on your terminal that build files have been written






Folder "particle_operations" includes any operations with charged particles not included in the standard time stepping procedure for fields, or collisions. Example does not use this, but an example operation that can be included is injection from a wall. This would be in a file "wall_injection.inp", and an example for electrons is:

-------------------------------------
e : particle_name
0 : wall node to release from
50.0 : current density (A/m^2)
0.0 : x component of velocity (m/s)
0.0 : y component of velocity (m/s)
0.0 : z component of velocity (m/s)
0.1 : Flux temperature (eV)
-----------------------------------

END

The particle name must match the particle name given in the particle inputs.
 





