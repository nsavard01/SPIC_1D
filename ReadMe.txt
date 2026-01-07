This is a 1D electrostatic Particle-in-Cell (PIC) code. The S stands for Savard, because I am creative. This is a re-written version of Fortran code written in C++ due to its inherent ability to use Object-oriented-program (Fortran does have this, but it is not a foundation of the language), and its continued popularity, which I hope will somewhat future-proof the code. It includes explicit formulations (momentum-conserving and energy-conserving), as well as implicit energy-and-charge conserving to the lowest degree of interpolation (nearest-grid-point field to particle, linear particle to charge density). The following describes a rough outline of what is needed to run the code in its current iteration. Now that I no longer live in the academic bubble, it may be a while until this code is touched, but the hope is that it will be an appropriate base for graduate students in the future. Furthermore, time is now limited for working on this, including documentation. So while I appreciate that documentation here is total crap, I simply do not have time or motivation to make it professional-grade.

-------------
COMPILATION
-------------

You will need the Intel oneAPI HPC Toolkit installed on a Linux platform. I have found that everything discussed here also works well with WSL on windows (Linux subsystem). The toolkit includes the icpx and mpiicpx compilers. It also includes mkl and openmp libraries used for some mathematical functions, which can be conveniently linked with -qmkl and -qopenmp. The mpiicpx compiler (the default) allows the program to use the mpi libraries. You will also need CMake for compilation of the source files.

One can use other C++ compilers and presumably change the CMakeLists.txt, but then you'll need to have appropriate linking to all libraries. I am not well versed in this area and therefore choose the convenient method that allows me to remain ignorant of the inner workings of compilation. 

Apparently a build folder is a good idea to separate all your compilation files from the source code, so I will follow this standard naively. First, download the SPIC_1D files at "https://github.com/nsavard01/SPIC_1D" locally to your computer. Now enter the SPIC_1D directory:

$ cd ./(intermediate-paths-to-SPIC_1D)/SPIC_1D

Now create a build folder within this directory, and enter it:

$ mkdir build
$ cd build

Now you let CMake do its thing:

$ cmake ..

You should get notification on your terminal that build files have been written. Now its time to make the compiled program. I believe make is a command that all linux systems come with, but if this is not the case, make sure you download it.

$ make

Now you should have in your build directory a file pic1d. If so, congrats, you have successfully compiled the program.


-------------
FILE INPUTS
-------------

Before getting into running the program, let's first define the inputs to the program, which are in the "inputs" folder. The "initial_setup.inp" file describes high-level parameters such as amount of openmp threads, simulation time, and number of diagnostic data dumps. "geometry.inp" is to setup the domain. "implicit_solver.inp" is only used if the implicit scheme is used, since it requires solving a non-linear equation. "charged_particles" directory is where you place .inp files which describe each charged particle within the simulation. "target_particles" describes each background particle used for background monte-carlo collisions. Files describing collision cross-sections between charged particles and target particles are described in "collisions/binary" folder, also as an .inp file. The format follows those used in LXCAT. Finally, there is an additional folder "particle_operations" which includes any operations with charged particles not included in the standard time stepping procedure for fields or MC collisions.

The example does not use this, but an example operation that can be included is injection from a wall. This would be in a file "wall_injection.inp", and an example for electrons is:

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

-------------
EXAMPLE RUN
-------------

The current SPIC_1D folder has inputs in order to run the classic Turner capacitive-couple RF plasma. Before running, two things that may need to be considered for your case in the "initial_setup.inp" file. It's a good rule of thumb (with some exceptions) to have the same number of openmp threads (line 1) as cpu cores. So change this line as required. Second, you need to have a directory to store the output files from the program, which is the line 8. This is a relative path from where you run the program, which is in the build folder. At the moment, I have SPIC_test directory two steps up the directory tree. You can choose to have the directory where ever you want, but make sure to have the right path to the folder.

Now to run the program while in the "build" directory, you will do:

$ mpirun -np N ./pic1D

N = number of mpi ranks. Another good rule of thumb: divide up the mpi ranks into number of total sockets x number of machines. This is because the program is designed so that you maximize computation within a rank using openmp threads (shared memory). Because transferring data between different sockets and machines comes with delay (no shared caches), you minimize the advantage of openmp threads shared across these components. So it's best to separate these components into different MPI ranks which perform parallel computations within each rank, and then transfer information between them when needed. Keep in mind line 1 in the "initial_setup.inp" file is the amount of openmp threads PER MPI RANK, so the total amount of openmp threads (which should be ~ number of cores) is N x # openmp threads.

If you want the program to run even if terminal is killed for whatever reason, use:

$ nohup mpirun -np N ./pic1D > output.txt &

This will run in the background, and place the outputs into a output.txt file.

Now after running, you should see a bunch of texts which describes the initial conditions, and then outputs for the diagnostics and averaging. Now within your "SPIC_test" directory you should have a directory "RF_Turner_test" which contains all the output files.


-------------
EXAMPLE ANALYSIS
-------------

In the "scripts" folder, there are python files. These describe classes used to collect/analyze the data, along with plotting functions. An example analysis file, named "example_analysis.py", is also contained in this file. Copy and place within your "SPIC_test" folder. Running this example will load the data into a data set object, and then show examples to plot the averaged densities, electric potential, and temperatures. Any errors may happen due to not having the right packages.
