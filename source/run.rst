Compiling and running iPIC3D
==============================

To install and run iPIC3D, you need:

* cmake (at least version 2.8);
* MPI;
* HDF5 (optional).

Follow these isntructions: ::

  #create a build directory and go there
  mkdir build && cd build
  cmake /path/to/where/CMakeLists.txt/located

  #compile
  #if successful, you will find an executable named iPIC3D in the build directory
  make

  #run correctness tests automatically
  ctest

  #run a test case manually
  #copy an inputfile named as testXXX.inp from ../inputfiles to the build directory
  #make sure you create a folder for the output data as specified in the input file
  #make sure no_of_proc = XLEN x YLEN x ZLEN as specified in the input file
  mpiexec -n no_of_proc ./iPIC3D ./inputfilename.inp

