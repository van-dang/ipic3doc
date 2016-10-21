Enabling OpenMP threading support in MPI in iPIC3D
==================================================

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

As the first step to port the original version of the iPIC3D code to the OpenMP threading version, the MPI initialisation routine has been changed: in utility/MPIdata.cpp the function MPI_Init() was replaced by MPI_Init_thread(), and now it looks as ::

   MPI_Init_thread( 0, 0, MPI_THREAD_MULTIPLE, &provided )

This allowed to initialise the MPI execution environment. MPI_Init_thread() has two additional parameters for the level of thread support required, and for the level of thread support provided by the MPI library implementation. The general C/C++ syntax of this function is: ::

   #include <mpi.h>
   int MPI_Init_thread(int *argc, char ***argv, int required, int *provided)

“required” specifies the requested level of thread support, and the actual level of support is returned into “provided”. After the initialization, a programmer may add OpenMP directives and runtime calls as long as he/she sticks to the level of thread safety that was specified in the call to MPI_Init_thread().

In order to get the support from MPI that allows any thread to call the MPI library at any time, one has to have MPI_THREAD_MULTIPLE for the “required” variable and then the value of the “provided” variable has to be equal to 3.

There are three another options for thread support that can be used in the MPI_Init_thread() function, they are described below.

* MPI_THREAD_SINGLE

  * Only one thread will execute. This represents a typical MPI-only application.

* MPI_THREAD_FUNNELED

  * The process may be multi-threaded, but only the master thread will make calls to the MPI library (all MPI calls are funneled to the main thread). The thread that calls MPI_Init_thread() is from now on the main thread. To check whether a thread is the main thread, one may call the MPI_Is_thread_main() routine.

* MPI_THREAD_SERIALIZED

  * The process may be multi-threaded, and multiple threads may make MPI calls, but only one at a time: MPI calls are not made concurrently from two distinct threads, thus all MPI calls are serialized.

With the MPI_THREAD_MULTIPLE option, any thread within any process is eligible to call MPI functions, and the MPI library is responsible for thread safety within that library and for libraries that it uses. Figure 1 illustrates two cases of the MPI and OpenMP programming models interoperability:

* without OpenMP thread support in an MPI library (on the left side of the figure);
* with enabling OpenMP thread support in an MPI library (on the right side of the figure).

.. figure:: threading.png

   Figure 1. OpenMP threading in an MPI code

In order to select a thread level higher than MPI_THREAD_SINGLE, a programmer should set the environment variable MPICH_MAX_THREAD_SAFETY, otherwise “provided” will return MPI_THREAD_SINGLE. For the case of MPI_THREAD_MULTIPLE, it should be as ::

   export MPICH_MAX_THREAD_SAFETY=multiple

In iPIC3D, halo exchange consists of three communication phases:

1. communication between faces;
2. communication between edges;
3. communication between corners.

The first phase (communication between faces) includes six sequential calls to MPI_Irecv() and six sequential calls to MPI_Isend(). With the usage of MPI_THREAD_MULTIPLE, it is possible to parallelise each set of these calls. As there are only six independent serial calls to the non-blocking MPI receive functions and six – to the non-blocking MPI send functions, the maximal possible level of parallelism here is six. And therefore, with taking into account that the number of CPUs in a node of a supercomputer is a power of two, there are several possible strategies for the parallelization of the region (with two, four and eight threads). Table below describes these approaches.

+----------------------+-----------------------------------------------------------------------------+
|No MPI_THREAD_MULTIPLE|6 MPI_Irecv() functions and 6 MPI_Isend() functions are called sequentially. |
+----------------------+-----------------------------------------------------------------------------+
|1 MPI process has     |Within 1 MPI process, 2 OpenMP threads call MPI_Irecv() functions in parallel|
|2 OpenMP threads      |and then MPI_Isend() functions in parallel.                                  |
+----------------------+-----------------------------------------------------------------------------+
|1 MPI process has     |Within 1 MPI process, 4 OpenMP threads call MPI_Irecv() functions in parallel|
|4 OpenMP threads      |and then MPI_Isend() functions in parallel, and thereafter 2 OpenMP threads  |
|                      |do the same with MPI_Isend() and MPI_Irecv() functions.                      |
+----------------------+-----------------------------------------------------------------------------+
|1 MPI process has     |Within 1 MPI process, every OpenMP thread (except 2) make calls to           |
|8 OpenMP threads      |MPI_Irecv() in parallel and then to MPI_Isend() in parallel.                 |
+----------------------+-----------------------------------------------------------------------------+

Listing below shows halo exchange in iPIC3D for the phase of communication between faces, without MPI_THREAD_MULTIPLE, so threads cannot call an MPI function, and therefore all MPI function calls have to be sequential ::

  if(left_neighborX != MPI_PROC_NULL && left_neighborX != myrank)
  {
    MPI_Irecv(&vector[0][1][1], 1, yzFacetype, left_neighborX,
    tag_XR, comm, &reqList[recvcnt++]);
    communicationCnt[0] = 1;
  }
  if(right_neighborX != MPI_PROC_NULL && right_neighborX != myrank)
  {
    MPI_Irecv(&vector[nx-1][1][1], 1, yzFacetype, right_neighborX,
    tag_XL, comm, &reqList[recvcnt++]);
    communicationCnt[1] = 1;
  }
  if(left_neighborY != MPI_PROC_NULL && left_neighborY != myrank)
  {
    MPI_Irecv(&vector[1][0][1], 1, xzFacetype, left_neighborY,
    tag_YR, comm, &reqList[recvcnt++]);
    communicationCnt[2] = 1;
  }
  if(right_neighborY != MPI_PROC_NULL && right_neighborY != myrank)
  {
    MPI_Irecv(&vector[1][ny-1][1], 1, xzFacetype, right_neighborY,
    tag_YL, comm, &reqList[recvcnt++]);
    communicationCnt[3] = 1;
  }
  if(left_neighborZ != MPI_PROC_NULL && left_neighborZ != myrank)
  {
    MPI_Irecv(&vector[1][1][0], 1, xyFacetype, left_neighborZ,
    tag_ZR, comm, &reqList[recvcnt++]);
    communicationCnt[4] = 1;
  }
  if(right_neighborZ != MPI_PROC_NULL&& right_neighborZ != myrank)
  {
    MPI_Irecv(&vector[1][1][nz-1], 1, xyFacetype, right_neighborZ,
    tag_ZL, comm, &reqList[recvcnt++]);
    communicationCnt[5] = 1;
  }
  sendcnt = recvcnt;
  int offset = (isCenterFlag ?0:1);
  if(communicationCnt[0] == 1)
  {
    MPI_Isend(&vector[1+offset][1][1], 1, yzFacetype, left_neighborX,
    tag_XL, comm, &reqList[sendcnt++]);
  }
  if(communicationCnt[1] == 1)
  {
    MPI_Isend(&vector[nx-2-offset][1][1], 1, yzFacetype, right_neighborX,
    tag_XR, comm, &reqList[sendcnt++]);
  }
  if(communicationCnt[2] == 1)
  {
    MPI_Isend(&vector[1][1+offset][1], 1, xzFacetype, left_neighborY,
    tag_YL, comm, &reqList[sendcnt++]);
  }
  if(communicationCnt[3] == 1)
  {
    MPI_Isend(&vector[1][ny-2-offset][1], 1, xzFacetype, right_neighborY,
    tag_YR, comm, &reqList[sendcnt++]);
  }
  if(communicationCnt[4] == 1)
  {
    MPI_Isend(&vector[1][1][1+offset], 1, xyFacetype, left_neighborZ,
    tag_ZL, comm, &reqList[sendcnt++]);
  }
  if(communicationCnt[5] == 1)
  {
    MPI_Isend(&vector[1][1][nz-2-offset], 1, xyFacetype, right_neighborZ,
    tag_ZR, comm, &reqList[sendcnt++]);
  }
  assert_eq(recvcnt,sendcnt-recvcnt);
  …

Another listing below illustrates the same case in iPIC3D, but with enabling OpenMP thread support for MPI for the case with four OpenMP threads, so with MPI_THREAD_MULTIPLE four threads call MPI functions in parallel. As there are six possible parallel calls to the MPI_Irecv() functions and another six – to the MPI_Isend() functions, four calls to MPI_Irecv() are executed by the threads in parallel, and then the remain two calls to MPI_Irecv() are done in parallel. Same strategy is applied to MPI_Isend(). ::

  #pragma omp parallel default(shared) private(id_thread,nthreads)
  {
    id_thread = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    …
    if (nthreads==4)
    {
      if(id_thread==0)
        if(left_neighborX != MPI_PROC_NULL && left_neighborX != myrank )
        {
          #pragma omp critical
          {
            MPI_Irecv(&vectorMPI_Isend(&vector[1][1][[0][1][1], 1,
            yzFacetype, left_neighborX,tag_XR, comm, &reqList[recvcnt]);
            recvcnt++;
            communicationCnt[0] = 1;
          }
        }
      if(id_thread==1)
        if(right_neighborX != MPI_PROC_NULL && right_neighborX != myrank )
        {
          #pragma omp critical
          {
            MPI_Irecv(&vector[nx-1][1][1], 1, yzFacetype,
            right_neighborX,tag_XL, comm, &reqList[recvcnt]);
            recvcnt++;
            communicationCnt[1] = 1;
          }
        }
      if(id_thread==2)
        if(left_neighborY != MPI_PROC_NULL && left_neighborY != myrank )
        {
          #pragma omp critical
          {
            MPI_Irecv(&vector[1][0][1], 1, xzFacetype, left_neighborY,
            tag_YR, comm, &reqList[recvcnt]);
            recvcnt++;
            communicationCnt[2] = 1;
          }
        }
      if(id_thread==3)
        if(right_neighborY != MPI_PROC_NULL && right_neighborY != myrank )
        {
          #pragma omp critical
          {
            MPI_Irecv(&vector[1][ny-1][1], 1, xzFacetype,
            right_neighborY,tag_YL, comm, &reqList[recvcnt]);
            recvcnt++;
            communicationCnt[3] = 1;
          }
        }
      if(id_thread==0)
        if(left_neighborZ != MPI_PROC_NULL && left_neighborZ != myrank )
        {
          #pragma omp critical
          {
            MPI_Irecv(&vector[1][1][0],1, xyFacetype,
            left_neighborZ,tag_ZR, comm, &reqList[recvcnt]);
            recvcnt++;
            communicationCnt[4] = 1;
          }
        }
      if(id_thread==1)
        if(right_neighborZ != MPI_PROC_NULL&& right_neighborZ != myrank )
        {
          #pragma omp critical
          {
            MPI_Irecv(&vector[1][1][nz-1], 1, xyFacetype,
            right_neighborZ,tag_ZL, comm, &reqList[recvcnt]);
            recvcnt++;
            communicationCnt[5] = 1;
          }
        }
      #pragma omp barrier
      if(id_thread==0)
      {
        sendcnt = recvcnt;
        offset = (isCenterFlag ?0:1);
      }
      #pragma omp barrier
      if(id_thread==0)
        if(communicationCnt[0] == 1)
        {
          #pragma omp critical
          {
            MPI_Isend(&vector[1+offset][1][1],1, yzFacetype,
            left_neighborX, tag_XL, comm, &reqList[sendcnt]);
            sendcnt++;
          }
        }
      if(id_thread==1)
        if(communicationCnt[1] == 1)
        {
          #pragma omp critical
          {
            MPI_Isend(&vector[nx-2-offset][1][1], 1, yzFacetype,
            right_neighborX,tag_XR, comm, &reqList[sendcnt]);
            sendcnt++;
          }
        }
      if(id_thread==2)
        if(communicationCnt[2] == 1)
        {
          #pragma omp critical
          {
            MPI_Isend(&vector[1][1+offset][1],1, xzFacetype,
            left_neighborY, tag_YL, comm, &reqList[sendcnt]);
            sendcnt++;
          }
        }
      if(id_thread==3)
        if(communicationCnt[3] == 1)
        {
          #pragma omp critical
          {
            MPI_Isend(&vector[1][ny-2-offset][1], 1, xzFacetype,
            right_neighborY,tag_YR, comm, &reqList[sendcnt]);
            sendcnt++;
          }
        }
      if(id_thread==0)
        if(communicationCnt[4] == 1)
        {
          #pragma omp critical
          {
            MPI_Isend(&vector[1][1][1+offset],1, xyFacetype,
            left_neighborZ, tag_ZL, comm, &reqList[sendcnt]);
            sendcnt++;
          }
        }
      if(id_thread==1)
        if(communicationCnt[5] == 1)
        {
          #pragma omp critical
          {
            MPI_Isend(&vector[1][1][nz-2-offset], 1, xyFacetype,
            right_neighborZ,tag_ZR, comm, &reqList[sendcnt]);
            sendcnt++;
          }
        }
      #pragma omp barrier
      if(id_thread==0)
      {
        assert_eq(recvcnt,sendcnt-recvcnt);
      }
      #pragma omp barrier
      …
    }
    …
  }

The phase of communication between edges also consists of the number of calls to the MPI_Irecv() functions and to the MPI_Isend() functions. The phase of communication between corners uses two calls to MPI_Irecv() and two calls to MPI_Isend(). In order to parallelise these two phases, exactly the same approach has been used as for the phase of communication between faces.

Scaling tests
^^^^^^^^^^^^^

In this section we will discuss the performance results from the weak scaling tests. Tests were performed on the Beskow supercomputer at the PDC Center for High Performance Computing at the KTH Royal Institute of Technology. To compare the original version of the iPIC3D code with the new version, based on MPI threading support, we used two standard simulation cases called GEM 3D and Magnetosphere 3D. In addition, we used two different data sizes/regimes for both simulation cases:

* Field solver dominated regime. Here, a relatively small number of particles (27 per cell) is used. The most computationally expensive part of the iPIC3D code results in the Maxwell field solver.
* Particle dominated regime. It is characterised by a large number of particles (1,000 per cell). The most computationally expensive part of the iPIC3D code results in the particle mover.

Thus, there were four different cases (two simulation versions with two regimes), and the input data for each case for 8 MPI processes are stated below.

The input data for the field solver dominated GEM 3D simulation: ::

   dt = 0.15 # dt = time step
   ncycles = 30000    # number of cycles
   Lx = 10 # Lx = simulation box length - X direction
   Ly = 10 # Ly = simulation box length - Y direction
   Lz = 10 # Lz = simulation box length - Z direction
   nxc = 20  # nxc = number of cells - X direction
   nyc = 20  # nyc = number of cells - Y direction
   nzc = 20  # nzc = number of cells - Z direction
   npcelx = 3 3 3 3  # npcelx = number of particles per cell – X direction
   npcely = 3 3 3 3  # npcely = number of particles per cell – Y direction
   npcelz = 3 3 3 3  # npcelz = number of particles per cell – Z direction
   CGtol = 1E-3	# CG solver stopping criterion tolerance
   GMREStol = 1E-3   # GMRES solver stopping criterion tolerance

The input data for the particle dominated GEM 3D simulation: ::

   dt = 0.15 # dt = time step
   ncycles = 1800     # number of cycles
   Lx = 10 # Lx = simulation box length - X direction
   Ly = 10 # Ly = simulation box length - Y direction
   Lz = 10 # Lz = simulation box length - Z direction
   nxc = 20  # nxc = number of cells - X direction
   nyc = 20  # nyc = number of cells - Y direction
   nzc = 20  # nzc = number of cells - Z direction
   npcelx = 10 10 10 10	    # npcelx = number of particles per cell – X direction
   npcely = 10 10 10 10	    # npcely = number of particles per cell – Y direction
   npcelz = 10 10 10 10	    # npcelz = number of particles per cell – Z direction
   CGtol = 1E-3	  # CG solver stopping criterion tolerance
   GMREStol = 1E-3  # GMRES solver stopping criterion tolerance

The input data for the field solver dominated Magnetosphere 3D simulation: ::

    dt = 0.15  # dt = time step
    ncycles = 20000 # number of cycles
    Lx = 10 # Lx = simulation box length - X direction
    Ly = 10 # Ly = simulation box length - Y direction
    Lz = 10 # Lz = simulation box length - Z direction
    nxc = 30  # nxc = number of cells - X direction
    nyc = 30  # nyc = number of cells - Y direction
    nzc = 30  # nzc = number of cells - Z direction
    npcelx = 3 3    # npcelx = number of particles per cell – X direction
    npcely = 3 3    # npcely = number of particles per cell – Y direction
    npcelz = 3 3    # npcelz = number of particles per cell – Z direction
    CGtol = 1E-3    # CG solver stopping criterion tolerance
    GMREStol = 1E-3 # GMRES solver stopping criterion tolerance

The input data for the particle dominated Magnetosphere 3D simulation: ::

    dt = 0.15  # dt = time step
    ncycles = 800   # number of cycles
    Lx = 10 # Lx = simulation box length - X direction
    Ly = 10 # Ly = simulation box length - Y direction
    Lz = 10 # Lz = simulation box length - Z direction
    nxc = 30  # nxc = number of cells - X direction
    nyc = 30  # nyc = number of cells - Y direction
    nzc = 30  # nzc = number of cells - Z direction
    npcelx = 10 10  # npcelx = number of particles per cell – X direction
    npcely = 10 10  # npcely = number of particles per cell – Y direction
    npcelz = 10 10  # npcelz = number of particles per cell – Z direction
    CGtol = 1E-3    # CG solver stopping criterion tolerance
    GMREStol = 1E-3 # GMRES solver stopping criterion tolerance

Each of four test cases (two version of codes for two simulation cases) has been executed on the increasing number of cores from 32, 64, 128 and up to 256 cores. The presented results here are for the cases with two, four and eight OpenMP threads. In order to ensure a fair comparison, the number of iterations in the linear solver was fixed to 20, although in a real simulation the number of iterations depends on the speed of convergence.

Figures 2, 3, 4, and 5 show the results of the weak scaling tests for each of the four cases. Three-dimensional decomposition of MPI processes on X-, Y- and Z-axes was used, resulting in different topologies of MPI processes, each having two, four, and then eight threads in addition. The total number of particles and cells in a simulation are calculated from nxc*nyc*nzc*npcelx*npcely*npcelz and nxc*nyc*nzc, respectively. Thus, for example, for the particle dominated Magnetosphere 3D simulation on 32 cores (2x2x2 MPI processes x 4 OpenMP threads), there were used 27x106 particles and 30x30x30 cells, and the simulation size increased proportionally to the number of processes.

.. figure:: gemf.png

   Figure 2. Weak scaling test for the field solver dominated GEM 3D simulation of the original and new hybrid versions of the iPIC3D code

.. figure:: gemp.png

   Figure 3. Weak scaling test for the particle dominated GEM 3D simulation of the original and new hybrid versions of the iPIC3D code

.. figure:: magf.png

   Figure 4. Weak scaling test for the field solver dominated Magnetosphere 3D simulation of the original and new hybrid versions of the iPIC3D code

.. figure:: magp.png

   Figure 5. Weak scaling test for the particle dominated Magnetosphere 3D simulation of the original and new hybrid versions of the iPIC3D code

These figures show that, for all four simulations, in most of the cases the original iPIC3D code with disabled OpenMP (meaning pure MPI) was slowest. In the field solver dominated regimes (the two left pictures in Figure 3), the new version of iPIC3D is slower than the original version with enabled OpenMP on smaller number of cores, while on 256 cores the new hybrid implementation shows the lowest execution time (for example, with eight threads per one MPI process, it is 18-20% faster). In the particle dominated regimes (the two right pictures in Figure 3), the original and new iPIC3D implementations with equal number of threads show almost identical performance. Also, for all test runs of the new version of iPIC3D with two, four and eight threads, the fastest was the new version with eight threads per one MPI process, thus showing the scalability. On 256 cores it shows 16-31% of speedup comparing to iPIC3D based on the pure MPI and 9-14% when comparing to the new version with two threads per one MPI process.