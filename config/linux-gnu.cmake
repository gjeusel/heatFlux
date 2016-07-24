set(CMAKE_BUILD_TYPE Debug CACHE STRING "")

# setting compilers
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")
set(CMAKE_C_COMPILER       gcc   CACHE STRING "")
set(CMAKE_CXX_COMPILER     g++  CACHE STRING "")

set(GFORTRAN_WARNINGS "-Warray-bounds -Wcharacter-truncation -Wconversion -Wimplicit-interface  -Wimplicit-procedure  -Wline-truncation -Wintrinsics-std  -Wsurprising  -Wno-tabs  -Wunderflow -Wunused-parameter")
set(CMAKE_Fortran_FLAGS_DEBUG "${GFORTRAN_WARNINGS} -g -pipe -fno-omit-frame-pointer  -fbounds-check -fbacktrace" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-pipe -O3 -ftree-vectorize" CACHE STRING "")

set(MPIRUN_CMD "mpiexec -n $NP" CACHE STRING "ex: mpiexec, mpirun.Ompi, mpiexec_mpt, mpirun -np $NP...")


set(CMAKE_Fortran_FLAGS_DEBWITHGPROF "-pg -g ${CMAKE_Fortran_FLAGS_DEBUG}" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELWITHGPROF "-pg -g ${CMAKE_Fortran_FLAGS_RELEASE}" CACHE STRING "")

set(CMAKE_BUILD_TYPE_STRINGS_EXTRA "DebWithGPROF" "RelWithGPROF")
