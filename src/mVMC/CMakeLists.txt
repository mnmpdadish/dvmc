# include guard
cmake_minimum_required(VERSION 2.8.0 )

add_definitions(-D_mVMC)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
add_definitions(-DMEXP=19937)

set(SOURCES_vmcmain
        vmcmain.c physcal_lanczos.c splitloop.c 
 )

set(SOURCES_sfmt
        ../sfmt/SFMT.c   
 )

add_library(StdFace STATIC ../StdFace/ChainLattice.c ../StdFace/HoneycombLattice.c ../StdFace/SquareLattice.c ../StdFace/StdFace_main.c ../StdFace/StdFace_ModelUtil.c ../StdFace/TriangularLattice.c ../StdFace/Ladder.c ../StdFace/Kagome.c ../StdFace/Orthorhombic.c ../StdFace/Pyrochlore.c ../StdFace/Wannier90.c ../StdFace/FCOrtho.c)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../pfapack)

add_executable(dvmcdry.out vmcdry.c)
target_link_libraries(dvmcdry.out StdFace m)
add_executable(dvmc.out ${SOURCES_vmcmain} ${SOURCES_sfmt})
target_link_libraries(dvmc.out StdFace)
target_link_libraries(dvmc.out pfapack StdFace ${LAPACK_LIBRARIES} m)

if(USE_SCALAPACK)
  string(REGEX REPLACE "-L[ ]+" "-L" sc_libs ${SCALAPACK_LIBRARIES})
  string(REGEX REPLACE "[ ]+" ";" sc_libs ${sc_libs})
  foreach(sc_lib IN LISTS sc_libs)
    target_link_libraries(dvmc.out ${sc_lib})
  endforeach(sc_lib)
  message(STATUS "SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")
endif(USE_SCALAPACK)

if(MPI_FOUND)
  target_link_libraries(dvmc.out ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
install(TARGETS dvmcdry.out RUNTIME DESTINATION bin)
install(TARGETS dvmc.out RUNTIME DESTINATION bin)
add_definitions(-D_mVMC)
