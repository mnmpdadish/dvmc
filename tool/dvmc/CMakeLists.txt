# include guard
cmake_minimum_required(VERSION 2.8.0 )
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_executable(dvmc_gnuplot src/dvmc_gnuplot.c)
add_library(dvmc_speedup SHARED src/dvmc_speedup.c)
set_property(TARGET dvmc_speedup PROPERTY POSITION_INDEPENDENT_CODE ON)

install(TARGETS dvmc_gnuplot RUNTIME DESTINATION bin)
install(TARGETS dvmc_speedup DESTINATION bin)

set(SCRIPTS averageOhgoe.py cluster.py dos_analysis.py dvmc_spectrum.py generate_template_gnuplot.py init_groundstate makeDH4.py makeExcitation.py mergeOutputBin.py postprocess_dVMC.sh preprocess_dVMC.sh selectExcitation.py select_kPath.py template_plot_Akw.gp template_plot_allAkw.gp template_plot_dos.gp)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPTS} DESTINATION .)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPTS} DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

