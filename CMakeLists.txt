cmake_minimum_required(VERSION 3.20)

project(2S-ESS 
  VERSION 2.4.2
  LANGUAGES Fortran)

option (ENABLE_SHARED            "Enable shared libraries" OFF)

enable_testing()


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
include_directories(${CMAKE_BINARY_DIR}/mod_files)
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
message(STATUS "Bin Dir " ${CMAKE_BINARY_DIR})
#=====================================================
# Compiler and Options
#=====================================================
# ... write stuff here
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-fpp -D IFORT -O3 -heap-arrays")
  set(CMAKE_Fortran_FLAGS_DEBUG "-fpp -D IFORT -g -warn all -check all -zero")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-std=legacy -mcmodel=large -O3 -ffree-form -ffree-line-length-none -fall-intrinsics -fimplicit-none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-std=legacy -mcmodel=large -g -O0 -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -Wall -fcheck=all -fbacktrace -fall-intrinsics -Werror=underflow")
endif()


#=====================================================
# Planet configuration modules
#=====================================================
set(PLANET_CONFIG_SOURCES
    sourcecode/planet_config_mod.f90
    sourcecode/planet_config_io.f90
)

#=====================================================
# twostr
#=====================================================
file(GLOB REG_2S_SRC sourcecode/2stream_2p4_feb15/*.f90)

set(2S_PATH sourcecode/2stream_2p4_optimized)
set(2S_SOL_PATH sourcecode/2stream_2p4_optimized/solar)
set(2S_OBS_PATH sourcecode/2stream_2p4_optimized/solar/obs)
set(FO_PATH_1P4 sourcecode/FO_1p4_optimized)

set(SOURCES_2S 
   ${2S_PATH}/2stream_Taylor_optimized.f90
   ${2S_PATH}/2stream_brdf_kernels.f90
   ${2S_OBS_PATH}/2stream_brdf_supplement_solar_obs.f90
   ${2S_PATH}/2stream_sleave_routines.f90
   ${2S_PATH}/2stream_sleave_supplement_optimized.f90
   ${2S_PATH}/2stream_geometry_optimized.f90
   ${2S_OBS_PATH}/2stream_geometry_solar_obs.f90
   ${2S_OBS_PATH}/2stream_inputs_solar_obs.f90
   ${2S_OBS_PATH}/2stream_writemodules_solar_obs.f90
   ${2S_PATH}/2stream_miscsetups_optimized.f90
   ${2S_OBS_PATH}/2stream_miscsetups_solar_obs.f90
   ${2S_OBS_PATH}/2stream_miscsetups_solar_obs_pp.f90
   ${2S_PATH}/2stream_solutions_optimized.f90 
   ${2S_SOL_PATH}/2stream_solutions_solar.f90 
   ${2S_PATH}/2stream_bvproblem_optimized.f90 
   ${2S_OBS_PATH}/2stream_intensity_solar_obs.f90
   ${2S_OBS_PATH}/2stream_converge_obs.f90
   ${2S_PATH}/2stream_fluxes.f90
)
#  masters
set(SOURCES_2S_MASTER
   ${SOURCES_2S}     
   ${2S_OBS_PATH}/2stream_master_solar_obs.f90
)

#=====================================================
# ESS
#=====================================================
file(GLOB REG_FO_SRC sourcecode/FO_1p4_sourcecode/*.f90)
file(GLOB OPT_FO_SRC sourcecode/FO_1p4_optimized/*.f90)

#=====================================================
# twostr library (static)
#=====================================================
add_library(TWOSTR STATIC ${REG_2S_SRC} ${PLANET_CONFIG_SOURCES})
add_library(TWOSTR_OPT STATIC ${SOURCES_2S_MASTER} ${PLANET_CONFIG_SOURCES})

install(
  TARGETS TWOSTR
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )

install(
  TARGETS TWOSTR_OPT
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )

#=====================================================
# ESS library (static)
#=====================================================
add_library(ESS STATIC ${REG_FO_SRC})
add_library(ESS_OPT STATIC ${OPT_FO_SRC})
install(
  TARGETS ESS
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )
install(
  TARGETS ESS_OPT
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
  )


#====================================================
# Add executables
#====================================================
set(EXE regtest)
add_executable(${EXE} TEST/test_regular.f90)
target_link_libraries(${EXE} TWOSTR ESS)
install(TARGETS ${EXE} RUNTIME DESTINATION bin)

set(EXE_OPT opttest)
add_executable(${EXE_OPT} TEST/test_opt.f90)
target_link_libraries(${EXE_OPT} TWOSTR_OPT ESS_OPT)
install(TARGETS ${EXE_OPT} RUNTIME DESTINATION bin)

# Planetary configuration test
set(EXE_PLANET planettest)
add_executable(${EXE_PLANET} TEST/test_planetary_config.f90)
target_link_libraries(${EXE_PLANET} TWOSTR ESS)
install(TARGETS ${EXE_PLANET} RUNTIME DESTINATION bin)

# Earth namelist test
set(EXE_EARTH earthtest)
add_executable(${EXE_EARTH} TEST/test_earth_namelist.f90)
target_link_libraries(${EXE_EARTH} TWOSTR ESS)
install(TARGETS ${EXE_EARTH} RUNTIME DESTINATION bin)

# Mars namelist test
set(EXE_MARS marstest)
add_executable(${EXE_MARS} TEST/test_mars_namelist.f90)
target_link_libraries(${EXE_MARS} TWOSTR ESS)
install(TARGETS ${EXE_MARS} RUNTIME DESTINATION bin)

# Earth namelist optimized test
set(EXE_EARTH_OPT earthtest_opt)
add_executable(${EXE_EARTH_OPT} TEST/test_earth_namelist_opt.f90)
target_link_libraries(${EXE_EARTH_OPT} TWOSTR_OPT ESS_OPT)
install(TARGETS ${EXE_EARTH_OPT} RUNTIME DESTINATION bin)


#====================================================
# Add test
#====================================================
add_test(
  NAME "regular run"
  COMMAND $<TARGET_FILE:${EXE}>
)
set_tests_properties("regular run" PROPERTIES LABELS "run;reg")

add_test(
  NAME "optimized run"
  COMMAND $<TARGET_FILE:${EXE_OPT}>
)
set_tests_properties("optimized run" PROPERTIES LABELS "run;opt")

add_test(
  NAME "regular output"
  COMMAND numdiff -S -r 1e-14 ${CMAKE_SOURCE_DIR}/TEST/expected_test_results.dat test_regular_output.dat 
)
set_tests_properties("regular output" PROPERTIES LABELS "test;reg")

add_test(
  NAME "optimized output"
  COMMAND numdiff -S -r 1e-14 ${CMAKE_SOURCE_DIR}/TEST/expected_test_results.dat test_opt_output.dat 
)
set_tests_properties("optimized output" PROPERTIES LABELS "test;opt")

add_test(
  NAME "planetary config"
  COMMAND $<TARGET_FILE:${EXE_PLANET}>
)
set_tests_properties("planetary config" PROPERTIES LABELS "test;planet")

add_test(
  NAME "earth namelist"
  COMMAND $<TARGET_FILE:${EXE_EARTH}>
)
set_tests_properties("earth namelist" PROPERTIES LABELS "test;earth")

add_test(
  NAME "mars namelist"
  COMMAND $<TARGET_FILE:${EXE_MARS}>
)
set_tests_properties("mars namelist" PROPERTIES LABELS "test;mars")

add_test(
  NAME "earth output"
  COMMAND numdiff -S -r 1e-12 ${CMAKE_SOURCE_DIR}/TEST/expected_test_results.dat test_earth_namelist_output.dat 
)
set_tests_properties("earth output" PROPERTIES LABELS "test;earth;compare")

add_test(
  NAME "earth namelist optimized"
  COMMAND $<TARGET_FILE:${EXE_EARTH_OPT}>
)
set_tests_properties("earth namelist optimized" PROPERTIES LABELS "test;earth;opt")

add_test(
  NAME "earth optimized output"
  COMMAND numdiff -S -r 1e-12 ${CMAKE_SOURCE_DIR}/TEST/expected_test_results.dat test_earth_namelist_opt_output.dat 
)
set_tests_properties("earth optimized output" PROPERTIES LABELS "test;earth;opt;compare")
