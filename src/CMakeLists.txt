# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

include(GNUInstallDirs)

# Option to select between shared and static libraries
option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones" OFF)

# Set library type (static or shared)
set(LIBRARY_TYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(LIBRARY_TYPE SHARED)
endif()

# Set module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})

# Build the library in two steps to ensure proper module dependency
# Step 1: Compile dualpol_op_mod.f90 first (base module)
add_library(${PROJECT_NAME}_base OBJECT dualpol_op_mod.f90)
set_target_properties(${PROJECT_NAME}_base PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})

# Step 2: Compile dualpol_op_tlad_mod.f90 (depends on base)
add_library(${PROJECT_NAME}_tlad OBJECT dualpol_op_tlad_mod.f90)
set_target_properties(${PROJECT_NAME}_tlad PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
target_include_directories(${PROJECT_NAME}_tlad PRIVATE ${CMAKE_BINARY_DIR}/${MODULE_DIR})
add_dependencies(${PROJECT_NAME}_tlad ${PROJECT_NAME}_base)

# Step 3: Combine into final library
add_library(${PROJECT_NAME} ${LIBRARY_TYPE} 
    $<TARGET_OBJECTS:${PROJECT_NAME}_base>
    $<TARGET_OBJECTS:${PROJECT_NAME}_tlad>
)

# Link dependencies (if needed, currently none)
if(OPENMP AND OpenMP_Fortran_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Set the Fortran module directory for build and install phases
set_target_properties(${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})

# Install Fortran modules into the correct destination
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${PPRO_INSTALL_PREFIX}/${MODULE_DIR})

# Set include directories for build and install separately
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>  # During the build
    $<INSTALL_INTERFACE:${MODULE_DIR}>)                   # After install

