# CMakeLists.txt for TCWA2 operator module

# TCWA2 forward operator module (independent, no dependencies on Zhang21)
add_library(${PROJECT_NAME}_tcwa2 OBJECT tcwa2_forward_mod.f90)

set_target_properties(${PROJECT_NAME}_tcwa2 PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})

# Make modules available
target_include_directories(${PROJECT_NAME}_tcwa2 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>)

