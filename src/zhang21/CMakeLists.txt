# CMakeLists.txt for Zhang21 operator modules

# Zhang21 core physics module (independent)
add_library(${PROJECT_NAME}_zhang21_core OBJECT zhang21_core_mod.f90)
set_target_properties(${PROJECT_NAME}_zhang21_core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})

# Zhang21 TL/AD module (depends on zhang21_core)
add_library(${PROJECT_NAME}_zhang21_tlad OBJECT zhang21_tlad_mod.f90)
set_target_properties(${PROJECT_NAME}_zhang21_tlad PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
target_include_directories(${PROJECT_NAME}_zhang21_tlad PRIVATE
    ${CMAKE_BINARY_DIR}/${MODULE_DIR})
add_dependencies(${PROJECT_NAME}_zhang21_tlad ${PROJECT_NAME}_zhang21_core)

