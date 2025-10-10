# CMakeLists.txt for PPRO examples

# Standalone test executable
add_executable(ppro_standalone_test standalone_test.f90)

# Link with ppro library
target_link_libraries(ppro_standalone_test PRIVATE ppro)

# Set Fortran module directory
set_target_properties(ppro_standalone_test PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/examples
)

# Installation (optional)
install(TARGETS ppro_standalone_test
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT examples
)

# Add to examples component
message(STATUS "PPRO standalone test example configured")

