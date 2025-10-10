cmake_minimum_required(VERSION 3.20)

project(ppro VERSION 1.0.0 LANGUAGES Fortran)

option(OPENMP "Build ppro with OpenMP support" ON)
option(FIX_FILE_PATH "Path to coefficient files (default: fix/)" OFF)
option(BUILD_EXAMPLES "Build standalone example programs" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})

# Include GNUInstallDirs to get the standard installation directory variables
include(GNUInstallDirs)

# Set custom installation paths for ppro only
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(PPRO_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/lib" CACHE PATH "ppro install path" FORCE)
else()
    set(PPRO_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

message("PPRO_INSTALL_PREFIX set to ${PPRO_INSTALL_PREFIX}")

# Set output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Function to get the current git commit hash
function(get_git_hash result)
    execute_process(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(${result} ${GIT_HASH} PARENT_SCOPE)
endfunction()

# Get the git hash
get_git_hash(PROJECT_HASH)

# Define ANSI color code for green
string(ASCII 27 ESCAPE)
set(RESET "${ESCAPE}[0m")
set(GREEN "${ESCAPE}[32m")

# Output the project name, version, and hash in green text
message(STATUS "${GREEN}[${PROJECT_NAME}] (${PROJECT_VERSION}) [${PROJECT_HASH}]${RESET}")

# Output the various install/output paths
message("CMAKE_INSTALL_PREFIX set to ${CMAKE_INSTALL_PREFIX}")
message("CMAKE_INSTALL_LIBDIR set to ${CMAKE_INSTALL_LIBDIR}")
message("CMAKE_LIBRARY_OUTPUT_DIRECTORY set to ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

## Configuration options
include(${PROJECT_NAME}_compiler_flags)
include(GNUInstallDirs)

## Dependencies
if(OPENMP)
  find_package(OpenMP COMPONENTS Fortran)

	# Check if OMP_NUM_THREADS is set in the environment
	if(DEFINED ENV{OMP_NUM_THREADS})
		set(OMP_NUM_THREADS $ENV{OMP_NUM_THREADS})
	else()
		# Set a default value if not set
		set(OMP_NUM_THREADS "1")
	endif()
	
	# Export OMP_NUM_THREADS to the environment
	set(ENV{OMP_NUM_THREADS} ${OMP_NUM_THREADS})
endif()

## Sources
add_subdirectory(src)

## Examples
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Generate and configure the project config file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/ppro-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ppro-config.cmake"
  @ONLY
)

# Installation for ppro
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}-config
        ARCHIVE DESTINATION ${PPRO_INSTALL_PREFIX}/lib
        LIBRARY DESTINATION ${PPRO_INSTALL_PREFIX}/lib)

# Generate the version file for the config package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ppro-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Install the configuration files
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/ppro-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/ppro-config-version.cmake"
  DESTINATION ${PPRO_INSTALL_PREFIX}/lib/cmake/ppro
)

# Install the export set for use with the install-tree
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${PPRO_INSTALL_PREFIX}/lib/cmake/ppro
)

# Install mod files
install(DIRECTORY ${CMAKE_BINARY_DIR}/module/
  DESTINATION ${PPRO_INSTALL_PREFIX}/include/${PROJECT_NAME}
)

