# PPRO: Parameterized Polarimetric Radar Operator Library

## Overview

PPRO is an independent, standalone library for computing dual-polarization radar variables from atmospheric model state variables. The library provides a modular, reusable implementation of the Parameterized Polarimetric Radar Operator (P-PRO) that can be used in any numerical weather prediction or data assimilation system.

This library can be used independently or integrated into data assimilation frameworks such as JEDI-UFO.

## Scientific Background

The P-PRO operator is based on the parameterized forward operators for polarimetric radar data described in:

> Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

The operator computes dual-polarization radar variables including:
- **Zhh**: Horizontal reflectivity (dBZ)
- **ZDR**: Differential reflectivity (dB)
- **KDP**: Specific differential phase (deg/km)
- **ρhv**: Copolar correlation coefficient

## Features

- Support for multiple microphysics schemes:
  - Thompson (double-moment for rain)
  - WSM6 (single-moment)
  - NSSL (double-moment with hail)
- S-band and C-band radar frequencies
- Melting layer treatment following Zhang et al. (2024)
- Pure and melting hydrometeor categories (rain, snow, graupel, hail)

## Directory Structure

```
ppro/
├── CMakeLists.txt          # Main CMake build configuration
├── VERSION.cmake           # Version information
├── ppro-config.cmake.in    # CMake package config template
├── README.md              # This file
├── LICENSE                # License information
├── cmake/                 # CMake modules
│   ├── ppro_compiler_flags.cmake
│   ├── compiler_flags_GNU_Fortran.cmake
│   ├── compiler_flags_Intel_Fortran.cmake
│   └── compiler_flags_IntelLLVM_Fortran.cmake
├── src/                   # Source code
│   ├── CMakeLists.txt
│   ├── dualpol_op_mod.f90       # Core dual-pol operator module
│   └── dualpol_op_tlad_mod.f90  # Tangent-linear and adjoint module
├── examples/              # Standalone example programs
│   ├── CMakeLists.txt
│   ├── standalone_test.f90      # Comprehensive test/demo program
│   └── README.md                # Examples documentation
└── fix/                   # Coefficient data files
    ├── README.txt
    ├── sband_rain_coefs.txt
    ├── sband_snow_coefs.txt
    ├── sband_graupel_coefs.txt
    ├── sband_hail_coefs.txt
    ├── cband_rain_coefs.txt
    ├── cband_snow_coefs.txt
    ├── cband_graupel_coefs.txt
    └── cband_hail_coefs.txt
```

## Building

### Prerequisites

- CMake 3.20 or higher
- Fortran compiler (GNU, Intel, or IntelLLVM)
- Optional: OpenMP for parallel computation

### Build Instructions

```bash
# Clone or copy the ppro directory
cd ppro

# Create build directory
mkdir build
cd build

# Configure with CMake
cmake ..

# Build the library
make

# Install (optional)
make install
```

### CMake Options

- `OPENMP`: Enable OpenMP support (default: ON)
- `BUILD_SHARED_LIBS`: Build shared libraries instead of static (default: OFF)
- `BUILD_EXAMPLES`: Build standalone example programs (default: ON)
- `CMAKE_INSTALL_PREFIX`: Installation prefix

Example with custom options:
```bash
cmake -DOPENMP=OFF -DBUILD_EXAMPLES=OFF -DCMAKE_INSTALL_PREFIX=/path/to/install ..
```

### Running Examples

After building, run the standalone test:
```bash
cd build/bin
./ppro_standalone_test
```

See `examples/README.md` for more details on examples and how to write your own programs using PPRO.

## Usage

### Standalone Usage

The library can be used independently in any Fortran-based numerical weather prediction model:

```fortran
use dualpol_op_mod, only: ppro_init_coefs, ppro_compute_point

! Initialize coefficients (call once)
call ppro_init_coefs()

! Compute dual-pol variables for a grid point
call ppro_compute_point(iband, 'THOMPSON', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, &
                        nr=nr)
```

### Integration with JEDI-UFO

The library can also be integrated into JEDI-UFO for data assimilation applications:

#### Option 1: Using find_package (Recommended for production)

After installing ppro, configure UFO with:
```cmake
find_package(PPRO REQUIRED)
target_link_libraries(ufo PUBLIC ppro)
```

#### Option 2: Using add_subdirectory (Recommended for development)

Place ppro in your JEDI bundle and use:
```cmake
ecbuild_bundle(PROJECT ppro SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/ppro)
target_link_libraries(ufo PUBLIC ppro)
```

## Coefficient Files

The library requires coefficient files for S-band and C-band radars. These files should be placed in the `fix/` directory. See `fix/README.txt` for details on the required files.

**Note**: Coefficient files must be obtained separately and are not included in this repository.

## Module Structure

### dualpol_op_mod.f90

Core physics module containing:
- `init_coefs()`: Initialize lookup table coefficients
- `dualpol_op_rain()`: Rain dual-pol variables
- `dualpol_op_icephase()`: Ice-phase dual-pol variables
- `dualpol_op_total()`: Aggregate total dual-pol variables
- `melting_scheme_zhang24()`: Melting layer partitioning
- `dm_z_2moment()`: Mean diameter and reflectivity for 2-moment schemes
- `dm_z_wsm6()`: Mean diameter and reflectivity for WSM6

### dualpol_op_tlad_mod.f90

Tangent-linear and adjoint module for data assimilation applications.

## API Example

```fortran
use dualpol_op_mod

! Initialize coefficients
call init_coefs()

! Compute dual-pol variables for a grid point
call ufo_PPRO_sim1obs(iband, density_air, temp_air, &
                      qr, qs, qg, zh, zdr, kdp, phv, &
                      nr=nr)  ! Optional: for double-moment schemes
```

## Version History

- **v1.0.0** (2025): Initial modularized release
  - Extracted from JEDI-UFO
  - Support for Thompson, WSM6, and NSSL microphysics
  - S-band and C-band radar support

## Authors

- **Zhiquan (Jake) Liu** (NCAR/MMM) - Original P-PRO implementation
- **Hejun Xie** - Integration of P-PRO operator into JEDI-UFO framework
- **Tao Sun** - Adding hail categories and extended microphysics support
- **Rong Kong** (NCAR/MMM) - Bug fixes, operator tuning and testing, modularization and external library development

## License

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

## References

1. Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

2. Zhang, G., et al., 2024: (Melting scheme reference)

## Contact

For questions, issues, or contributions:
- **Technical Lead**: Rong Kong (rkong@ucar.edu)
- **Issues**: Please open an issue in the repository

## Acknowledgments

This work was developed at the National Center for Atmospheric Research (NCAR), Mesoscale and Microscale Meteorology Laboratory (MMM).

PPRO is designed as a general-purpose, standalone dual-polarization radar operator library. While it is compatible with data assimilation frameworks such as JEDI, it can be used independently in any numerical weather prediction system or forward operator application.

