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

- **Multiple Polarimetric Operators:**
  - **Zhang21**: Lookup table based (Zhang et al. 2021)
  - **TCWA2**: Analytical gamma distribution formulation
  
- **Microphysics Schemes:**
  - Thompson (double-moment for rain)
  - WSM6 (single-moment)
  - NSSL (double-moment with hail)
  
- **Radar Frequencies:**
  - S-band and C-band support
  
- **Advanced Physics:**
  - Melting layer treatment (Zhang21)
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
│   ├── dualpol_op_mod.f90       # Main interface (operator dispatch)
│   ├── zhang21/                 # Zhang21 operator
│   │   ├── zhang21_core_mod.f90   # Zhang21 core physics
│   │   ├── zhang21_tlad_mod.f90   # Zhang21 TL/AD
│   │   └── CMakeLists.txt
│   └── tcwa2/                   # TCWA2 operator
│       ├── tcwa2_core_mod.f90     # TCWA2 core physics
│       └── CMakeLists.txt
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

**Using Zhang21 operator (requires coefficient files):**
```fortran
use dualpol_op_mod, only: ppro_init_coefs, ppro_compute_point

! Initialize with Zhang21 operator
call ppro_init_coefs(operator_name='Zhang21')

! Compute dual-pol variables
call ppro_compute_point(iband, 'THOMPSON', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, nr=nr)
```

**Using TCWA2 operator (analytical formulation, no coefficient files needed):**
```fortran
use dualpol_op_mod, only: ppro_init_coefs, ppro_compute_point

! Initialize with TCWA2 operator (no coefficient files required)
call ppro_init_coefs(operator_name='TCWA2')

! Compute dual-pol variables (TCWA2 requires additional parameters)
! Note: TCWA2 needs cloud ice (qi, ni), cloud water (qc), and melting fractions (smlf, gmlf)
call ppro_compute_point(iband, 'THOMPSON', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, &
                        nr=nr, ns=ns, ng=ng, ni=ni, qi=qi, qc=qc, &
                        smlf=smlf, gmlf=gmlf)
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

**Zhang21 operator** requires coefficient files for S-band and C-band radars. These files should be placed in the `fix/` directory. See `fix/README.txt` for details on the required files.

**TCWA2 operator** does not require coefficient files as it uses analytical formulations.

**Note**: Coefficient files for Zhang21 must be obtained separately and are not included in this repository.

## Module Structure

### Main Interface (`dualpol_op_mod.f90`)

Unified interface providing:
- `ppro_init_coefs()`: Initialize library with operator selection
- `ppro_compute_point()`: Main computation (dispatches to selected operator)
- `ppro_set_operator()`: Switch between operators
- `ppro_finalize()`: Cleanup

### Zhang21 Operator (`zhang21/`)

**zhang21_core_mod.f90** - Lookup table based operator:
- `zhang21_init_coefs()`: Load coefficient lookup tables
- `zhang21_compute_point()`: Compute using lookup tables
- `melting_scheme_zhang24()`: Melting layer treatment
- `dualpol_op_rain/icephase/total()`: Individual hydrometeor calculations

**zhang21_tlad_mod.f90** - Tangent-linear and adjoint for DA applications

### TCWA2 Operator (`tcwa2/`)

**tcwa2_core_mod.f90** - Analytical gamma distribution based:
- `dualpol_op_rain_tcwa2()`: Rain using analytical formulas
- `dualpol_op_ice_tcwa2()`: Ice using analytical formulas
- `dualpol_op_snow_tcwa2()`: Snow with melting
- `dualpol_op_graup_tcwa2()`: Graupel with melting
- `GAMLN()`: Gamma function computation

## API Example

**Zhang21 operator:**
```fortran
use dualpol_op_mod

! Initialize Zhang21 operator
call ppro_init_coefs(operator_name='Zhang21')

! Compute dual-pol variables
call ppro_compute_point(iband, 'THOMPSON', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, nr=nr)
```

**TCWA2 operator:**
```fortran
use dualpol_op_mod

! Initialize TCWA2 operator
call ppro_init_coefs(operator_name='TCWA2')

! Compute dual-pol variables
call ppro_compute_point(iband, 'THOMPSON', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, &
                        nr=nr, ns=ns, ng=ng, ni=ni, qi=qi, qc=qc, &
                        smlf=smlf, gmlf=gmlf)
```

## Version History

- **v1.0.0** (2025): Initial modularized release
  - Extracted from JEDI-UFO as independent library
  - Multi-operator architecture (Zhang21 + TCWA2)
  - Modular directory structure (zhang21/ and tcwa2/ subdirectories)
  - Support for Thompson, WSM6, and NSSL microphysics
  - S-band and C-band radar support
  - Standalone examples and comprehensive documentation

## Authors

- **Zhiquan (Jake) Liu** (NCAR/MMM) - Zhang21 operator implementation
- **Tzu-Chin Tsai** - TCWA2 operator implementation
- **Hejun Xie** - Integration into JEDI-UFO framework, TL/AD development
- **Tao Sun** - Adding hail categories and extended microphysics support
- **Rong Kong** (NCAR/MMM) - Bug fixes, operator tuning and testing, modularization and multi-operator architecture

## License

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

## References

1. Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

2. Liu, et al., 2024: (Melting scheme reference): A New Melting Model and Its Implementation in Parameterized Forward Operators for Polarimetric Radar Data Simulation With Double Moment Microphysics Schemes, JGR, 129.

## Contact

For questions, issues, or contributions:
- **Technical Lead**: Rong Kong (rkong@ucar.edu)
- **Issues**: Please open an issue in the repository

## Acknowledgments

This work was developed at the National Center for Atmospheric Research (NCAR), Mesoscale and Microscale Meteorology Laboratory (MMM).

PPRO is designed as a general-purpose, standalone dual-polarization radar operator library. While it is compatible with data assimilation frameworks such as JEDI, it can be used independently in any numerical weather prediction system or forward operator application.

