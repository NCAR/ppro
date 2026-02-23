# PPRO: Parameterized Polarimetric Radar Operator Library

## Overview

PPRO is an independent, standalone library for computing dual-polarization radar variables from atmospheric model state variables. The library provides a modular, reusable implementation of the Parameterized Polarimetric Radar Operator (P-PRO) that can be used in any numerical weather prediction or data assimilation system.

This library can be used independently or integrated into data assimilation frameworks such as JEDI-UFO.

## Scientific Background

The P-PRO operator is based on the parameterized forward operators for polarimetric radar data described in:

> Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

**Zhang21 Method**: Polynomial functions of mean diameter (Dm) and water fraction (fx), with coefficients fitted from extensive T-matrix scattering simulations. At runtime, only polynomial evaluation is performed—**no table lookup**.

**TCWA2 Method**: Empirically-fitted polynomial formulations based on gamma distribution PSD parameters (λ, α). The polynomials were derived from T-matrix scattering calculations—**no coefficient files needed at runtime**.

The operator computes dual-polarization radar variables including:
- **Zhh**: Horizontal reflectivity (dBZ)
- **ZDR**: Differential reflectivity (dB)
- **KDP**: Specific differential phase (deg/km)
- **ρhv**: Copolar correlation coefficient

## Features

- **Multiple Polarimetric Operators:**
  - **Zhang21**: Polynomial functions fitted from T-matrix simulations (Zhang et al. 2021)
  - **TCWA2**: Empirically-fitted polynomial formulations based on gamma DSD (Tsai et al. 2025)
  
- **Microphysics Schemes:**
  - **Thompson**: Double-moment for rain (with Zhang21 operator)
  - **WSM6**: Single-moment (with Zhang21 operator)
  - **NSSL**: Double-moment with hail (with Zhang21 operator)
  - **TCWA2**: Taiwan CWA scheme with fitted polynomials (with TCWA2 operator)
  
- **Radar Frequencies:**
  - S-band (11 cm wavelength) and C-band (5.3 cm wavelength) support
  
- **Advanced Physics:**
  - Melting layer treatment (Zhang21: QC-based; TCWA2: melted fraction based)
  - Pure and melting hydrometeor categories (rain, snow, graupel, hail, ice)
  - Gamma PSD (Particle Size Distribution) parameterization (TCWA2)

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

**Using Zhang21 operator (polynomial-based from T-matrix):**
```fortran
use dualpol_op_mod, only: ppro_init_coefs, ppro_compute_point

! Initialize with Zhang21 operator (reads polynomial coefficient files)
call ppro_init_coefs(operator_name='Zhang21')

! Compute dual-pol variables using polynomial functions of Dm and fx
call ppro_compute_point(iband, 'THOMPSON', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, nr=nr)
```

**Using TCWA2 operator (fitted polynomials based on gamma PSD, no coefficient files needed):**
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

**TCWA2 operator** does not require coefficient files as it uses built-in fitted polynomial formulations.

**Note**: Coefficient files for Zhang21 must be obtained separately and are not included in this repository.

## Module Structure

### Main Interface (`dualpol_op_mod.f90`)

Unified interface providing:
- `ppro_init_coefs()`: Initialize library with operator selection
- `ppro_compute_point()`: Main computation (dispatches to selected operator)
- `ppro_set_operator()`: Switch between operators
- `ppro_finalize()`: Cleanup

### Zhang21 Operator (`zhang21/`)

**zhang21_core_mod.f90** - Polynomial-based operator:
- `zhang21_init_coefs()`: Load polynomial coefficients (S-band/C-band)
- `zhang21_compute_point()`: Compute using polynomial functions of Dm and fx
- `melting_scheme_zhang24()`: Melting layer treatment
- `dualpol_op_rain/icephase/total()`: Individual hydrometeor calculations

**zhang21_tlad_mod.f90** - Tangent-linear and adjoint for DA applications

### TCWA2 Operator (`tcwa2/`)

**tcwa2_core_mod.f90** - Empirically-fitted polynomial formulations based on gamma PSD:
- `dualpol_op_rain_tcwa2()`: Rain using fitted polynomial formulas
- `dualpol_op_ice_tcwa2()`: Ice using fitted polynomial formulas
- `dualpol_op_snow_tcwa2()`: Snow with melting treatment
- `dualpol_op_graup_tcwa2()`: Graupel with melting treatment
- `GAMLN()`: Gamma function computation

## Detailed Operator Usage

### Zhang21 Operator with Multiple Microphysics Schemes

**Method**: Polynomial functions of Dm (mean diameter) and fx (water fraction) fitted from T-matrix simulations

**Compatible Schemes**: WSM6, Thompson, NSSL

**WSM6 (Single-moment, 7 variables):**
```fortran
use dualpol_op_mod

call ppro_init_coefs(operator_name='Zhang21')
call ppro_compute_point(iband, 'WSM6', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv)
```

**Thompson (Double-moment for rain, 8 variables):**
```fortran
use dualpol_op_mod

call ppro_init_coefs(operator_name='Zhang21')
call ppro_compute_point(iband, 'THOMPSON', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, nr=nr)
```

**NSSL (Double-moment with hail, 14 variables):**
```fortran
use dualpol_op_mod

call ppro_init_coefs(operator_name='Zhang21')
call ppro_compute_point(iband, 'NSSL', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, &
                        qh=qh, nr=nr, ns=ns, ng=ng, nh=nh, vg=vg, vh=vh)
```

### TCWA2 Operator

**Method**: Empirically-fitted polynomial formulations based on gamma DSD

**Compatible Scheme**: TCWA2 only

**TCWA2 (15 variables with melting fractions):**
```fortran
use dualpol_op_mod

! Initialize TCWA2 operator (no coefficient files needed)
call ppro_init_coefs(operator_name='TCWA2')

! Compute dual-pol variables with TCWA2-specific inputs
! Note: Requires qi, ni (ice), qc (cloud water), smlf, gmlf (melting fractions)
call ppro_compute_point(iband, 'TCWA2', density_air, temp_air, &
                        qr, qs, qg, zh, zdr, kdp, phv, &
                        nr=nr, ns=ns, ng=ng, ni=ni, qi=qi, qc=qc, &
                        smlf=smlf, gmlf=gmlf)
```

**Key Differences**:
- **Zhang21**: 
  - Initialization: Read polynomial coefficients from files (one-time)
  - Runtime: Polynomial evaluation only—**no table lookup**
  - Melting: Temperature-based QC applied internally
  - Schemes: WSM6, Thompson, NSSL
  
- **TCWA2**: 
  - Initialization: No coefficient files needed
  - Runtime: Empirically-fitted polynomial evaluation (based on gamma PSD parameters)
  - Melting: Explicit melted fraction fields (smlf, gmlf) from model
  - Schemes: TCWA2 only

## Version History

- **v1.0.0** (2025): Initial modularized release
  - Extracted from JEDI-UFO as independent library
  - Multi-operator architecture (Zhang21 + TCWA2)
  - Modular directory structure (zhang21/ and tcwa2/ subdirectories)
  - Support for Thompson, WSM6, and NSSL microphysics
  - S-band and C-band radar support
  - Standalone examples and comprehensive documentation

## Authors

**Core Development:**
- **Rong Kong** (NCAR/MMM) - Library architecture, modularization, multi-operator framework, testing and maintenance

**Operator Implementation:**
- **Zhiquan (Jake) Liu** (NCAR/MMM) - Zhang21 operator implementation
- **Tzu-Chin Tsai** (CWA) - TCWA2 operator implementation

**Extensions and Integration:**
- **Hejun Xie** - JEDI-UFO integration
- **Tao Sun** - Hail categories and extended microphysics support

## License

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

## Operator Comparison: Zhang21 vs TCWA2

| Feature | Zhang21 | TCWA2 |
|---------|---------|-------|
| **Method** | Polynomial functions fitted from T-matrix simulations | Empirically-fitted polynomial formulations |
| **Basis** | Polynomials of Dm and fx (water fraction) | Fitted polynomials of gamma PSD parameters (λ, α) |
| **Coefficient files** | Required (fix/*.txt) | Not required |
| **Microphysics** | Thompson, WSM6, NSSL | TCWA2 (Taiwan CWA) |
| **Melting treatment** | Temperature-based QC | Melted fraction fields |
| **Required inputs** | qr, qs, qg [+ qh, nr, ns, ng, nh, vg, vh for NSSL] | qr, qs, qg, qc, qi, nr, ns, ng, ni, smlf, gmlf |
| **Ice crystal** | Not used | Explicit (qi, ni) |
| **Cloud water** | Not used | Explicit (qc) |
| **Validation** | Extensive | CWA operational |

**Key Algorithmic Differences**:
- **Zhang21**: Dm and fx → Polynomial evaluation → Zh, ZDR, KDP, ρhv
- **TCWA2**: Gamma PSD (λ, α) → Fitted polynomial evaluation → Zh, ZDR, KDP

**Recommendation**: 
- Use **Zhang21** for standard WRF microphysics schemes (Thompson, WSM6, NSSL)
- Use **TCWA2** for Taiwan CWA microphysics with explicit melting information

## References

1. **Zhang21 Operator**:  
   Zhang, G., J. Gao, and M. Du, 2021: Parameterized forward operators for simulation and assimilation of polarimetric radar data with numerical weather predictions. *Adv. Atmos. Sci.*, **38**(5), 737−754.

2. **Zhang21 Melting Scheme**:  
   Liu, et al., 2024: A New Melting Model and Its Implementation in Parameterized Forward Operators for Polarimetric Radar Data Simulation With Double Moment Microphysics Schemes. *JGR Atmospheres*, **129**.

3. **TCWA2 Operator**:  
   Tsai, T.-C., J.-P. Chen, Z. Liu, S.-Y. Jiang, R. Kong, Y.-J. Wu, J. Ban, L.-F. Hsiao, Y.-S. Tang, P.-L. Chang, and J.-S. Hong, 2025: Development of the TCWA2 Bulk Cloud Microphysics Scheme and Its Integration with a Dual-Polarization Radar Operator for Forecasting Applications. *J. Adv. Model. Earth Syst.*, submitted.

## Contact

For questions, issues, or contributions:
- **Technical Lead**: Rong Kong (rkong@ucar.edu)
- **Issues**: Please open an issue in the repository

## Acknowledgments

This work was developed at the National Center for Atmospheric Research (NCAR), Mesoscale and Microscale Meteorology Laboratory (MMM).

PPRO is designed as a general-purpose, standalone dual-polarization radar operator library. While it is compatible with data assimilation frameworks such as JEDI, it can be used independently in any numerical weather prediction system or forward operator application.



