# PPRO Examples

This directory contains standalone example programs demonstrating how to use the PPRO library independently.

## Available Examples

### `standalone_test.f90`

A comprehensive test program that demonstrates:
- How to initialize the PPRO library
- Computing dual-pol variables for different microphysics schemes (WSM6, Thompson, NSSL)
- Temperature sensitivity analysis (melting layer effects)
- Proper usage of the API for both single-moment and double-moment schemes

## Building Examples

Examples are built by default when you build the PPRO library:

```bash
cd ppro
mkdir build && cd build
cmake ..
make
```

To disable building examples:

```bash
cmake -DBUILD_EXAMPLES=OFF ..
```

## Running the Standalone Test

After building, the executable will be in `build/bin/`:

```bash
cd build/bin
./ppro_standalone_test
```

**Important**: Make sure coefficient files are available in one of these locations:
- `ppro/fix/` (relative to where you run the program)
- `./fix/` (current directory)
- Path specified by `PPRO_COEF_PATH` environment variable

Example:
```bash
export PPRO_COEF_PATH=/path/to/ppro/fix
./ppro_standalone_test
```

## Expected Output

The test program will output:
1. Initialization confirmation
2. Test results for WSM6 scheme
3. Test results for Thompson scheme
4. Test results for NSSL scheme
5. Temperature sensitivity table

Example output:
```
==============================================
PPRO Standalone Test Program
==============================================

Initializing PPRO coefficients...
Coefficients initialized successfully!

----------------------------------------------
Test Case 1: WSM6 (Single-moment)
----------------------------------------------
  Temperature:      280.000 K
  Density:            1.000 kg/m^3
  Qr (rain):          1.0000 g/kg
  Qs (snow):          0.5000 g/kg
  Qg (graupel):       0.3000 g/kg

  Zhh:              42.35 dBZ
  ZDR:               1.234 dB
  KDP:               0.5678 deg/km
  ρhv:               0.9850
...
```

## Using PPRO in Your Own Program

Here's a minimal example:

```fortran
program my_program
  use dualpol_op_mod, only: ppro_init_coefs, ppro_compute_point
  implicit none
  
  real(kind=8) :: density_air, temp_air, qr, qs, qg
  real(kind=8) :: zh, zdr, kdp, phv
  
  ! Initialize (call once at program start)
  call ppro_init_coefs()
  
  ! Set your atmospheric state
  density_air = 1.0_8    ! kg/m^3
  temp_air = 280.0_8     ! K
  qr = 1.0_8             ! g/kg
  qs = 0.5_8             ! g/kg
  qg = 0.3_8             ! g/kg
  
  ! Compute dual-pol variables
  ! iband: 1=S-band, 2=C-band
  call ppro_compute_point(1, 'WSM6', density_air, temp_air, &
                          qr, qs, qg, zh, zdr, kdp, phv)
  
  ! Convert to dBZ and dB for output
  zh = 10.0_8 * log10(zh + 1.0e-10_8)
  zdr = 10.0_8 * log10(zdr + 1.0e-10_8)
  
  print*, 'Zhh =', zh, ' dBZ'
  print*, 'ZDR =', zdr, ' dB'
  print*, 'KDP =', kdp, ' deg/km'
  print*, 'ρhv =', phv
  
end program my_program
```

## Troubleshooting

### "Cannot open coefficient file"

Ensure coefficient files are in the correct location. The library searches for files in this order:
1. Current directory (`./fix/`)
2. Relative path from current directory
3. Path from `PPRO_COEF_PATH` environment variable

### Linking errors

If you get linking errors when compiling your own program:
```bash
gfortran my_program.f90 -I/path/to/ppro/include -L/path/to/ppro/lib -lppro
```

Or use CMake:
```cmake
find_package(PPRO REQUIRED)
target_link_libraries(my_program PRIVATE ppro)
```

