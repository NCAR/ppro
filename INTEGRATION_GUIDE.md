# PPRO-LIB Integration Guide for JEDI-UFO

This guide describes how to integrate the modularized ppro-lib into JEDI-UFO.

## Overview of Modularization

The P-PRO radar operator has been refactored following the CRTM integration pattern:

**Before (Monolithic)**:
```
ufo/src/ufo/operators/radarreflectivity/ppro/
├── ObsPPRO.cc/.h                    (C++ interface)
├── ObsPPRO.interface.F90/.h         (Fortran-C++ bridge)
├── ufo_PPRO_mod.F90                 (High-level operator)
├── dualpol_op_mod.f90              (Core physics) ← MOVED
└── dualpol_op_tlad_mod.f90         (TL/AD)        ← MOVED
```

**After (Modular)**:
```
ppro-lib/                           (External library)
└── src/
    ├── dualpol_op_mod.f90         (Core physics)
    └── dualpol_op_tlad_mod.f90    (TL/AD)

ufo/src/ufo/operators/radarreflectivity/ppro/
├── ObsPPRO.cc/.h                   (C++ interface)
├── ObsPPRO.interface.F90/.h        (Fortran-C++ bridge)
└── ufo_PPRO_mod.F90                (High-level operator - links to ppro-lib)
```

## Integration Steps

### Step 1: Build ppro-lib

```bash
cd ppro-lib
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make
make install
```

### Step 2: Modify UFO CMakeLists.txt

Add the following to `ufo/CMakeLists.txt` or `ufo/src/ufo/CMakeLists.txt`:

#### Option A: Using find_package (Production)

```cmake
# Find ppro library
find_package(PPRO REQUIRED)

# Link UFO with ppro
target_link_libraries(ufo PUBLIC ppro)
```

Then configure UFO with:
```bash
cmake -DPPRO_ROOT=/path/to/ppro/install ..
```

#### Option B: Using add_subdirectory (Development)

Place ppro-lib in your workspace and add to UFO's CMakeLists.txt:

```cmake
# Add ppro as a subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/../ppro-lib ppro-lib)

# Link UFO with ppro
target_link_libraries(ufo PUBLIC ppro)
```

### Step 3: Update Module Paths

The Fortran modules from ppro-lib will be automatically available through the CMake target interface.

If you need to explicitly set paths, add:
```cmake
target_include_directories(ufo PUBLIC 
    ${PPRO_ROOT}/include/ppro
)
```

### Step 4: Verify the Integration

Build UFO:
```bash
cd ufo/build
cmake ..
make
```

Check that:
1. The ppro library is found
2. Fortran modules `dualpol_op_mod` and `dualpol_op_tlad_mod` are accessible
3. UFO compiles without errors

## Code Changes Required in UFO

### In ufo_PPRO_mod.F90

The module should already have:
```fortran
use dualpol_op_mod
```

This will now link to the external ppro-lib instead of a local module.

**No code changes are needed** in `ufo_PPRO_mod.F90` if the module names remain the same.

### In ufo_PPRO_tlad_mod.F90

Similarly, the TL/AD module uses:
```fortran
use dualpol_op_tlad_mod
```

This will automatically link to ppro-lib.

## Directory Structure After Integration

```
workspace/
├── ppro-lib/                      # External library
│   ├── CMakeLists.txt
│   ├── src/
│   │   ├── dualpol_op_mod.f90
│   │   └── dualpol_op_tlad_mod.f90
│   └── fix/                      # Coefficient files
│
└── ufo/
    ├── CMakeLists.txt            # Modified to link ppro
    └── src/ufo/operators/radarreflectivity/ppro/
        ├── CMakeLists.txt        # Modified (removed core modules)
        ├── ObsPPRO.cc/.h
        ├── ObsPPRO.interface.F90
        └── ufo_PPRO_mod.F90      # Uses ppro-lib modules
```

## Coefficient Files

The ppro-lib requires coefficient files for S-band and C-band radars.

### Option 1: Absolute Path

Set the coefficient file path at runtime or in the YAML configuration:
```yaml
obs operator:
  name: PPRO
  coefficient_path: /path/to/ppro-lib/fix
```

### Option 2: Environment Variable

```bash
export PPRO_COEF_PATH=/path/to/ppro-lib/fix
```

### Option 3: Relative Path

Place coefficient files in the run directory or set paths in `dualpol_op_mod.f90`:
```fortran
call read_coefs_rain('/path/to/sband_rain_coefs.txt', sband_rain_coefs)
```

## Testing

After integration, test with:

1. **Build test**: Ensure UFO compiles with ppro-lib
2. **Unit test**: Run existing P-PRO operator tests
3. **Integration test**: Run JEDI applications using P-PRO

Example test command:
```bash
ctest -R ppro
```

## Troubleshooting

### Module Not Found Error

```
Error: Can't open module file 'dualpol_op_mod.mod'
```

**Solution**: Ensure ppro-lib is properly installed and the module path is included:
```cmake
target_include_directories(ufo PUBLIC ${PPRO_ROOT}/include/ppro)
```

### Linking Error

```
undefined reference to `dualpol_op_mod_mp_init_coefs_'
```

**Solution**: Ensure ppro library is linked:
```cmake
target_link_libraries(ufo PUBLIC ppro)
```

### Coefficient Files Not Found

```
Error opening file: sband_rain_coefs.txt
```

**Solution**: 
1. Check that coefficient files exist in `ppro-lib/fix/`
2. Verify the file paths in `dualpol_op_mod.f90`
3. Run from the correct directory or use absolute paths

## Comparison with CRTM Integration

The ppro-lib integration follows the same pattern as CRTM in UFO:

| Aspect | CRTM | PPRO |
|--------|------|------|
| External library | `crtm-lib` | `ppro-lib` |
| UFO interface | `ufo_radiancecrtm_mod.F90` | `ufo_PPRO_mod.F90` |
| Core modules | CRTM_Module | dualpol_op_mod |
| Integration | `find_package(CRTM)` | `find_package(PPRO)` |
| Coefficient files | `crtm-lib/fix/` | `ppro-lib/fix/` |

## Benefits of Modularization

1. **Separation of Concerns**: Physics code separate from UFO interface
2. **Reusability**: ppro-lib can be used in other projects
3. **Maintainability**: Easier to update core physics independently
4. **Testing**: Core algorithms can be tested independently
5. **Version Control**: ppro-lib can have its own version/release cycle

## Future Enhancements

1. Add Python bindings for standalone testing
2. Create regression test suite for ppro-lib
3. Add more microphysics schemes (e.g., Morrison, P3)
4. Optimize performance with vectorization/GPU support
5. Add automated coefficient file generation tools

## Contact

For questions about integration, contact:
rkong@ucar.edu

