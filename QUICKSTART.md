# PPRO-LIB Quick Start Guide

## 5-Minute Setup

### Prerequisites
- CMake 3.20+
- Fortran compiler (gfortran, ifort, or ifx)

### Build ppro

```bash
cd ppro
mkdir build && cd build
cmake ..
make
```

This creates `libppro.a` in `ppro/build/lib/`

### Use with UFO

#### Option 1: Development Mode (Recommended)

Place ppro next to UFO:
```
workspace/
├── ppro/
└── ufo/
```

In `ufo/CMakeLists.txt`, add:
```cmake
add_subdirectory(${CMAKE_SOURCE_DIR}/../ppro ppro)
target_link_libraries(ufo PUBLIC ppro)
```

Build UFO:
```bash
cd ufo/build
cmake ..
make
```

#### Option 2: Installed Library

Install ppro:
```bash
cd ppro/build
make install  # Installs to CMAKE_INSTALL_PREFIX
```

Configure UFO:
```bash
cd ufo/build
cmake -DPPRO_ROOT=/path/to/ppro/install ..
make
```

## Testing

Run UFO tests:
```bash
cd ufo/build
ctest -R ppro -V
```

## What's Next?

1. **Add coefficient files**: Copy `*_coefs.txt` to `ppro/fix/`
2. **Read the full guide**: See `INTEGRATION_GUIDE.md` for details
3. **Check documentation**: See `README.md` for API reference

## Troubleshooting

**Problem**: Module not found error
```
Error: Can't open module file 'dualpol_op_mod.mod'
```
**Solution**: Ensure ppro is linked properly in UFO's CMakeLists.txt

**Problem**: Library not found
```
undefined reference to `__dualpol_op_mod_MOD_init_coefs`
```
**Solution**: Add `target_link_libraries(ufo PUBLIC ppro)` to UFO's CMakeLists.txt

## Support

- **Technical Lead**: Rong Kong (rkong@ucar.edu)
- **Issues**: Open a GitHub issue
- **Documentation**: See `README.md` and `INTEGRATION_GUIDE.md`


