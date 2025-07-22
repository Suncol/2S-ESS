# 2S-ESS Build & Development Guide

## Quick Navigation
- [🚀 Quick Start](#quick-start)
- [🏗️ Build Systems](#build-systems)
- [🧰 Development Setup](#development-setup)
- [🧪 Testing](#testing)
- [⚙️ Configuration](#configuration)
- [🔧 Troubleshooting](#troubleshooting)

---

## Quick Start

### Prerequisites
```bash
# Required
- Fortran compiler (gfortran or ifort)
- CMake 3.20+ (recommended) or Make
- numdiff (for test validation)

# System packages (Ubuntu/Debian)
sudo apt-get install gfortran cmake numdiff

# System packages (RHEL/CentOS)
sudo yum install gcc-gfortran cmake numdiff
```

### 5-Minute Build
```bash
# Clone and build
git clone <repository-url> 2S-ESS
cd 2S-ESS

# Option 1: CMake (recommended)
cmake -B build
cd build
cmake --build .
ctest

# Option 2: Traditional makefiles
cd TEST
make -f makefile
./regtest
```

---

## Build Systems

### CMake Build (Recommended)

#### Basic Build
```bash
# Configure
cmake -B build

# Build
cmake --build build

# Install
cmake --install build
```

#### Advanced Configuration
```bash
# Debug build
cmake -B build -DCMAKE_BUILD_TYPE=Debug

# Release build with specific compiler
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=ifort

# Shared libraries
cmake -B build -DENABLE_SHARED=ON

# Custom installation prefix
cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local
```

#### Build Targets
```bash
# Build specific targets
cmake --build build --target TWOSTR          # Regular 2STREAM library
cmake --build build --target TWOSTR_OPT      # Optimized 2STREAM library
cmake --build build --target ESS             # Regular ESS library
cmake --build build --target ESS_OPT         # Optimized ESS library
cmake --build build --target regtest         # Regular test executable
cmake --build build --target opttest         # Optimized test executable
```

### Traditional Makefiles

#### Location and Usage
```bash
cd TEST/

# Build regular version
make -f makefile

# Build optimized version
make -f makefile_opt

# Clean build artifacts
make -f makefile clean
make -f makefile_opt clean

# Run tests
./regtest                    # Regular test
./opttest                    # Optimized test
```

#### Makefile Structure
```bash
# Main targets in TEST/makefile
regtest: test_regular.f90 + TWOSTR + ESS libraries
opttest: test_opt.f90 + TWOSTR_OPT + ESS_OPT libraries
```

---

## Development Setup

### Compiler Configuration

#### GNU Fortran (gfortran)
```bash
# Standard flags
gfortran -std=legacy -mcmodel=large -O3 -ffree-form -ffree-line-length-none

# Debug flags
gfortran -std=legacy -mcmodel=large -g -O0 -fbounds-check -frange-check \
         -ffpe-trap=invalid,zero,overflow -Wall -fcheck=all -fbacktrace

# OpenMP support
gfortran -fopenmp [other flags]
```

#### Intel Fortran (ifort)
```bash
# Standard flags
ifort -fpp -D IFORT -O3 -heap-arrays

# Debug flags  
ifort -fpp -D IFORT -g -warn all -check all -zero

# OpenMP support
ifort -openmp [other flags]
```

### Environment Setup
```bash
# Set compiler
export FC=gfortran          # or ifort
export CMAKE_Fortran_COMPILER=gfortran

# Large model support
export FCFLAGS="-mcmodel=large"

# CMake configuration
export CMAKE_BUILD_TYPE=Release    # or Debug
```

### IDE Integration

#### VS Code Setup
```json
// .vscode/settings.json
{
    "cmake.configureSettings": {
        "CMAKE_Fortran_COMPILER": "gfortran"
    },
    "cmake.buildDirectory": "${workspaceFolder}/build"
}
```

#### Modern Fortran Extension
- Syntax highlighting
- Code formatting
- Symbol navigation
- Error detection

---

## Testing

### Test Structure
```
TEST/
├── makefile                      # Traditional build (regular)
├── makefile_opt                  # Traditional build (optimized) 
├── test_regular.f90              # Regular version test program
├── test_opt.f90                  # Optimized version test program
└── expected_test_results.dat     # Reference results
```

### CMake Testing
```bash
# Run all tests
ctest

# Run specific test types
ctest -L reg                      # Regular tests only
ctest -L opt                      # Optimized tests only
ctest -L run                      # Runtime tests only
ctest -L test                     # Output validation tests

# Verbose output
ctest --verbose

# Run specific test
ctest -R "regular run"
ctest -R "optimized output"
```

### Test Configuration
**Atmospheric Parameters**:
- **Layers**: 114 (identical optical properties)
- **Optical Depth**: 1.14 total (0.01 per layer)  
- **Single Scattering Albedo**: 0.5
- **Asymmetry Parameter**: 0.9 (Henyey-Greenstein)
- **Surface Albedo**: 0.3
- **Geometry**: SZA=30°, VZA=0°, AZM=0°
- **Layer Height**: 1 km (TOA at 114 km)
- **Spherical Correction**: Enhanced mode enabled

### Manual Testing
```bash
# Build and run tests manually
cd TEST
make -f makefile
./regtest > test_output.dat

# Compare with reference
numdiff -S -r 1e-14 expected_test_results.dat test_output.dat
```

### Custom Test Cases
```fortran
! Create custom test programs using the template:
program my_test
    use twostream_master_m
    ! [Setup parameters as needed]
    ! [Call TWOSTREAM_MASTER]
    ! [Output results]
end program
```

---

## Configuration

### CMake Options
```cmake
# Available options in CMakeLists.txt
option(ENABLE_SHARED "Enable shared libraries" OFF)

# Configuration variables
CMAKE_BUILD_TYPE          # Debug, Release, RelWithDebInfo
CMAKE_Fortran_COMPILER    # Compiler selection
CMAKE_INSTALL_PREFIX      # Installation directory
```

### Compiler-Specific Settings

#### Performance Optimization
```bash
# GNU Fortran - Maximum performance
FCFLAGS="-O3 -march=native -funroll-loops -ffast-math"

# Intel Fortran - Maximum performance  
FCFLAGS="-O3 -xHost -ipo -fast"

# Memory optimization
FCFLAGS="-mcmodel=large"          # GNU
FCFLAGS="-heap-arrays"            # Intel
```

#### Debug Configuration
```bash
# GNU Fortran - Full debugging
FCFLAGS="-g -O0 -fbounds-check -frange-check -ffpe-trap=all -Wall -Wextra"

# Intel Fortran - Full debugging
FCFLAGS="-g -O0 -check all -warn all -traceback -fpe0"
```

### Directory Structure Configuration
```bash
# CMake output directories
CMAKE_LIBRARY_OUTPUT_DIRECTORY = ${CMAKE_BINARY_DIR}/lib
CMAKE_ARCHIVE_OUTPUT_DIRECTORY = ${CMAKE_BINARY_DIR}/lib  
CMAKE_Fortran_MODULE_DIRECTORY = ${CMAKE_BINARY_DIR}/mod_files

# Installation structure
bin/           # Executables (regtest, opttest)
lib/           # Static libraries (TWOSTR, TWOSTR_OPT, ESS, ESS_OPT)
include/       # Fortran module files (.mod)
```

---

## Troubleshooting

### Common Build Issues

#### Compiler Not Found
```bash
# Error: No Fortran compiler found
# Solution: Install Fortran compiler and set FC environment variable
sudo apt-get install gfortran     # Ubuntu/Debian
export FC=gfortran
```

#### CMake Version Too Old
```bash
# Error: CMake 3.20 or higher required
# Solution: Update CMake
pip install cmake --upgrade
# or download from https://cmake.org/download/
```

#### Module File Issues
```bash
# Error: Module files not found
# Solution: Check module directory and compiler compatibility
export CMAKE_Fortran_MODULE_DIRECTORY=${PWD}/build/mod_files
```

### Runtime Issues

#### Segmentation Fault
```bash
# Possible causes and solutions:
1. Stack overflow: Use -heap-arrays (Intel) or increase stack size
   ulimit -s unlimited

2. Array bounds: Enable bounds checking in debug builds
   FCFLAGS="-fbounds-check"

3. Uninitialized variables: Use debug flags
   FCFLAGS="-finit-real=nan"
```

#### Numerical Issues
```bash
# Floating point exceptions
FCFLAGS="-ffpe-trap=invalid,zero,overflow"    # GNU
FCFLAGS="-fpe0"                               # Intel

# Numerical precision
# Use consistent real precision throughout (double precision recommended)
```

#### Memory Issues
```bash
# Large dataset handling
FCFLAGS="-mcmodel=large"                      # GNU
FCFLAGS="-heap-arrays -mcmodel=large"         # Intel

# OpenMP memory
export OMP_STACKSIZE=64M
```

### Test Failures

#### Output Mismatch
```bash
# Check numerical differences
numdiff -S -r 1e-14 expected_test_results.dat actual_output.dat

# Common causes:
1. Compiler differences (acceptable within tolerance)
2. Optimization level differences  
3. OpenMP race conditions
4. Platform-specific floating point behavior
```

#### Performance Issues
```bash
# Profile execution time
time ./regtest
time ./opttest

# Memory profiling
valgrind --tool=massif ./regtest
```

### Getting Help

#### Debugging Steps
1. **Clean build**: `rm -rf build/` and rebuild
2. **Debug build**: Use `CMAKE_BUILD_TYPE=Debug`
3. **Compiler warnings**: Enable all warnings
4. **Runtime checks**: Use debug compiler flags
5. **Version check**: Verify compiler and CMake versions

#### Log Collection
```bash
# CMake configuration log
cmake -B build --debug-output > cmake_debug.log 2>&1

# Build log with verbose output
cmake --build build --verbose > build.log 2>&1

# Test output with details
ctest --verbose > test.log 2>&1
```

#### Contact Information
- **Primary Contact**: vijay.natraj@jpl.nasa.gov
- **Issues**: Check compiler versions and provide full error logs
- **Performance**: Include timing comparisons and system specifications

---

## Advanced Topics

### Cross-Platform Development
```bash
# Windows (MinGW)
cmake -G "MinGW Makefiles" -B build

# macOS
cmake -B build
# May need: brew install gcc cmake

# HPC clusters
module load intel/19.0 cmake/3.20
export FC=ifort
```

### Integration with Other Projects
```fortran
! Using 2S-ESS as a library
use twostream_master_m
use fo_scalars_rtcalcs_m

! Link against libraries
# -lTWOSTR -lESS (regular)
# -lTWOSTR_OPT -lESS_OPT (optimized)
```

*For more detailed technical documentation, see `/docs/2STREAM_v2.4_Technical_Documentation.md`*