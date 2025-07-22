# 2S-ESS Project Documentation Index

## Overview
**2S-ESS** (Two-Stream + Exact Single Scattering) is a scientific computational library for atmospheric radiative transfer modeling, combining two major components for high-accuracy atmospheric radiation calculations.

### Quick Navigation
- [📖 Core Documentation](#core-documentation)
- [🏗️ Architecture](#architecture)
- [🧰 API Reference](#api-reference)
- [🚀 Getting Started](#getting-started)
- [🔧 Build System](#build-system)
- [🧪 Testing](#testing)
- [📚 Examples](#examples)
- [👥 Contributing](#contributing)

---

## Core Documentation

### Primary Documents
| Document | Purpose | Path |
|----------|---------|------|
| **README.md** | Project overview and quick start | `/README.md` |
| **CLAUDE.md** | Codebase documentation for AI tools | `/CLAUDE.md` |
| **Technical Documentation** | Detailed 2STREAM API reference | `/docs/2STREAM_v2.4_Technical_Documentation.md` |
| **LICENSE** | GNU GPL v3.0 license | `/LICENSE` |

### Scientific Papers
| Paper | Topic | Path |
|-------|-------|------|
| **JQSRT Paper** | Two-stream discrete ordinate method | `/docs/Spurr_Natraj_JQSRT_2stream_2011.pdf` |
| **First Order Guide** | ESS implementation guide | `/docs/first_order_1p3_userguide_v5_09apr2013.pdf` |

---

## Architecture

### Library Structure
```
2S-ESS/
├── 📁 sourcecode/           # Core implementation
│   ├── 2stream_2p4_feb15/   # Original 2STREAM v2.4
│   ├── 2stream_2p4_optimized/ # Optimized 2STREAM v2.4
│   │   ├── solar/          # Solar-specific modules
│   │   │   ├── lat/        # Latitude-dependent
│   │   │   └── obs/        # Observation geometry
│   │   └── thermal/        # Thermal emission
│   ├── FO_1p4_sourcecode/  # Original ESS v1.4
│   └── FO_1p4_optimized/   # Optimized ESS v1.4
├── 📁 TEST/                 # Test programs
├── 📁 docs/                 # Documentation
└── 📁 build/                # Build artifacts (generated)
```

### Core Components

#### 1. Two-Stream Radiative Transfer (2STREAM v2.4)
**Purpose**: Discrete ordinate method for atmospheric radiative transfer
**Authors**: Robert J. D. Spurr (RT Solutions, Inc.), Vijay Natraj (JPL/Caltech)

##### Key Features:
- Two-stream discrete ordinate radiative transfer model
- Solar and thermal radiation support
- BRDF (Bidirectional Reflectance Distribution Function) modeling
- Linearized Jacobian calculations for sensitivity analysis
- Enhanced spherical correction for improved accuracy
- Optimized pentadiagonal matrix solver

##### Module Classification:
| Category | Modules | Purpose |
|----------|---------|---------|
| **Master Interfaces** | `*master*.f90` | Main entry points |
| **Computational Kernels** | `*solutions*.f90`, `*intensity*.f90` | Core calculations |
| **Boundary Value** | `*bvproblem*.f90` | Matrix solvers |
| **Physical Processes** | `*brdf*.f90`, `*thermal*.f90` | Surface & thermal |
| **Linearization** | `*lp_*.f90`, `*lc_*.f90` | Jacobian calculations |
| **Support** | `*inputs*.f90`, `*miscsetups*.f90` | Utilities |

#### 2. First Order Exact Single Scattering (ESS v1.4)
**Purpose**: First-order single scattering calculations with high accuracy

##### Key Features:
- Exact single scattering calculations
- Solar and thermal source support
- Observational geometry configurations
- Spherical and pseudo-spherical corrections
- Optimized computational algorithms

---

## API Reference

### Primary Interfaces

#### TWOSTREAM_MASTER
**Location**: `sourcecode/2stream_2p4_feb15/2stream_master.f90`
**Purpose**: Main interface for two-stream calculations
**Signature**: 
```fortran
SUBROUTINE TWOSTREAM_MASTER(
    ! Dimensions, control flags, geometry, optical properties,
    ! surface properties, outputs, error handling
)
```

#### Linearization Interfaces
| Interface | Purpose | Location |
|-----------|---------|----------|
| **TWOSTREAM_LPS_MASTER** | Profile Jacobians (∂I/∂x per layer) | `*lps_master*.f90` |
| **TWOSTREAM_LCS_MASTER** | Column Jacobians (∂I/∂X integrated) | `*lcs_master*.f90` |

#### First Order Interfaces
| Interface | Purpose | Location |
|-----------|---------|----------|
| **FO_ScalarSS_RTCalcs** | Single scattering radiances | `FO_*RTCalcs*.f90` |
| **FO_Thermal_RTCalcs** | Thermal single scattering | `FO_Thermal*.f90` |
| **FO_Geometry_Master** | Geometry preprocessing | `FO_*geometry*.f90` |

### Key Parameters

#### Dimensional Parameters
- **MAXLAYERS**: Maximum atmospheric layers (114)
- **MAXTOTAL**: Total BVP dimension (2 × MAXLAYERS)
- **MAX_GEOMETRIES**: Maximum observation geometries

#### Control Flags
- **DO_UPWELLING/DO_DNWELLING**: Calculation direction
- **DO_SOLAR_SOURCES**: Include solar radiation
- **DO_THERMAL_EMISSION**: Include thermal emission
- **DO_BRDF_SURFACE**: Use BRDF surface model
- **DO_USER_OBSGEOMS**: Observational geometry mode

#### Optical Properties
- **DELTAU_INPUT**: Layer optical depths
- **OMEGA_INPUT**: Single scattering albedos
- **ASYMM_INPUT**: Asymmetry parameters

---

## Getting Started

### System Requirements
- **Fortran Compiler**: Intel Fortran (ifort) or GNU Fortran (gfortran)
- **Build System**: CMake 3.20+ or Make
- **Dependencies**: numdiff (for testing)
- **Standards**: Fortran 90/95 with legacy extensions

### Quick Build
```bash
# Using CMake (recommended)
cmake -B build
cd build
cmake --build .
ctest

# Using traditional makefiles
cd TEST
make -f makefile        # Regular version
make -f makefile_opt    # Optimized version
```

### Compiler Configuration
```bash
# Intel Fortran
ifort -fpp -D IFORT -O3 -heap-arrays

# GNU Fortran
gfortran -std=legacy -mcmodel=large -O3 -ffree-form -ffree-line-length-none
```

---

## Build System

### CMake Configuration
**File**: `CMakeLists.txt`
**Version**: CMake 3.20+

#### Build Targets
| Target | Type | Source | Description |
|--------|------|--------|-------------|
| **TWOSTR** | Static Library | `2stream_2p4_feb15/` | Original 2STREAM |
| **TWOSTR_OPT** | Static Library | `2stream_2p4_optimized/` | Optimized 2STREAM |
| **ESS** | Static Library | `FO_1p4_sourcecode/` | Original ESS |
| **ESS_OPT** | Static Library | `FO_1p4_optimized/` | Optimized ESS |
| **regtest** | Executable | `TEST/test_regular.f90` | Regular version test |
| **opttest** | Executable | `TEST/test_opt.f90` | Optimized version test |

#### Compiler Flags
- **Intel**: `-fpp -D IFORT -O3 -heap-arrays`
- **GNU**: `-std=legacy -mcmodel=large -O3 -ffree-form`
- **Debug**: Adds bounds checking and runtime verification

### Traditional Makefiles
**Location**: `TEST/makefile*`
**Usage**: Simplified build for testing

---

## Testing

### Test Programs
| Executable | Library | Purpose |
|------------|---------|---------|
| **regtest** | TWOSTR + ESS | Tests original implementation |
| **opttest** | TWOSTR_OPT + ESS_OPT | Tests optimized implementation |

### Test Configuration
**Atmospheric Setup**:
- **Layers**: 114 (1 km each, TOA at 114 km)
- **Optical Properties**: τ=0.01/layer, ω=0.5, g=0.9
- **Surface**: Lambertian albedo = 0.3
- **Geometry**: SZA=30°, VZA=0° (nadir), AZM=0° (principal plane)
- **Features**: Enhanced spherical correction enabled

### Validation
Tests compare output against expected results in `TEST/expected_test_results.dat` using numerical difference checking with tolerance `1e-14`.

### Test Commands
```bash
# Run all tests
ctest

# Run specific tests
ctest -L reg      # Regular tests only
ctest -L opt      # Optimized tests only
ctest --verbose   # Detailed output
```

---

## Examples

### Basic Usage Pattern
```fortran
program example
    ! 1. Setup dimensions and control flags
    ! 2. Configure geometry and optical properties
    ! 3. Call TWOSTREAM_MASTER
    ! 4. Check error status
    ! 5. Process results
end program
```

### Key Usage Examples
- **Basic Solar Calculation**: Standard atmospheric radiative transfer
- **Profile Jacobian**: Sensitivity analysis for retrievals
- **BRDF Surface**: Complex surface reflection modeling
- **Thermal Emission**: Infrared radiation calculations

---

## File Reference

### Source Code Organization

#### Core Fortran Modules (114 files total)
```
Regular Implementation (2stream_2p4_feb15/):
├── 2stream_master.f90              # Main interface
├── 2stream_solutions.f90           # Homogeneous solutions
├── 2stream_bvproblem.f90           # Boundary value problem
├── 2stream_intensity.f90           # Radiance calculations
├── 2stream_inputs.f90              # Input validation
├── 2stream_miscsetups.f90          # Preprocessing
└── [Additional modules...]

Optimized Implementation (2stream_2p4_optimized/):
├── solar/obs/                      # Observation geometry
├── solar/lat/                      # Latitude configurations
├── thermal/                        # Thermal emission
└── [Optimized modules...]

First Order (FO_1p4_*/):
├── FO_ScalarSS_RTCalcs_*.f90      # Single scattering
├── FO_geometry_*.f90              # Geometry handling
└── [FO modules...]
```

### Performance Comparison
| Version | Performance | Memory | Accuracy |
|---------|-------------|--------|----------|
| **Regular** | Baseline | Standard | Reference |
| **Optimized** | 3-5x faster | Reduced | Identical |

---

## Contributing

### Development Workflow
1. **Environment Setup**: Configure Fortran compiler and CMake
2. **Build Testing**: Verify both regular and optimized versions
3. **Code Standards**: Follow existing Fortran 90/95 conventions
4. **Testing**: Run full test suite with numerical validation
5. **Documentation**: Update relevant documentation

### Contact Information
- **Primary Contact**: Vijay Natraj (vijay.natraj@jpl.nasa.gov)
- **Original Author**: Robert J. D. Spurr (RT Solutions, Inc.)
- **Institution**: Jet Propulsion Laboratory, California Institute of Technology

### License
**GNU General Public License v3.0** - See `LICENSE` file for full terms.

---

## Version Information
- **2S-ESS Version**: 2.4.2
- **2STREAM Version**: 2.4 (February 2015)
- **ESS Version**: 1.4
- **Documentation Generated**: 2025-01-22

*This index provides comprehensive navigation for the 2S-ESS atmospheric radiative transfer library. For detailed technical information, refer to the linked documentation and source code.*