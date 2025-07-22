# 2STREAM v2.4 Technical Documentation

## Table of Contents
1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Core Interfaces](#core-interfaces)
4. [Module Reference](#module-reference)
5. [Calling Patterns](#calling-patterns)
6. [Mathematical Framework](#mathematical-framework)
7. [Implementation Guide](#implementation-guide)
8. [Performance Optimization](#performance-optimization)
9. [Error Handling](#error-handling)
10. [Examples](#examples)

---

## Overview

### Library Purpose
The 2STREAM v2.4 library implements a two-stream discrete ordinate radiative transfer model for atmospheric applications. It provides accurate and efficient computation of solar and thermal radiation through multi-layered atmospheric systems with support for:

- **Solar and thermal radiation sources**
- **Bidirectional Reflectance Distribution Functions (BRDF)**
- **Surface leaving radiance**
- **Linearized Jacobian calculations**
- **Spherical and plane-parallel geometries**

### Key Features
- **High Performance**: Optimized pentadiagonal matrix solver
- **Numerical Stability**: Taylor series expansions for optically thin layers
- **Linearization Support**: Profile and column Jacobian matrices
- **Flexible Geometry**: Observational geometry and lattice configurations
- **Multiple Sources**: Solar, thermal emission, and surface leaving terms

### Scientific Applications
- Satellite radiance simulation and validation
- Atmospheric parameter retrieval algorithms
- Climate model radiation schemes
- Ground-based remote sensing analysis

---

## Architecture

### Design Philosophy
The library follows a hierarchical modular design with clear separation of concerns:

```
Application Layer
    ↓
Master Interface Layer (TWOSTREAM_MASTER)
    ↓
Computational Kernel Layer (Solutions, BVP, Intensity)
    ↓
Physical Process Layer (BRDF, Thermal, Surface)
    ↓
Utility Layer (Inputs, Miscsetups, Writemodules)
```

### Module Classification

#### Core Computational Modules
| Module | Primary Function | Key Subroutines |
|--------|------------------|-----------------|
| `2stream_master.f90` | Main interface | `TWOSTREAM_MASTER`, `TWOSTREAM_FOURIER_MASTER` |
| `2stream_solutions.f90` | Homogeneous solutions | `TWOSTREAM_HOM_SOLUTION`, `TWOSTREAM_GBEAM_SOLUTION` |
| `2stream_bvproblem.f90` | Boundary value problem | `TWOSTREAM_BVP_SOLUTION_PENTADIAG` |
| `2stream_intensity.f90` | Radiance computation | `TWOSTREAM_UPUSER_INTENSITY`, `TWOSTREAM_DNUSER_INTENSITY` |

#### Linearization Modules
| Module | Purpose | Jacobian Type |
|--------|---------|---------------|
| `2stream_lps_master.f90` | Profile linearization | ∂I/∂(τ,ω,g) per layer |
| `2stream_lcs_master.f90` | Column linearization | ∂I/∂(total column) |
| `2stream_lp_jacobians.f90` | Profile Jacobian computation | Layer-by-layer sensitivities |
| `2stream_lc_jacobians.f90` | Column Jacobian computation | Integrated column sensitivities |

#### Physical Process Modules
| Module | Physical Process | Implementation |
|--------|------------------|----------------|
| `2stream_brdf_kernels.f90` | Surface reflection | Lambertian, Rahman, Cox-Munk, Hapke |
| `2stream_thermalsup.f90` | Thermal emission | Planck function, thermal sources |
| `2stream_sleave_*.f90` | Surface leaving | Ocean color, fluorescence |

#### Support Modules
| Module | Function | Key Features |
|--------|----------|--------------|
| `2stream_inputs.f90` | Input validation | Dimension checking, parameter limits |
| `2stream_miscsetups.f90` | Preprocessing | Geometry, optical properties, transmittances |
| `2stream_writemodules.f90` | Debug output | Formatted parameter dumps |

---

## Core Interfaces

### Primary Interface: TWOSTREAM_MASTER

#### Function Signature
```fortran
SUBROUTINE TWOSTREAM_MASTER ( &
    ! === Dimensions ===
    MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES, &
    MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS, &
    
    ! === Control Flags ===
    DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT, &
    DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, &
    DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, &
    DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS, &
    DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE, &
    
    ! === Numerical Control ===
    BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF, &
    
    ! === Geometry ===
    NLAYERS, NTOTAL, STREAM_VALUE, &
    N_USER_OBSGEOMS, USER_OBSGEOMS, &
    N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS, &
    FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID, &
    
    ! === Optical Properties ===
    DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING, &
    THERMAL_BB_INPUT, &
    
    ! === Surface Properties ===
    LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, &
    EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0, &
    
    ! === Outputs ===
    INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA, &
    RADLEVEL_UP, RADLEVEL_DN, N_GEOMETRIES, &
    
    ! === Error Handling ===
    STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS, &
    STATUS_EXECUTION, E_MESSAGE, E_TRACE_1, E_TRACE_2, &
    
    ! === Performance ===
    geom_timer &
)
```

#### Key Parameters

##### Dimensional Parameters
- **MAXLAYERS**: Maximum atmospheric layers (typically 114)
- **MAXTOTAL**: Total BVP dimension = 2 × MAXLAYERS
- **MAX_GEOMETRIES**: Maximum geometry configurations
- **MAX_USER_OBSGEOMS**: Maximum observational geometries

##### Control Flags
- **DO_UPWELLING/DO_DNWELLING**: Direction control
- **DO_PLANE_PARALLEL**: Geometry mode (plane-parallel vs. spherical)
- **DO_2S_LEVELOUT**: Output radiances at all atmospheric levels
- **DO_SOLAR_SOURCES**: Include solar radiation
- **DO_THERMAL_EMISSION**: Include thermal emission
- **DO_BRDF_SURFACE**: Use BRDF surface model
- **DO_USER_OBSGEOMS**: Observational geometry mode

##### Optical Properties
- **DELTAU_INPUT(MAXLAYERS)**: Layer optical depths
- **OMEGA_INPUT(MAXLAYERS)**: Single scattering albedos
- **ASYMM_INPUT(MAXLAYERS)**: Asymmetry parameters (Henyey-Greenstein)
- **D2S_SCALING(MAXLAYERS)**: Delta-2Stream scaling factors

##### Surface Properties
- **LAMBERTIAN_ALBEDO**: Lambertian surface albedo
- **BRDF_F_0(0:1,MAXBEAMS)**: BRDF Fourier components (solar incidence)
- **BRDF_F(0:1)**: BRDF Fourier components (quadrature streams)
- **EMISSIVITY**: Surface emissivity for thermal radiation

##### Outputs
- **INTENSITY_TOA(MAX_GEOMETRIES)**: Top-of-atmosphere radiances
- **INTENSITY_BOA(MAX_GEOMETRIES)**: Bottom-of-atmosphere radiances
- **FLUXES_TOA(MAXBEAMS,2)**: TOA fluxes (upward/downward)
- **FLUXES_BOA(MAXBEAMS,2)**: BOA fluxes (upward/downward)
- **RADLEVEL_UP/DN(MAX_GEOMETRIES,0:MAXLAYERS)**: Level-by-level radiances

### Linearization Interfaces

#### Profile Jacobian Interface: TWOSTREAM_LPS_MASTER
```fortran
SUBROUTINE TWOSTREAM_LPS_MASTER ( &
    ! All TWOSTREAM_MASTER parameters plus:
    DO_PROFILEWFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
    N_TOTALPROFILE_WFS, &
    ! Additional outputs:
    PROFILEWF_TOA, PROFILEWF_BOA, &
    FLUXPROFILEWF_TOA, FLUXPROFILEWF_BOA &
)
```

**Purpose**: Compute ∂I/∂x for layer-wise parameters x ∈ {τ, ω, g}

#### Column Jacobian Interface: TWOSTREAM_LCS_MASTER
```fortran
SUBROUTINE TWOSTREAM_LCS_MASTER ( &
    ! All TWOSTREAM_MASTER parameters plus:
    DO_COLUMNWFS, N_TOTALCOLUMN_WFS, &
    ! Additional outputs:
    COLUMNWF_TOA, COLUMNWF_BOA, &
    FLUXCOLUMNWF_TOA, FLUXCOLUMNWF_BOA &
)
```

**Purpose**: Compute ∂I/∂X for column-integrated parameters X

---

## Module Reference

### Core Solution Module: 2stream_solutions.f90

#### Key Subroutines

##### TWOSTREAM_HOM_SOLUTION
```fortran
SUBROUTINE TWOSTREAM_HOM_SOLUTION ( &
    MAXLAYERS, N, FOURIER, STREAM_VALUE, PXSQ, &
    OMEGA, ASYMM, DELTAU_VERT, &
    SAB, DAB, EIGENVALUE, EIGENTRANS, &
    XPOS, XNEG, NORM_SAVED &
)
```
**Function**: Solves the homogeneous RTE to find eigenvalues and eigenvectors
**Mathematics**: 
- Computes eigenvalues λ from the characteristic equation
- Determines eigenvectors for upward/downward streams
- Applies normalization for numerical stability

##### TWOSTREAM_GBEAM_SOLUTION
```fortran
SUBROUTINE TWOSTREAM_GBEAM_SOLUTION ( &
    MAXLAYERS, N, FOURIER, PI4, STREAM_VALUE, &
    FLUX_FACTOR, DO_SOLAR_SOURCES, NBEAMS, &
    X0, OMEGA, ASYMM, DELTAU_VERT, &
    AVERAGE_SECANT, INITIAL_TRANS, &
    WUPPER, WLOWER &
)
```
**Function**: Computes particular integral for solar beam source
**Mathematics**: Green's function approach for inhomogeneous term

#### Mathematical Framework
The two-stream approximation reduces the RTE to a system:
```
μ dI⁺/dτ = I⁺ - ω/2[I⁺ + I⁻ + 2πδ(μ-μ₀)F₀e^(-τ/μ₀)]
-μ dI⁻/dτ = I⁻ - ω/2[I⁺ + I⁻ + 2πδ(μ-μ₀)F₀e^(-τ/μ₀)]
```

### Boundary Value Problem Module: 2stream_bvproblem.f90

#### Key Subroutines

##### TWOSTREAM_BVP_SOLUTION_PENTADIAG
```fortran
SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG ( &
    BVPINDEX, DO_PENTADIAG_INVERSE, &
    MAXLAYERS, MAXTOTAL, NSTR2, N_SUPDIAG, N_SUBDIAG, &
    NLAYERS, NTOTAL, NBEAMS, &
    SURFACE_FACTOR, ALBEDO, BRDF_F, DIRECT_BEAM, &
    DELTAU_VERT, L_DELTAU_VERT, SOLARBEAM_BOATRANS, &
    XPOS, XNEG, WUPPER, WLOWER, &
    BANDMAT2, SMAT2, IPIVOT, SIPIVOT, &
    COL2, SCOL2, LCON, MCON, STATUS, MESSAGE, TRACE &
)
```
**Function**: Solves the boundary value problem using pentadiagonal matrix methods
**Algorithm**: 
1. Constructs pentadiagonal matrix from boundary conditions
2. Applies LU decomposition with pivoting
3. Back-substitution to find integration constants

#### Boundary Conditions
- **Top boundary**: No downward diffuse radiation at TOA
- **Bottom boundary**: Surface reflection (Lambertian + BRDF)
- **Layer interfaces**: Continuity of radiance

### Intensity Calculation Module: 2stream_intensity.f90

#### Key Subroutines

##### TWOSTREAM_UPUSER_INTENSITY
```fortran
SUBROUTINE TWOSTREAM_UPUSER_INTENSITY ( &
    MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, &
    DO_THERMAL_TRANSONLY, &
    NLAYERS, N_USER_STREAMS, NBEAMS, &
    LAYER_PIS_CUTOFF, T_DELT_USERM, &
    GAMMA_M, GAMMA_P, SIGMA_P, &
    INITIAL_TRANS, ITRANS_USERM, T_DELT_MUBAR, &
    LCON_XVEC, MCON_XVEC, &
    WUPPER, WLOWER, LCON, MCON, &
    U_XPOS, U_XNEG, U_WPOS1, U_WNEG1, HMULT_1, HMULT_2, &
    EMULT_UP, LAYER_TSUP_UP, LAYER_TSUP_UTUP, &
    INTENSITY_F_UP, CUMSOURCE_UP &
)
```
**Function**: Computes upward radiance at user-defined angles
**Method**: Post-processing using method of particular solutions

##### TWOSTREAM_FLUXES
```fortran
SUBROUTINE TWOSTREAM_FLUXES ( &
    MAXBEAMS, NLAYERS, NBEAMS, &
    STREAM_VALUE, FLUX_FACTOR, &
    LCON, MCON, XPOS, XNEG, &
    FLUXES_TOA, FLUXES_BOA &
)
```
**Function**: Computes hemispherical fluxes
**Mathematics**: Integration over solid angle using quadrature

---

## Calling Patterns

### Standard Execution Flow

#### Phase 1: Initialization and Validation
```fortran
! 1. Dimension checking
CALL TWOSTREAM_CHECK_INPUT_DIMS(...)

! 2. Basic parameter validation  
CALL TWOSTREAM_CHECK_INPUTS_BASIC(...)

! 3. Optical property validation
CALL TWOSTREAM_CHECK_INPUTS_OPTICAL(...)
```

#### Phase 2: Geometric Preprocessing
```fortran
! 1. Chapman function calculation for spherical geometry
DO IB = 1, NBEAMS
    CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE(...)
ENDDO

! 2. User angle preprocessing
! 3. Delta-2S scaling application
! 4. Singularity busting for extreme optical properties
```

#### Phase 3: Optical Preprocessing  
```fortran
! 1. Quasi-spherical attenuation
CALL TWOSTREAM_QSPREP(...)

! 2. Transmittance factors
CALL TWOSTREAM_PREPTRANS(...)

! 3. Thermal setup (if enabled)
IF (DO_THERMAL_EMISSION) THEN
    CALL TWOSTREAM_THERMALSETUP(...)
ENDIF

! 4. Beam source multipliers
IF (DO_SOLAR_SOURCES .AND. DO_POSTPROCESSING) THEN
    CALL TWOSTREAM_EMULTMASTER(...)
ENDIF
```

#### Phase 4: Fourier Loop
```fortran
DO FOURIER = 0, N_FOURIERS
    ! Azimuth factor calculation
    
    ! Main computational engine
    CALL TWOSTREAM_FOURIER_MASTER(...)
    
    ! Convergence testing
    IF (DO_USER_OBSGEOMS) THEN
        CALL TWOSTREAM_CONVERGE_OBSGEO(...)
    ELSE
        CALL TWOSTREAM_CONVERGE(...)
    ENDIF
ENDDO
```

### FOURIER_MASTER Internal Flow

#### Solar Source Processing (m = 0)
```fortran
! 1. Direct beam attenuation
CALL TWOSTREAM_DIRECTBEAM(...)

! 2. Auxiliary geometry setup
CALL TWOSTREAM_AUXGEOM(...)

! 3. Homogeneous solutions
CALL TWOSTREAM_HOM_SOLUTION(...)
IF (DO_POSTPROCESSING) THEN
    CALL TWOSTREAM_HOM_USERSOLUTION(...)
ENDIF

! 4. Homogeneous multipliers
CALL TWOSTREAM_HMULT_MASTER(...)

! 5. BVP matrix setup
CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG(...)

! 6. Thermal processing (if enabled)
IF (DO_THERMAL_EMISSION) THEN
    CALL TWOSTREAM_THERMALGFSOLUTION(...)
    CALL TWOSTREAM_THERMALSTERMS(...)
ENDIF

! 7. BVP solution
CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG(...)

! 8. Radiance calculation
IF (DO_UPWELLING .AND. DO_POSTPROCESSING) THEN
    CALL TWOSTREAM_UPUSER_INTENSITY(...)
ENDIF
IF (DO_DNWELLING .AND. DO_POSTPROCESSING) THEN  
    CALL TWOSTREAM_DNUSER_INTENSITY(...)
ENDIF

! 9. Flux calculation
IF (.NOT. DO_MVOUT_ONLY) THEN
    CALL TWOSTREAM_FLUXES(...)
ENDIF
```

#### Higher Fourier Components (m ≥ 1)
```fortran
! Solar beam particular solution
CALL TWOSTREAM_GBEAM_SOLUTION(...)

! Repeat BVP solution and intensity calculation
! (Thermal terms not needed for m ≥ 1)
```

### Linearization Call Patterns

#### Profile Jacobian Flow
```fortran
! Standard calculation first
CALL TWOSTREAM_FOURIER_MASTER(...)

! Then linearized calculation
IF (DO_PROFILEWFS) THEN
    ! Profile-specific preprocessing
    CALL TWOSTREAM_LP_QSPREP(...)
    CALL TWOSTREAM_LP_EMULTMASTER(...)
    
    ! Linearized homogeneous solutions
    CALL TWOSTREAM_L_HOM_SOLUTION(...)
    CALL TWOSTREAM_L_HOM_USERSOLUTION(...)
    
    ! Linearized BVP solution
    CALL TWOSTREAM_BVP_LP_SOLUTION_MASTER(...)
    
    ! Linearized intensity calculation
    CALL TWOSTREAM_UPUSER_PROFILEWF(...)
    CALL TWOSTREAM_DNUSER_PROFILEWF(...)
    CALL TWOSTREAM_FLUXES_PROFILEWF(...)
ENDIF
```

---

## Mathematical Framework

### Two-Stream Approximation

#### Radiative Transfer Equation
The azimuthally-averaged intensity I(τ,μ) satisfies:
```
μ ∂I(τ,μ)/∂τ = I(τ,μ) - S(τ,μ)
```

where the source function is:
```
S(τ,μ) = ω(τ)/4π ∫_{-1}^{1} P(τ,μ,μ') I(τ,μ') dμ' + ω(τ)/4 P(τ,μ,μ₀) F₀ e^{-τ/μ₀}
```

#### Two-Stream Discretization
Approximate with two streams: I⁺ (upward) and I⁻ (downward)
```
I⁺ ≈ I(τ,+μ)  (μ > 0)
I⁻ ≈ I(τ,-μ)  (μ < 0)
```

#### Discrete Ordinates System
```
+μ dI⁺/dτ = I⁺ - ω/2[(1+g)I⁺ + (1-g)I⁻] - ωF₀/(4μ₀)P₁₁e^{-τ/μ₀}
-μ dI⁻/dτ = I⁻ - ω/2[(1-g)I⁺ + (1+g)I⁻] - ωF₀/(4μ₀)P₁₋₁e^{-τ/μ₀}
```

where g is the asymmetry parameter and P₁ₘ are phase function elements.

### Homogeneous Solution

#### Eigenvalue Problem
The homogeneous system yields the characteristic equation:
```
λ² = (1-ω)/μ² + ω(1-g²)/(4μ²)
```

#### Eigenvectors
```
k⁺ = √λ,  k⁻ = -√λ
e⁺ = [1, γ],  e⁻ = [γ, 1]
```

where γ = (1-ω/μ + √λ)/(ω(1-g)/(2μ))

### Particular Solution

#### Green's Function Method
For solar source F₀δ(μ-μ₀), the particular solution is:
```
I_p(τ) = G(τ,μ₀) F₀ e^{-τ/μ₀}
```

where G(τ,μ₀) satisfies the inhomogeneous equation.

### Boundary Value Problem

#### Matrix Formulation
The system of boundary conditions yields:
```
[A] {x} = {b}
```

where {x} contains integration constants and {b} contains source terms.

#### Pentadiagonal Structure
The matrix [A] has pentadiagonal structure due to layer-coupling:
```
[A] = [
  d₁  e₁  f₁   0   0  ...
  c₂  d₂  e₂  f₂   0  ...
  b₃  c₃  d₃  e₃  f₃  ...
   0  b₄  c₄  d₄  e₄  ...
  ...
]
```

### Linearization Theory

#### Profile Jacobian
For parameter p in layer k:
```
∂I/∂p = ∂I/∂{LCON,MCON} × ∂{LCON,MCON}/∂p
```

#### Chain Rule Application
```
∂{LCON,MCON}/∂p = -[A]⁻¹ × ∂[A]/∂p × {LCON,MCON} + [A]⁻¹ × ∂{b}/∂p
```

---

## Implementation Guide

### Basic Setup

#### Required Parameters
```fortran
! Dimensions
integer, parameter :: MAXLAYERS = 114
integer, parameter :: MAXBEAMS = 1  
integer, parameter :: MAX_USER_STREAMS = 1
integer, parameter :: MAX_USER_RELAZMS = 1
integer, parameter :: MAX_USER_OBSGEOMS = 1
integer, parameter :: MAXTOTAL = 2 * MAXLAYERS
integer, parameter :: MAX_GEOMETRIES = MAXBEAMS * MAX_USER_STREAMS * MAX_USER_RELAZMS

! Control flags  
logical :: DO_UPWELLING = .TRUE.
logical :: DO_DNWELLING = .TRUE.
logical :: DO_SOLAR_SOURCES = .TRUE.
logical :: DO_THERMAL_EMISSION = .FALSE.
logical :: DO_BRDF_SURFACE = .FALSE.
logical :: DO_USER_OBSGEOMS = .TRUE.
```

#### Geometry Configuration
```fortran
! Layer structure
integer :: NLAYERS = 114
real(dp) :: HEIGHT_GRID(0:MAXLAYERS)
real(dp) :: EARTH_RADIUS = 6371.0_dp

! Solar geometry
integer :: NBEAMS = 1
real(dp) :: BEAM_SZAS(MAXBEAMS) = [30.0_dp]  ! Solar zenith angle

! Viewing geometry  
integer :: N_USER_OBSGEOMS = 1
real(dp) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3)
USER_OBSGEOMS(1,:) = [30.0_dp, 0.0_dp, 0.0_dp]  ! SZA, VZA, AZM
```

#### Optical Properties
```fortran
! Atmospheric optical properties
real(dp) :: DELTAU_INPUT(MAXLAYERS) = 0.01_dp    ! Optical depth per layer
real(dp) :: OMEGA_INPUT(MAXLAYERS) = 0.5_dp      ! Single scattering albedo
real(dp) :: ASYMM_INPUT(MAXLAYERS) = 0.9_dp      ! Asymmetry parameter

! Surface properties
real(dp) :: LAMBERTIAN_ALBEDO = 0.3_dp
```

#### Numerical Control
```fortran
! BVP solver configuration
integer :: BVPINDEX = 1                    ! Use pentadiagonal solver
logical :: DO_PENTADIAG_INVERSE = .FALSE.  ! Bottom-to-top solution
real(dp) :: BVPSCALEFACTOR = 1.0_dp       ! Scaling factor

! Taylor series parameters
integer :: TAYLOR_ORDER = 2
real(dp) :: TAYLOR_SMALL = 1.0E-4_dp
real(dp) :: TCUTOFF = 1.0E-8_dp           ! Optical depth cutoff

! Stream value
real(dp) :: STREAM_VALUE = 0.5_dp         ! Gaussian quadrature point
```

### Advanced Configuration

#### BRDF Surface Setup
```fortran
logical :: DO_BRDF_SURFACE = .TRUE.

! BRDF Fourier components (typically from separate BRDF supplement)
real(dp) :: BRDF_F_0(0:1,MAXBEAMS)       ! Solar incident
real(dp) :: BRDF_F(0:1)                  ! Quadrature incident  
real(dp) :: UBRDF_F(0:1,MAX_USER_STREAMS) ! User streams

! Example: Lambertian + Rahman kernel combination
BRDF_F_0(0,:) = 0.1_dp   ! Isotropic component
BRDF_F_0(1,:) = 0.05_dp  ! First Fourier component
```

#### Thermal Emission Setup  
```fortran
logical :: DO_THERMAL_EMISSION = .TRUE.
logical :: DO_SURFACE_EMISSION = .TRUE.

! Thermal black-body inputs (Kelvin)
real(dp) :: THERMAL_BB_INPUT(0:MAXLAYERS)
THERMAL_BB_INPUT = 250.0_dp  ! Typical atmospheric temperature

! Surface thermal properties
real(dp) :: SURFBB = 290.0_dp      ! Surface temperature
real(dp) :: EMISSIVITY = 0.95_dp   ! Surface emissivity
```

#### Linearization Setup
```fortran
! Profile Jacobian configuration
logical :: DO_PROFILEWFS = .TRUE.
logical :: LAYER_VARY_FLAG(MAXLAYERS) = .TRUE.  ! Vary all layers
integer :: LAYER_VARY_NUMBER(MAXLAYERS) = 3     ! τ, ω, g for each layer
integer :: N_TOTALPROFILE_WFS = NLAYERS * 3

! Column Jacobian configuration  
logical :: DO_COLUMNWFS = .TRUE.
integer :: N_TOTALCOLUMN_WFS = 3  ! Total τ, ω̄, ḡ
```

### Error Handling

#### Input Validation
```fortran
! Check return status
if (STATUS_INPUTCHECK /= 0) then
    print*, 'Input validation failed:'
    do i = 1, C_NMESSAGES
        print*, C_MESSAGES(i)
        print*, 'Action:', C_ACTIONS(i)
    enddo
    stop
endif
```

#### Execution Monitoring
```fortran
! Check execution status
if (STATUS_EXECUTION /= 0) then
    print*, 'Execution failed:'
    print*, 'Error message:', trim(E_MESSAGE)
    print*, 'Trace 1:', trim(E_TRACE_1)  
    print*, 'Trace 2:', trim(E_TRACE_2)
    stop
endif
```

---

## Performance Optimization

### Numerical Solver Selection

#### BVP Solver Comparison
| Solver | BVPINDEX | Performance | Accuracy | Memory |
|--------|----------|-------------|----------|---------|
| LAPACK | 0 | Standard | High | Standard |
| Pentadiagonal #1 | 1 | 3-5x faster | High | Low |
| Pentadiagonal #2 | 2 | 3-5x faster | High | Low |

**Recommendation**: Use BVPINDEX = 1 for optimal performance

#### Pentadiagonal Advantages
- **Banded structure**: Exploits sparse matrix properties
- **O(N) complexity**: Linear scaling with layer number
- **Cache efficiency**: Better memory access patterns
- **Numerical stability**: Specialized algorithms for banded systems

### Geometry Optimization

#### Observational Geometry Mode
```fortran
logical :: DO_USER_OBSGEOMS = .TRUE.
```
**Benefits**:
- Eliminates redundant azimuth calculations
- Reduces memory footprint
- Simplifies post-processing loops
- 40-60% speed improvement for single-geometry calculations

#### Post-processing Control
```fortran
logical :: DO_MVOUT_ONLY = .TRUE.    ! Flux-only output
logical :: DO_POSTPROCESSING = .FALSE. ! Skip user-angle calculations
```

### Memory Management

#### Compiler Optimizations
```bash
# Intel Fortran
ifort -O3 -heap-arrays -mcmodel=large -fpp

# GNU Fortran  
gfortran -O3 -mcmodel=large -ffree-form -ffree-line-length-none
```

#### Array Dimensioning
- Use exact dimensions when possible
- Minimize MAX_GEOMETRIES for single-angle calculations
- Consider memory layout for cache efficiency

### Taylor Series Optimization

#### Optically Thin Layer Handling
```fortran
integer :: TAYLOR_ORDER = 3           ! Higher order for accuracy
real(dp) :: TAYLOR_SMALL = 1.0E-5_dp  ! Smaller threshold for stability
real(dp) :: TCUTOFF = 1.0E-9_dp      ! Skip extremely thin layers
```

**Benefits**:
- Prevents numerical overflow in exponential calculations
- Maintains accuracy for tau << 1 cases
- Reduces computational overhead for negligible layers

---

## Error Handling

### Common Error Scenarios

#### Input Validation Errors
```fortran
! Typical error messages and solutions:

! "NLAYERS > MAXLAYERS"
! Solution: Increase MAXLAYERS or reduce NLAYERS

! "Single scattering albedo out of range"  
! Solution: Ensure 0 ≤ OMEGA_INPUT ≤ 1

! "Asymmetry parameter out of range"
! Solution: Ensure -1 ≤ ASYMM_INPUT ≤ 1

! "Solar zenith angle too large"
! Solution: Ensure BEAM_SZAS < 89°
```

#### Numerical Stability Issues
```fortran
! Singularity busting activates automatically:
! - OMEGA clamped to [1e-9, 0.999999999]
! - ASYMM clamped to [-0.999999999, 0.999999999]
! - Small values set to ±1e-9 to avoid underflow
```

#### Convergence Problems
```fortran
! Fourier series convergence issues:
! - Increase FOURIER_TOL (in convergence routines)
! - Check for extreme optical properties
! - Verify surface BRDF parameters
```

### Debugging Strategies

#### Debug Output Activation
```fortran
logical :: DO_DEBUG_INPUT = .TRUE.  ! In master routines
```

#### Write Module Usage
```fortran
! Enable detailed parameter dumps
CALL TWOSTREAM_WRITE_STD_INPUT(...)
CALL TWOSTREAM_WRITE_SUP_BRDF_INPUT(...)
```

#### Timing Analysis
```fortran
real(sp) :: geom_timer
! Timer incremented during geometry calculations
! Use for performance profiling
```

---

## Examples

### Example 1: Basic Solar Calculation

```fortran
program basic_solar
    use twostream_master_m
    implicit none
    
    ! Parameters
    integer, parameter :: dp = selected_real_kind(15)
    integer, parameter :: MAXLAYERS = 50
    integer, parameter :: MAXBEAMS = 1
    integer, parameter :: MAX_USER_STREAMS = 1
    integer, parameter :: MAX_USER_RELAZMS = 1
    integer, parameter :: MAX_USER_OBSGEOMS = 1
    integer, parameter :: MAXTOTAL = 2 * MAXLAYERS
    integer, parameter :: MAX_GEOMETRIES = 1
    integer, parameter :: MAXMESSAGES = 25
    
    ! Control variables
    logical :: DO_UPWELLING = .TRUE.
    logical :: DO_DNWELLING = .TRUE.
    logical :: DO_PLANE_PARALLEL = .FALSE.
    logical :: DO_2S_LEVELOUT = .FALSE.
    logical :: DO_MVOUT_ONLY = .FALSE.
    logical :: DO_ADDITIONAL_MVOUT = .FALSE.
    logical :: DO_SOLAR_SOURCES = .TRUE.
    logical :: DO_THERMAL_EMISSION = .FALSE.
    logical :: DO_SURFACE_EMISSION = .FALSE.
    logical :: DO_D2S_SCALING = .FALSE.
    logical :: DO_BRDF_SURFACE = .FALSE.
    logical :: DO_USER_OBSGEOMS = .TRUE.
    logical :: DO_SURFACE_LEAVING = .FALSE.
    logical :: DO_SL_ISOTROPIC = .FALSE.
    logical :: DO_PENTADIAG_INVERSE = .FALSE.
    
    ! Numerical control
    integer :: BVPINDEX = 1
    real(dp) :: BVPSCALEFACTOR = 1.0_dp
    integer :: TAYLOR_ORDER = 2
    real(dp) :: TAYLOR_SMALL = 1.0E-4_dp
    real(dp) :: TCUTOFF = 1.0E-8_dp
    
    ! Geometry
    integer :: NLAYERS = 20
    integer :: NTOTAL = 2 * NLAYERS
    real(dp) :: STREAM_VALUE = 0.5_dp
    integer :: N_USER_OBSGEOMS = 1
    real(dp) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3)
    integer :: N_USER_STREAMS = 1
    real(dp) :: USER_ANGLES(MAX_USER_STREAMS)
    integer :: N_USER_RELAZMS = 1
    real(dp) :: USER_RELAZMS(MAX_USER_RELAZMS)
    real(dp) :: FLUX_FACTOR = 1.0_dp
    integer :: NBEAMS = 1
    real(dp) :: BEAM_SZAS(MAXBEAMS)
    real(dp) :: EARTH_RADIUS = 6371.0_dp
    real(dp) :: HEIGHT_GRID(0:MAXLAYERS)
    
    ! Optical properties
    real(dp) :: DELTAU_INPUT(MAXLAYERS)
    real(dp) :: OMEGA_INPUT(MAXLAYERS)
    real(dp) :: ASYMM_INPUT(MAXLAYERS)
    real(dp) :: D2S_SCALING(MAXLAYERS)
    real(dp) :: THERMAL_BB_INPUT(0:MAXLAYERS)
    
    ! Surface properties
    real(dp) :: LAMBERTIAN_ALBEDO = 0.3_dp
    real(dp) :: BRDF_F_0(0:1,MAXBEAMS)
    real(dp) :: BRDF_F(0:1)
    real(dp) :: UBRDF_F(0:1,MAX_USER_STREAMS)
    real(dp) :: EMISSIVITY = 0.95_dp
    real(dp) :: SURFBB = 290.0_dp
    real(dp) :: SLTERM_ISOTROPIC(MAXBEAMS)
    real(dp) :: SLTERM_F_0(0:1,MAXBEAMS)
    
    ! Outputs
    real(dp) :: INTENSITY_TOA(MAX_GEOMETRIES)
    real(dp) :: INTENSITY_BOA(MAX_GEOMETRIES)
    real(dp) :: FLUXES_TOA(MAXBEAMS,2)
    real(dp) :: FLUXES_BOA(MAXBEAMS,2)
    real(dp) :: RADLEVEL_UP(MAX_GEOMETRIES,0:MAXLAYERS)
    real(dp) :: RADLEVEL_DN(MAX_GEOMETRIES,0:MAXLAYERS)
    integer :: N_GEOMETRIES
    
    ! Error handling
    integer :: STATUS_INPUTCHECK, STATUS_EXECUTION
    integer :: C_NMESSAGES
    character*100 :: C_MESSAGES(0:MAXMESSAGES)
    character*100 :: C_ACTIONS(0:MAXMESSAGES)
    character*100 :: E_MESSAGE, E_TRACE_1, E_TRACE_2
    
    ! Timing
    real :: geom_timer = 0.0
    
    ! Local variables
    integer :: i
    
    ! Setup geometry
    USER_OBSGEOMS(1,:) = [30.0_dp, 0.0_dp, 0.0_dp]  ! SZA=30°, VZA=0°, AZM=0°
    BEAM_SZAS(1) = 30.0_dp
    
    ! Setup height grid (1 km per layer)
    do i = 0, NLAYERS
        HEIGHT_GRID(i) = real(NLAYERS - i, dp)  ! km
    enddo
    
    ! Setup optical properties (uniform atmosphere)
    DELTAU_INPUT = 0.05_dp    ! Optical depth per layer
    OMEGA_INPUT = 0.8_dp      ! Single scattering albedo
    ASYMM_INPUT = 0.7_dp      ! Asymmetry parameter
    D2S_SCALING = 0.0_dp      ! No delta scaling
    
    ! Initialize arrays
    BRDF_F_0 = 0.0_dp
    BRDF_F = 0.0_dp
    UBRDF_F = 0.0_dp
    SLTERM_ISOTROPIC = 0.0_dp
    SLTERM_F_0 = 0.0_dp
    THERMAL_BB_INPUT = 250.0_dp
    
    ! Call 2STREAM
    CALL TWOSTREAM_MASTER ( &
        MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES, &
        MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS, &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT, &
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, &
        DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS, &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE, &
        BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF, &
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS, &
        N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS, &
        FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID, &
        DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING, &
        THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, &
        EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0, &
        INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA, &
        RADLEVEL_UP, RADLEVEL_DN, N_GEOMETRIES, &
        STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS, &
        STATUS_EXECUTION, E_MESSAGE, E_TRACE_1, E_TRACE_2, &
        geom_timer &
    )
    
    ! Check for errors
    if (STATUS_INPUTCHECK /= 0) then
        print*, 'Input check failed'
        stop
    endif
    
    if (STATUS_EXECUTION /= 0) then
        print*, 'Execution failed:', trim(E_MESSAGE)
        stop
    endif
    
    ! Print results
    print*, 'TOA Radiance:', INTENSITY_TOA(1)
    print*, 'BOA Radiance:', INTENSITY_BOA(1)
    print*, 'TOA Flux (up/down):', FLUXES_TOA(1,:)
    print*, 'BOA Flux (up/down):', FLUXES_BOA(1,:)
    print*, 'Geometry timer (s):', geom_timer
    
end program basic_solar
```

### Example 2: Profile Jacobian Calculation

```fortran
program profile_jacobian
    use twostream_lps_master_m
    implicit none
    
    ! [Include all declarations from Example 1]
    
    ! Additional linearization variables
    logical :: DO_PROFILEWFS = .TRUE.
    logical :: LAYER_VARY_FLAG(MAXLAYERS)
    integer :: LAYER_VARY_NUMBER(MAXLAYERS)
    integer :: N_TOTALPROFILE_WFS
    
    ! Jacobian outputs
    real(dp) :: PROFILEWF_TOA(MAX_GEOMETRIES,MAXLAYERS,3)
    real(dp) :: PROFILEWF_BOA(MAX_GEOMETRIES,MAXLAYERS,3)
    real(dp) :: FLUXPROFILEWF_TOA(MAXBEAMS,2,MAXLAYERS,3)
    real(dp) :: FLUXPROFILEWF_BOA(MAXBEAMS,2,MAXLAYERS,3)
    
    ! [Include all setup from Example 1]
    
    ! Setup linearization
    LAYER_VARY_FLAG = .TRUE.        ! Vary all layers
    LAYER_VARY_NUMBER = 3           ! τ, ω, g for each layer
    N_TOTALPROFILE_WFS = NLAYERS * 3
    
    ! Call linearized 2STREAM
    CALL TWOSTREAM_LPS_MASTER ( &
        ! [All parameters from TWOSTREAM_MASTER call]
        DO_PROFILEWFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        N_TOTALPROFILE_WFS, &
        PROFILEWF_TOA, PROFILEWF_BOA, &
        FLUXPROFILEWF_TOA, FLUXPROFILEWF_BOA &
    )
    
    ! Print Jacobian results
    print*, 'dI_TOA/d_tau(layer1):', PROFILEWF_TOA(1,1,1)
    print*, 'dI_TOA/d_omega(layer1):', PROFILEWF_TOA(1,1,2)
    print*, 'dI_TOA/d_g(layer1):', PROFILEWF_TOA(1,1,3)
    
end program profile_jacobian
```

---

## Conclusion

The 2STREAM v2.4 library provides a comprehensive and efficient solution for atmospheric radiative transfer modeling. Its modular design, robust numerical methods, and extensive linearization capabilities make it well-suited for a wide range of scientific applications.

**Key Strengths**:
- High computational efficiency through optimized algorithms
- Comprehensive physical modeling capabilities
- Robust error handling and numerical stability
- Extensive linearization support for inverse problems
- Well-documented and maintainable code structure

**Recommended Usage Patterns**:
- Use pentadiagonal solver (BVPINDEX=1) for optimal performance
- Enable observational geometry mode for single-angle calculations
- Apply Taylor series for optically thin layers
- Implement proper error checking and validation
- Consider linearization capabilities for parameter retrieval applications

This documentation provides the foundation for effective utilization of the 2STREAM v2.4 library in research and operational applications.