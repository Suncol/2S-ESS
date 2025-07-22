# 2S-ESS API Reference

## Quick Navigation
- [Core Interfaces](#core-interfaces)
- [Two-Stream API](#two-stream-api)
- [First Order ESS API](#first-order-ess-api)
- [Parameter Reference](#parameter-reference)
- [Error Handling](#error-handling)
- [Usage Examples](#usage-examples)

---

## Core Interfaces

### TWOSTREAM_MASTER
**Primary interface for two-stream radiative transfer calculations**

**Location**: `sourcecode/2stream_2p4_feb15/2stream_master.f90`

#### Signature
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
| Parameter | Type | Purpose |
|-----------|------|---------|
| **MAXLAYERS** | integer | Maximum atmospheric layers |
| **DO_SOLAR_SOURCES** | logical | Enable solar radiation |
| **DO_THERMAL_EMISSION** | logical | Enable thermal emission |
| **DELTAU_INPUT** | real(dp) array | Layer optical depths |
| **OMEGA_INPUT** | real(dp) array | Single scattering albedos |
| **INTENSITY_TOA** | real(dp) array | Top-of-atmosphere radiances |

---

## Two-Stream API

### Linearization Interfaces

#### TWOSTREAM_LPS_MASTER
**Profile linearization (layer-by-layer Jacobians)**

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

**Purpose**: Computes ∂I/∂x for x ∈ {τ, ω, g} in each layer

#### TWOSTREAM_LCS_MASTER
**Column linearization (integrated Jacobians)**

```fortran
SUBROUTINE TWOSTREAM_LCS_MASTER ( &
    ! All TWOSTREAM_MASTER parameters plus:
    DO_COLUMNWFS, N_TOTALCOLUMN_WFS, &
    ! Additional outputs:
    COLUMNWF_TOA, COLUMNWF_BOA, &
    FLUXCOLUMNWF_TOA, FLUXCOLUMNWF_BOA &
)
```

**Purpose**: Computes ∂I/∂X for column-integrated parameters X

### Core Solution Modules

#### TWOSTREAM_HOM_SOLUTION
**Homogeneous solution calculation**

**Location**: `sourcecode/2stream_2p4_feb15/2stream_solutions.f90`

```fortran
SUBROUTINE TWOSTREAM_HOM_SOLUTION ( &
    MAXLAYERS, N, FOURIER, STREAM_VALUE, PXSQ, &
    OMEGA, ASYMM, DELTAU_VERT, &
    SAB, DAB, EIGENVALUE, EIGENTRANS, &
    XPOS, XNEG, NORM_SAVED &
)
```

**Purpose**: Solves eigenvalue problem for homogeneous RTE

#### TWOSTREAM_BVP_SOLUTION_PENTADIAG
**Boundary value problem solver**

**Location**: `sourcecode/2stream_2p4_feb15/2stream_bvproblem.f90`

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

**Purpose**: Solves pentadiagonal matrix system for integration constants

---

## First Order ESS API

### Scalar Single Scattering

#### FO_ScalarSS_RTCalcs_I
**Single scattering radiance calculations**

**Location**: `sourcecode/FO_1p4_sourcecode/FO_ScalarSS_RTCalcs_I.f90`

```fortran
SUBROUTINE FO_ScalarSS_RTCalcs_I ( &
    ! Geometry and optical inputs
    ! Single scattering calculations
    ! Radiance outputs
)
```

### Thermal Single Scattering

#### FO_Thermal_RTCalcs_I
**Thermal single scattering calculations**

**Location**: `sourcecode/FO_1p4_sourcecode/FO_Thermal_RTCalcs_I.f90`

```fortran
SUBROUTINE FO_Thermal_RTCalcs_I ( &
    ! Thermal source inputs
    ! Temperature profiles
    ! Thermal radiance outputs
)
```

---

## Parameter Reference

### Dimensional Parameters
```fortran
integer, parameter :: MAXLAYERS = 114          ! Maximum atmospheric layers
integer, parameter :: MAXTOTAL = 2 * MAXLAYERS ! BVP matrix dimension
integer, parameter :: MAXBEAMS = 1             ! Solar beam count
integer, parameter :: MAX_USER_STREAMS = 1     ! User viewing angles
integer, parameter :: MAX_USER_OBSGEOMS = 1    ! Observation geometries
```

### Control Flags
```fortran
logical :: DO_UPWELLING = .TRUE.        ! Calculate upward radiance
logical :: DO_DNWELLING = .TRUE.        ! Calculate downward radiance
logical :: DO_SOLAR_SOURCES = .TRUE.    ! Include solar radiation
logical :: DO_THERMAL_EMISSION = .FALSE. ! Include thermal emission
logical :: DO_BRDF_SURFACE = .FALSE.    ! Use BRDF surface model
logical :: DO_USER_OBSGEOMS = .TRUE.    ! Observational geometry mode
```

### Numerical Control
```fortran
integer :: BVPINDEX = 1                  ! Pentadiagonal solver
integer :: TAYLOR_ORDER = 2             ! Taylor series order
real(dp) :: TAYLOR_SMALL = 1.0E-4_dp    ! Taylor expansion threshold
real(dp) :: TCUTOFF = 1.0E-8_dp        ! Optical depth cutoff
real(dp) :: STREAM_VALUE = 0.5_dp      ! Gaussian quadrature point
```

### Optical Properties
```fortran
real(dp) :: DELTAU_INPUT(MAXLAYERS)     ! Layer optical depths
real(dp) :: OMEGA_INPUT(MAXLAYERS)      ! Single scattering albedos
real(dp) :: ASYMM_INPUT(MAXLAYERS)      ! Asymmetry parameters (g)
```

### Geometry Setup
```fortran
real(dp) :: BEAM_SZAS(MAXBEAMS)         ! Solar zenith angles
real(dp) :: USER_OBSGEOMS(MAX_OBSGEOMS,3) ! [SZA, VZA, AZM]
real(dp) :: HEIGHT_GRID(0:MAXLAYERS)    ! Layer height boundaries
real(dp) :: EARTH_RADIUS = 6371.0_dp    ! Earth radius (km)
```

---

## Error Handling

### Input Validation
```fortran
integer :: STATUS_INPUTCHECK             ! Input validation status
integer :: C_NMESSAGES                   ! Number of messages
character*100 :: C_MESSAGES(0:MAXMESSAGES) ! Error messages
character*100 :: C_ACTIONS(0:MAXMESSAGES)  ! Recommended actions
```

### Execution Monitoring
```fortran
integer :: STATUS_EXECUTION             ! Execution status
character*100 :: E_MESSAGE              ! Error message
character*100 :: E_TRACE_1              ! Error trace level 1
character*100 :: E_TRACE_2              ! Error trace level 2
```

### Common Error Checks
```fortran
! Check input validation
if (STATUS_INPUTCHECK /= 0) then
    print*, 'Input validation failed'
    do i = 1, C_NMESSAGES
        print*, trim(C_MESSAGES(i))
        print*, 'Action:', trim(C_ACTIONS(i))
    enddo
    stop
endif

! Check execution status
if (STATUS_EXECUTION /= 0) then
    print*, 'Execution failed:', trim(E_MESSAGE)
    print*, 'Trace:', trim(E_TRACE_1)
    stop
endif
```

---

## Usage Examples

### Basic Solar Calculation
```fortran
program basic_example
    implicit none
    
    ! Declare all required parameters
    ! [Parameter declarations...]
    
    ! Setup geometry
    USER_OBSGEOMS(1,:) = [30.0_dp, 0.0_dp, 0.0_dp]  ! SZA=30°, VZA=0°
    
    ! Setup optical properties
    DELTAU_INPUT = 0.01_dp    ! Uniform optical depth
    OMEGA_INPUT = 0.5_dp      ! Single scattering albedo
    ASYMM_INPUT = 0.9_dp      ! Asymmetry parameter
    
    ! Surface properties
    LAMBERTIAN_ALBEDO = 0.3_dp
    
    ! Call main interface
    CALL TWOSTREAM_MASTER( [all parameters] )
    
    ! Check errors and process results
    if (STATUS_EXECUTION == 0) then
        print*, 'TOA Radiance:', INTENSITY_TOA(1)
        print*, 'BOA Radiance:', INTENSITY_BOA(1)
    endif
    
end program
```

### Profile Jacobian Example
```fortran
program jacobian_example
    ! Setup linearization
    DO_PROFILEWFS = .TRUE.
    LAYER_VARY_FLAG = .TRUE.        ! All layers
    LAYER_VARY_NUMBER = 3           ! τ, ω, g
    N_TOTALPROFILE_WFS = NLAYERS * 3
    
    ! Call linearized interface
    CALL TWOSTREAM_LPS_MASTER( [parameters including jacobian setup] )
    
    ! Access jacobians
    print*, 'dI/dτ (layer 1):', PROFILEWF_TOA(1,1,1)
    print*, 'dI/dω (layer 1):', PROFILEWF_TOA(1,1,2)
    print*, 'dI/dg (layer 1):', PROFILEWF_TOA(1,1,3)
end program
```

---

## Performance Notes

### Optimization Recommendations
1. **Use pentadiagonal solver**: `BVPINDEX = 1` (3-5x faster)
2. **Enable observational geometry**: `DO_USER_OBSGEOMS = .TRUE.`
3. **Optimize for single angles**: Minimize `MAX_GEOMETRIES`
4. **Compiler flags**: Use `-O3` optimization

### Memory Considerations
- Regular version: Standard memory footprint
- Optimized version: Reduced memory usage with identical accuracy
- Use `mcmodel=large` for extensive atmospheric datasets

*Complete API documentation with all parameters and detailed descriptions available in `/docs/2STREAM_v2.4_Technical_Documentation.md`*