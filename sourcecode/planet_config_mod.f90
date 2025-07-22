module planet_config_mod
  implicit none
  private
  
  integer, parameter :: dp = kind(1.0d0)
  
  ! Planet configuration type
  type :: planet_config_type
    ! Planetary geometry
    real(dp) :: radius_km = 6371.0_dp
    real(dp) :: radius_min_km = 6320.0_dp
    real(dp) :: radius_max_km = 6420.0_dp
    
    ! Ocean parameters (if applicable)
    logical  :: has_oceans = .true.
    real(dp) :: salinity_reference_ppt = 34.3_dp
    real(dp) :: refractive_correction = 0.006_dp
    
    ! Planet identification
    character(len=20) :: name = 'Earth'
  end type planet_config_type
  
  ! Global configuration instance
  type(planet_config_type), save :: current_planet_config
  
  ! Public interface
  public :: planet_config_type, current_planet_config
  public :: set_planet_config, get_planet_config
  public :: set_earth_config, set_mars_config, set_venus_config

contains

  subroutine set_planet_config(config)
    type(planet_config_type), intent(in) :: config
    current_planet_config = config
  end subroutine set_planet_config
  
  function get_planet_config() result(config)
    type(planet_config_type) :: config
    config = current_planet_config
  end function get_planet_config
  
  subroutine set_earth_config()
    current_planet_config%name = 'Earth'
    current_planet_config%radius_km = 6371.0_dp
    current_planet_config%radius_min_km = 6320.0_dp
    current_planet_config%radius_max_km = 6420.0_dp
    current_planet_config%has_oceans = .true.
    current_planet_config%salinity_reference_ppt = 34.3_dp
    current_planet_config%refractive_correction = 0.006_dp
  end subroutine set_earth_config
  
  subroutine set_mars_config()
    current_planet_config%name = 'Mars'
    current_planet_config%radius_km = 3389.5_dp
    current_planet_config%radius_min_km = 3380.0_dp
    current_planet_config%radius_max_km = 3400.0_dp
    current_planet_config%has_oceans = .false.
  end subroutine set_mars_config
  
  subroutine set_venus_config()
    current_planet_config%name = 'Venus'
    current_planet_config%radius_km = 6051.8_dp
    current_planet_config%radius_min_km = 6040.0_dp
    current_planet_config%radius_max_km = 6060.0_dp
    current_planet_config%has_oceans = .false.
  end subroutine set_venus_config

end module planet_config_mod