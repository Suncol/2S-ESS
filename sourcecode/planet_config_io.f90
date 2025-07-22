module planet_config_io
  use planet_config_mod
  implicit none
  
  private
  public :: load_planet_config_from_file
  
  ! Import dp parameter from planet_config_mod
  integer, parameter :: dp = kind(1.0d0)
  
contains

  subroutine load_planet_config_from_file(filename)
    character(len=*), intent(in) :: filename
    
    ! Namelist variables
    character(len=20) :: planet_name = 'Earth'
    real(dp) :: planet_radius_km = 6371.0_dp
    real(dp) :: planet_radius_min_km = 6320.0_dp
    real(dp) :: planet_radius_max_km = 6420.0_dp
    logical  :: has_oceans = .true.
    real(dp) :: salinity_reference_ppt = 34.3_dp
    real(dp) :: refractive_correction = 0.006_dp
    
    namelist /PLANET_CONFIG/ planet_name, planet_radius_km, &
                            planet_radius_min_km, planet_radius_max_km, &
                            has_oceans, salinity_reference_ppt, &
                            refractive_correction
    
    type(planet_config_type) :: config
    integer :: unit_num, ios
    
    ! Read configuration
    open(newunit=unit_num, file=filename, status='old', iostat=ios)
    if (ios == 0) then
      read(unit_num, nml=PLANET_CONFIG, iostat=ios)
      close(unit_num)
      
      if (ios == 0) then
        ! Set configuration
        config%name = planet_name
        config%radius_km = planet_radius_km
        config%radius_min_km = planet_radius_min_km
        config%radius_max_km = planet_radius_max_km
        config%has_oceans = has_oceans
        config%salinity_reference_ppt = salinity_reference_ppt
        config%refractive_correction = refractive_correction
        
        call set_planet_config(config)
        print *, 'Loaded planetary configuration for: ', trim(planet_name)
      else
        print *, 'Warning: Error reading planet configuration from ', trim(filename)
        print *, 'Using default Earth configuration'
        call set_earth_config()
      endif
    else
      print *, 'Warning: Could not open planet configuration file: ', trim(filename)
      print *, 'Using default Earth configuration'
      call set_earth_config()
    endif
  end subroutine load_planet_config_from_file

end module planet_config_io