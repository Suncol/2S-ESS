program test_planetary_config
  use planet_config_mod
  use planet_config_io
  implicit none
  
  ! Define dp parameter
  integer, parameter :: dp = kind(1.0d0)
  
  ! Test Earth configuration (default)
  call test_earth_config()
  
  ! Test Mars configuration
  call test_mars_config()
  
  ! Test Venus configuration
  call test_venus_config()
  
  ! Test configuration file loading
  call test_config_file_loading()
  
  print *, 'All planetary configuration tests passed!'

contains

  subroutine test_earth_config()
    type(planet_config_type) :: config
    call set_earth_config()
    config = get_planet_config()
    
    if (config%radius_km /= 6371.0_dp) then
      stop 'Earth radius test failed!'
    endif
    
    if (.not. config%has_oceans) then
      stop 'Earth ocean test failed!'
    endif
    
    if (config%salinity_reference_ppt /= 34.3_dp) then
      stop 'Earth salinity test failed!'
    endif
    
    print *, 'Earth configuration test: PASSED'
  end subroutine test_earth_config
  
  subroutine test_mars_config()
    type(planet_config_type) :: config
    call set_mars_config()
    config = get_planet_config()
    
    if (config%radius_km /= 3389.5_dp) then
      stop 'Mars radius test failed!'
    endif
    
    if (config%has_oceans) then
      stop 'Mars ocean test failed!'
    endif
    
    if (config%name /= 'Mars') then
      stop 'Mars name test failed!'
    endif
    
    print *, 'Mars configuration test: PASSED'
  end subroutine test_mars_config
  
  subroutine test_venus_config()
    type(planet_config_type) :: config
    call set_venus_config()
    config = get_planet_config()
    
    if (config%radius_km /= 6051.8_dp) then
      stop 'Venus radius test failed!'
    endif
    
    if (config%has_oceans) then
      stop 'Venus ocean test failed!'
    endif
    
    if (config%name /= 'Venus') then
      stop 'Venus name test failed!'
    endif
    
    print *, 'Venus configuration test: PASSED'
  end subroutine test_venus_config
  
  subroutine test_config_file_loading()
    type(planet_config_type) :: config
    
    ! Test loading Mars from configuration file
    call load_planet_config_from_file('INPUT/mars_config.nml')
    config = get_planet_config()
    
    if (config%radius_km /= 3389.5_dp) then
      stop 'Mars config file test failed!'
    endif
    
    if (config%name /= 'Mars') then
      stop 'Mars config file name test failed!'
    endif
    
    ! Test loading Earth from configuration file  
    call load_planet_config_from_file('INPUT/earth_config.nml')
    config = get_planet_config()
    
    if (config%radius_km /= 6371.0_dp) then
      stop 'Earth config file test failed!'
    endif
    
    if (.not. config%has_oceans) then
      stop 'Earth config file ocean test failed!'
    endif
    
    print *, 'Configuration file loading test: PASSED'
  end subroutine test_config_file_loading

end program test_planetary_config