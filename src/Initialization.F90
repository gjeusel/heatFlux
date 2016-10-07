module Initialization_mod
#include "config.h"

  use Common_Variables_mod
  use IO_mod

  implicit none
  private

  public :: initialize_default_values
  public :: read_nml

  contains

  subroutine initialize_default_values() !{{{

  ! ---> Material properties of steel
  therm_conduct = 44 ! (W/m.K)
  cp = 1.81d6        ! (J/K.kg)
  rho = 1.d0         ! (kg/m³)

  ! ---> Mesh prop
  length = 0.01
  width = 0.01
IN_3D thickness = 0.001

  nx = 16
  ny = 16
  IN_3D nz = 5

IN_2D  n_mesh = nx*ny                     ! global variable
IN_3D  n_mesh = nx*ny*nz                  ! global variable

  delta_x = length / real(nx)
  delta_y = width / real(ny)
IN_3D delta_z = thickness / real(nz)

  ! ---> Temporal integration
  tMax = 1.895
  dt = 0.005
  ntime_step_max = 10

  ! ----> Boundary conditions
  boundary_conditions_type = "Dirichlet"

  ! ---> Output
  noutput_file = 10
  o_base_name = "temperature_field"
  o_format = "csv"
  o_field_area = .false.

  end subroutine initialize_default_values !}}}

  subroutine read_nml(namelist_path) !{{{
    character(len=STRING_LENGTH), intent(in) :: namelist_path
    integer :: iunit_nml
    logical :: file_exists
    logical :: nml_exists

    namelist /material_prop/ therm_conduct, cp, rho
!IN_2D  namelist /mesh_prop/ length, width, nx, ny
IN_2D  namelist /mesh_prop/ nx, ny, length, width
IN_3D  namelist /mesh_prop/ length, width, thickness, nx, ny, nz
    namelist /temporal_integration/ tMax, dt, ntime_step_max, noutput_file
    namelist /boundary_conditions/ boundary_conditions_type, heat_Generation_file
    namelist /output_files/ noutput_file, o_base_name, o_format, o_field_area, o_field_area_z

    call opening_file(namelist_path, iunit_nml)

    ! ---> Material properties
    read(iunit_nml, nml=material_prop)

    ! ---> Mesh prop
    read(iunit_nml, nml=mesh_prop)
    delta_x = length / real(nx)
    delta_y = width / real(ny)
IN_3D delta_z = thickness / real(nz)

IN_2D  n_mesh = nx*ny                     ! global variable
IN_3D  n_mesh = nx*ny*nz                  ! global variable

    ! Allocation of global arrays
    allocate(Temp(n_mesh), Temp_dt(n_mesh))
    allocate(heat_surface_gen(nx*ny))

    ! ---> Temporal integration
    read(iunit_nml, nml=temporal_integration)

    ! ---> Boundary conditions
    read(iunit_nml, nml=boundary_conditions)

    if (boundary_conditions_type=="Dirichlet") then
      heat_surface_gen(:) = 0
    endif

    ! if Dirichlet condition choosed but path to heat_Generation_file given,
    ! set Warning and choose Dirichlet conditions by default prioritary.
    if(boundary_conditions_type=="Dirichlet" .and. len_trim(heat_Generation_file) /= STRING_LENGTH) then
      write(*,*) "Warning nml : boundary_conditions_type set to Dirichlet, and path for heat_Generation_file given."
      write(*,*) "Boundary conditions set to Dirichlet by default and continue..."
    endif

    ! If heat_Generation_file is altered by the namelist :
    if (len_trim(heat_Generation_file) /= STRING_LENGTH) then
      inquire(file=heat_Generation_file, exist=file_exists)
      if (file_exists) then
        open(newunit=iunit_heatGen, file=heat_Generation_file, action = 'read')
      else ! -> without heat source :
        write(*,*) "Error nml path : ", trim(adjustl(heat_Generation_file)), " does not exist."
        stop
      endif
    endif

    ! ---> Output
    read(iunit_nml, nml=output_files)
#ifdef heatFlux_CALC_3D
    if (o_field_area_z > thickness) then
      write(*,"(A33,F12.8,A23)") "Error nml : The plan o_field_area_z =",o_field_area_z," is out of the geometry"
      stop
    endif
#endif

  end subroutine read_nml !}}}

end module Initialization_mod
