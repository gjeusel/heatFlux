module Common_Variables_mod
#include "config.h"

  use Parameters_mod

  implicit none
  !everything public

  integer :: i,j,k

  ! Physical variables :
  real(RP), dimension(:), allocatable :: Temp ! Tijk = Temp(j + ny*(i-1) + nx*ny*(k-1))
  real(RP), dimension(:), allocatable :: Temp_dt ! T_dt ijk = Temp_dt(j + ny*(i-1) + nx*ny*(k-1))

  ! Material properties
  real(RP) :: therm_conduct ! thermal conductivity
  real(RP) :: cp            ! specific heat
  real(RP) :: rho           ! density

  ! Spatial decomposition
  real(RP) :: length, width
IN_3D real(RP) :: thickness
  integer :: nx, ny
IN_3D integer :: nz
  integer :: n_mesh
  real(RP) :: delta_x, delta_y
IN_3D real(RP) :: delta_z

  ! Temporal integration
  real(RP) :: tMax, dt
  integer :: ntime_step_max

  ! Output
  integer :: noutput_file
  character(len=STRING_LENGTH) :: o_base_name
  character(len=STRING_LENGTH) :: o_format
  logical :: o_field_area
  real(RP) :: o_field_area_z

  ! Boundary conditions
  integer :: iunit_heatGen
  character(len=STRING_LENGTH) :: heat_Generation_file
  real(RP), dimension(:), allocatable :: heat_surface_gen ! heat_surface_gen ij = T(j + ny*(i-1))
  character(len=STRING_LENGTH) :: boundary_conditions_type

  ! Initial conditions

end module Common_Variables_mod
