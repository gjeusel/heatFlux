module Common_Variables_mod
#include "config.h"

  use Parameters_mod

  implicit none
  !everything public

  integer :: i,j,k

  ! Material properties
  real(RP) :: therm_conduct ! thermal conductivity
  real(RP) :: cp ! specific heat
  real(RP) :: rho ! density

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
  integer :: noutput_file

end module Common_Variables_mod
