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
  integer :: nx, ny
  integer :: n_mesh
  real(RP) :: delta_x, delta_y

  ! Temporal integration
  real(RP) :: tMax, dt
  integer :: ntime_step_max

end module Common_Variables_mod
