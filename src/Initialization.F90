module Initialization_mod
#include "config.h"

  use Common_Variables_mod

  implicit none
  private

  public :: Initialize_Temperature

  contains

  subroutine Initialize_Temperature(Temp)
    real(RP), dimension(:), intent(inout) :: Temp

    do k=1, n_mesh
      Temp(k) = 293.d0
    enddo

    !do j=1, ny
    !  Temp(j) = 373.d0
    !enddo

#ifdef heatFlux_CALC_2D

    Temp(ny/2 + ny*(nx/2-1)) = 373.d0

#elif defined heatFlux_CALC_3D

    do k=1,nz
      Temp(ny/2 + ny*(nx/2-1) + ny*nx*(k-1)) = 373.d0
    enddo

#endif

  end subroutine Initialize_Temperature


end module Initialization_mod
