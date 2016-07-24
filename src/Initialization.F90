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

    Temp(ny/2 + ny*(nx/2-1)) = 373.d0

  end subroutine Initialize_Temperature


end module Initialization_mod
