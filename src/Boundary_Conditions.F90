module Boundary_Conditions_mod
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  use IO_mod

  implicit none
  private

  public :: update_boundary_conditions

  contains

  subroutine update_boundary_conditions()

  select case (boundary_conditions_type)
    case ("Dirichlet") ! Regarding the temperature :
      call set_Temp_center(Temp, 393.d0)
    case ("Neumann") ! Regarding the heat flux :
      call read_heatGenFile_perTimeStep(iunit_heatGen, heat_surface_gen)
    case default
      stop
  end select

  end subroutine update_boundary_conditions

  subroutine set_Temp_center(Temp, T_center) !{{{
    real(RP), dimension(:), intent(inout) :: Temp
    real(RP), intent(in) :: T_center

#ifdef heatFlux_CALC_2D
    Temp(ny/2 + ny*(nx/2-1)) = T_center
#elif defined heatFlux_CALC_3D
    Temp(ny/2   + ny*(nx/2  -1) + ny*nx*(nz-1)) = T_center
    Temp(ny/2+1 + ny*(nx/2  -1) + ny*nx*(nz-1)) = T_center
    Temp(ny/2-1 + ny*(nx/2  -1) + ny*nx*(nz-1)) = T_center
    Temp(ny/2   + ny*(nx/2+1-1) + ny*nx*(nz-1)) = T_center
    Temp(ny/2   + ny*(nx/2-1-1) + ny*nx*(nz-1)) = T_center
    !do k=1,nz
    !  Temp(ny/2 + ny*(nx/2-1) + ny*nx*(k-1)) = T_center
    !enddo
#endif

  end subroutine set_Temp_center !}}}

end module Boundary_Conditions_mod
