module update_Temp_mod
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  implicit none
  private

  public :: update_Temp_dt

  contains

  subroutine update_Temp_dt(Temp, Temp_dt)
    real(RP), dimension(:), intent(in) :: Temp
    real(RP), dimension(:), intent(inout) :: Temp_dt

    real(RP) :: beta

    integer :: kP, kN, kE, kS, kW

    Temp_dt = 0

    beta = therm_conduct/(rho*cp*delta_x*delta_y)

    call update_Temp_dt_boundaries_isolated(Temp, Temp_dt, beta)

    do i = 2, nx-1
      do j = 2, ny-1

        kP = j   + ny*(i-1)
        kN = j   + ny*(i-1-1)
        kE = j+1 + ny*(i-1)
        kS = j   + ny*(i-1+1)
        kW = j-1 + ny*(i-1)

        Temp_dt(kP) = beta*(delta_y/delta_x*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                           +delta_x/delta_y*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                           )

      enddo
    enddo

  end subroutine update_Temp_dt


! WITH BOUNDARIES FLUX EQUAL ZERO
  subroutine update_Temp_dt_boundaries_isolated(Temp, Temp_dt, beta)
    real(RP), dimension(:), intent(in) :: Temp
    real(RP), dimension(:), intent(inout) :: Temp_dt
    real(RP), intent(in) :: beta

    integer :: kP, kN, kE, kS, kW

    ! /// corners ///

    ! i=1 , j=1
    kP = 1   + ny*(1-1)
    kE = 1+1 + ny*(1-1)
    kS = 1   + ny*(1-1+1)

    Temp_dt(kP) = beta*(delta_y/delta_x*(Temp(kE) - Temp(kP)) &
                       +delta_x/delta_y*(-1)*(Temp(kP) - Temp(kS)) &
                       )

    ! i=1 , j=ny
    kP = ny   + ny*(1-1)
    kS = ny   + ny*(1-1+1)
    kW = ny-1 + ny*(1-1)

    Temp_dt(kP) = beta*(delta_y/delta_x*(-1)*(Temp(kP) - Temp(kW)) &
                       +delta_x/delta_y*(-1)*(Temp(kP) - Temp(kS)) &
                       )

    ! i=nx , j=1
    kP = 1   + ny*(nx-1)
    kN = 1   + ny*(nx-1-1)
    kE = 1+1 + ny*(nx-1)

    Temp_dt(kP) = beta*(delta_y/delta_x*(Temp(kE) - Temp(kP)) &
                       +delta_x/delta_y*(Temp(kN) - Temp(kP)) &
                       )

    ! i=nx , j=ny
    kP = ny   + ny*(nx-1)
    kN = ny   + ny*(nx-1-1)
    kW = ny-1 + ny*(nx-1)

    Temp_dt(kP) = beta*(delta_y/delta_x*(-1)*(Temp(kP) - Temp(kW)) &
                       +delta_x/delta_y*(Temp(kN) - Temp(kP)) &
                       )


    ! /// Between the corners ///

    ! i=2..nx-1 , j=1
    do i=2, nx-1
      kP = 1   + ny*(i-1)
      kN = 1   + ny*(i-1-1)
      kE = 1+1 + ny*(i-1)
      kS = 1   + ny*(i-1+1)

      Temp_dt(kP) = beta*(delta_y/delta_x*(Temp(kE) - Temp(kP)) &
                         +delta_x/delta_y*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         )
    enddo

    ! i=2..nx-1 , j=ny
    do i=2, nx-1
      kP = ny   + ny*(i-1)
      kN = ny   + ny*(i-1-1)
      kS = ny   + ny*(i-1+1)
      kW = ny-1 + ny*(i-1)

      Temp_dt(kP) = beta*(delta_y/delta_x*(-1)*(Temp(kP) - Temp(kW)) &
                         +delta_x/delta_y*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         )
    enddo

    ! i=1 , j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(1-1)
      kE = j+1 + ny*(1-1)
      kS = j   + ny*(1-1+1)
      kW = j-1 + ny*(1-1)

      Temp_dt(kP) = beta*(delta_y/delta_x*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +delta_x/delta_y*(-1)*(Temp(kP) - Temp(kS)) &
                         )
    enddo

    ! i=nx , j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(i-1)
      kN = j   + ny*(i-1-1)
      kE = j+1 + ny*(i-1)
      kW = j-1 + ny*(i-1)

      Temp_dt(kP) = beta*(delta_y/delta_x*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +delta_x/delta_y*(Temp(kN) - Temp(kP)) &
                         )
    enddo

  end subroutine update_Temp_dt_boundaries_isolated

end module update_Temp_mod
