module Finite_Element_Resolution_mod
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  implicit none
  private

  public :: explicit_Temp_dt_compute

  contains

#ifdef heatFlux_CALC_2D

  subroutine explicit_Temp_dt_compute(Temp, Temp_dt) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    real(RP), dimension(:), intent(inout) :: Temp_dt

    real(RP) :: alpha

    integer :: kP, kN, kE, kS, kW

    Temp_dt = 0

    alpha = therm_conduct/(rho*cp*delta_x*delta_y)

    call explicit_Temp_dt_compute_boundaries_isolated(Temp, Temp_dt, alpha)

    do i = 2, nx-1
      do j = 2, ny-1

        kP = j   + ny*(i-1)
        kN = j   + ny*(i-1-1)
        kE = j+1 + ny*(i-1)
        kS = j   + ny*(i-1+1)
        kW = j-1 + ny*(i-1)


        Temp_dt(kP) = alpha*(delta_y/delta_x*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                           +delta_x/delta_y*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                           )

      enddo
    enddo

  end subroutine explicit_Temp_dt_compute !}}}

#elif defined heatFlux_CALC_3D

  subroutine explicit_Temp_dt_compute(Temp, Temp_dt) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    real(RP), dimension(:), intent(inout) :: Temp_dt

    real(RP) :: alpha

    integer :: kP, kN, kE, kS, kW
    integer :: kT, kB

    Temp_dt = 0

    ! alpha = 1/alpha with alpha the thermal diffusivity
    alpha = therm_conduct/(rho*cp)

    call explicit_Temp_dt_compute_boundaries_isolated(Temp, Temp_dt, alpha)

     do k = 2, nz-1
       do i = 2, nx-1
         do j = 2, ny-1

           kP = j   + ny*(i-1)   + nx*ny*(k-1)
           kN = j+1 + ny*(i-1)   + nx*ny*(k-1)
           kS = j-1 + ny*(i-1)   + nx*ny*(k-1)
           kE = j   + ny*(i-1+1) + nx*ny*(k-1)
           kW = j   + ny*(i-1-1) + nx*ny*(k-1)
           kT = j   + ny*(i-1)   + nx*ny*(k-1+1)
           kB = j   + ny*(i-1)   + nx*ny*(k-1-1)

           Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                               +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                               +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                               )

         enddo
       enddo
     enddo

  end subroutine explicit_Temp_dt_compute !}}}

#endif


! WITH BOUNDARIES FLUX EQUAL ZERO

#ifdef heatFlux_CALC_2D

  subroutine explicit_Temp_dt_compute_boundaries_isolated(Temp, Temp_dt, alpha) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    real(RP), dimension(:), intent(inout) :: Temp_dt
    real(RP), intent(in) :: alpha

    integer :: kP, kN, kE, kS, kW

    ! /// corners ///

    ! i=1 , j=1
    kP = 1   + ny*(1-1)
    kE = 1+1 + ny*(1-1)
    kS = 1   + ny*(1-1+1)

    Temp_dt(kP) = alpha*(delta_y/delta_x*(Temp(kE) - Temp(kP)) &
                       +delta_x/delta_y*(-1)*(Temp(kP) - Temp(kS)) &
                       )

    ! i=1 , j=ny
    kP = ny   + ny*(1-1)
    kS = ny   + ny*(1-1+1)
    kW = ny-1 + ny*(1-1)

    Temp_dt(kP) = alpha*(delta_y/delta_x*(-1)*(Temp(kP) - Temp(kW)) &
                       +delta_x/delta_y*(-1)*(Temp(kP) - Temp(kS)) &
                       )

    ! i=nx , j=1
    kP = 1   + ny*(nx-1)
    kN = 1   + ny*(nx-1-1)
    kE = 1+1 + ny*(nx-1)

    Temp_dt(kP) = alpha*(delta_y/delta_x*(Temp(kE) - Temp(kP)) &
                       +delta_x/delta_y*(Temp(kN) - Temp(kP)) &
                       )

    ! i=nx , j=ny
    kP = ny   + ny*(nx-1)
    kN = ny   + ny*(nx-1-1)
    kW = ny-1 + ny*(nx-1)

    Temp_dt(kP) = alpha*(delta_y/delta_x*(-1)*(Temp(kP) - Temp(kW)) &
                       +delta_x/delta_y*(Temp(kN) - Temp(kP)) &
                       )


    ! /// Between the corners ///

    ! i=2..nx-1 , j=1
    do i=2, nx-1
      kP = 1   + ny*(i-1)
      kN = 1   + ny*(i-1-1)
      kE = 1+1 + ny*(i-1)
      kS = 1   + ny*(i-1+1)

      Temp_dt(kP) = alpha*(delta_y/delta_x*(Temp(kE) - Temp(kP)) &
                         +delta_x/delta_y*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         )
    enddo

    ! i=2..nx-1 , j=ny
    do i=2, nx-1
      kP = ny   + ny*(i-1)
      kN = ny   + ny*(i-1-1)
      kS = ny   + ny*(i-1+1)
      kW = ny-1 + ny*(i-1)

      Temp_dt(kP) = alpha*(delta_y/delta_x*(-1)*(Temp(kP) - Temp(kW)) &
                         +delta_x/delta_y*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         )
    enddo

    ! i=1 , j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(1-1)
      kE = j+1 + ny*(1-1)
      kS = j   + ny*(1-1+1)
      kW = j-1 + ny*(1-1)

      Temp_dt(kP) = alpha*(delta_y/delta_x*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +delta_x/delta_y*(-1)*(Temp(kP) - Temp(kS)) &
                         )
    enddo

    ! i=nx , j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(i-1)
      kN = j   + ny*(i-1-1)
      kE = j+1 + ny*(i-1)
      kW = j-1 + ny*(i-1)

      Temp_dt(kP) = alpha*(delta_y/delta_x*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +delta_x/delta_y*(Temp(kN) - Temp(kP)) &
                         )
    enddo

  end subroutine explicit_Temp_dt_compute_boundaries_isolated !}}}

#elif defined heatFlux_CALC_3D

  subroutine explicit_Temp_dt_compute_boundaries_isolated(Temp, Temp_dt, alpha) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    real(RP), dimension(:), intent(inout) :: Temp_dt
    real(RP), intent(in) :: alpha

    integer :: kP, kN, kE, kS, kW, kB, kT

    ! /// corners /// {{{

    ! k=1, i=1, j=1
    kP = 1   + ny*(1-1)   + nx*ny*(1-1)
    kN = 1+1 + ny*(1-1)   + nx*ny*(1-1)
    kE = 1   + ny*(1-1+1) + nx*ny*(1-1)
    kT = 1   + ny*(1-1)   + nx*ny*(1-1+1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - Temp(kP)) &
                       +1/(delta_y**2)*(Temp(kN) - Temp(kP)) &
                       +1/(delta_z**2)*(Temp(kT) - Temp(kP)) &
                       )

    ! k=nz, i=1, j=1
    kP = 1   + ny*(1-1)   + nx*ny*(nz-1)
    kN = 1+1 + ny*(1-1)   + nx*ny*(nz-1)
    kE = 1   + ny*(1-1+1) + nx*ny*(nz-1)
    kB = 1   + ny*(1-1)   + nx*ny*(nz-1-1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - Temp(kP)           ) &
                       +1/(delta_y**2)*(Temp(kN) - Temp(kP)           ) &
                       +1/(delta_z**2)*(         - Temp(kP) + Temp(kB)) &
                       )

    ! k=nz, i=1, j=ny
    kP = ny   + ny*(1-1)   + nx*ny*(nz-1)
    kS = ny-1 + ny*(1-1)   + nx*ny*(nz-1)
    kE = ny   + ny*(1-1+1) + nx*ny*(nz-1)
    kB = ny   + ny*(1-1)   + nx*ny*(nz-1-1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - Temp(kP)           ) &
                       +1/(delta_y**2)*(         - Temp(kP) + Temp(kS)) &
                       +1/(delta_z**2)*(         - Temp(kP) + Temp(kB)) &
                       )

    ! k=1, i=1, j=ny
    kP = ny   + ny*(1-1)   + nx*ny*(1-1)
    kS = ny-1 + ny*(1-1)   + nx*ny*(1-1)
    kE = ny   + ny*(1-1+1) + nx*ny*(1-1)
    kT = ny   + ny*(1-1)   + nx*ny*(1-1+1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - Temp(kP)           ) &
                       +1/(delta_y**2)*(         - Temp(kP) + Temp(kS)) &
                       +1/(delta_z**2)*(Temp(kT) - Temp(kP)           ) &
                       )

    ! k=1, i=nx, j=1
    kP = 1   + ny*(nx-1)   + nx*ny*(1-1)
    kN = 1+1 + ny*(nx-1)   + nx*ny*(1-1)
    kW = 1   + ny*(nx-1-1) + nx*ny*(1-1)
    kT = 1   + ny*(nx-1)   + nx*ny*(1-1+1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(         - Temp(kP) + Temp(kW)) &
                       +1/(delta_y**2)*(Temp(kN) - Temp(kP)           ) &
                       +1/(delta_z**2)*(Temp(kT) - Temp(kP)           ) &
                       )

    ! k=nz, i=nx, j=1
    kP = 1   + ny*(nx-1)   + nx*ny*(nz-1)
    kN = 1+1 + ny*(nx-1)   + nx*ny*(nz-1)
    kW = 1   + ny*(nx-1-1) + nx*ny*(nz-1)
    kB = 1   + ny*(nx-1)   + nx*ny*(nz-1-1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(         - Temp(kP) + Temp(kW)) &
                       +1/(delta_y**2)*(Temp(kN) - Temp(kP)           ) &
                       +1/(delta_z**2)*(         - Temp(kP) + Temp(kB)) &
                       )

    ! k=1, i=nx, j=ny
    kP = ny   + ny*(nx-1)   + nx*ny*(1-1)
    kS = ny-1 + ny*(nx-1)   + nx*ny*(1-1)
    kW = ny   + ny*(nx-1-1) + nx*ny*(1-1)
    kT = ny   + ny*(nx-1)   + nx*ny*(1-1+1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(         - Temp(kP) + Temp(kW)) &
                       +1/(delta_y**2)*(         - Temp(kP) + Temp(kS)) &
                       +1/(delta_z**2)*(Temp(kT) - Temp(kP)           ) &
                       )

    ! k=nz, i=nx, j=ny
    kP = ny   + ny*(nx-1)   + nx*ny*(nz-1)
    kS = ny-1 + ny*(nx-1)   + nx*ny*(nz-1)
    kW = ny   + ny*(nx-1-1) + nx*ny*(nz-1)
    kB = ny   + ny*(nx-1)   + nx*ny*(nz-1-1)

    Temp_dt(kP) = alpha*(1/(delta_x**2)*(         - Temp(kP) + Temp(kW)) &
                       +1/(delta_y**2)*(         - Temp(kP) + Temp(kS)) &
                       +1/(delta_z**2)*(         - Temp(kP) + Temp(kB)) &
                       )

   !}}}

    ! /// edges /// !{{{

    ! k=2..nz-1, i=1, j=1
    do k=2, nz-1
      kP = 1   + ny*(1-1)   + nx*ny*(k-1)
      kN = 1+1 + ny*(1-1)   + nx*ny*(k-1)
      kE = 1   + ny*(1-1+1) + nx*ny*(k-1)
      kT = 1   + ny*(1-1)   + nx*ny*(k-1+1)
      kB = 1   + ny*(1-1)   + nx*ny*(k-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - Temp(kP)             ) &
                         +1/(delta_y**2)*(Temp(kN) - Temp(kP)             ) &
                         +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                         )
    enddo

    ! k=2..nz-1, i=1, j=ny
    do k=2, nz-1
      kP = ny   + ny*(1-1)   + nx*ny*(k-1)
      kS = ny-1 + ny*(1-1)   + nx*ny*(k-1)
      kE = ny   + ny*(1-1+1) + nx*ny*(k-1)
      kT = ny   + ny*(1-1)   + nx*ny*(k-1+1)
      kB = ny   + ny*(1-1)   + nx*ny*(k-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) -   Temp(kP)           ) &
                         +1/(delta_y**2)*(         -   Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                         )
    enddo

    ! k=2..nz-1, i=nx, j=ny
    do k=2, nz-1
      kP = ny   + ny*(nx-1)   + nx*ny*(k-1)
      kS = ny-1 + ny*(nx-1)   + nx*ny*(k-1)
      kW = ny   + ny*(nx-1-1) + nx*ny*(k-1)
      kT = ny   + ny*(nx-1)   + nx*ny*(k-1+1)
      kB = ny   + ny*(nx-1)   + nx*ny*(k-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(         -   Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(         -   Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                         )
    enddo

    ! k=2..nz-1, i=nx, j=1
    do k=2, nz-1
      kP = 1   + ny*(nx-1)   + nx*ny*(k-1)
      kN = 1+1 + ny*(nx-1)   + nx*ny*(k-1)
      kW = 1   + ny*(nx-1-1) + nx*ny*(k-1)
      kT = 1   + ny*(nx-1)   + nx*ny*(k-1+1)
      kB = 1   + ny*(nx-1)   + nx*ny*(k-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(         -   Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(Temp(kN) -   Temp(kP)           ) &
                         +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                         )
    enddo

    ! k=nz, i=1, j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(1-1)   + nx*ny*(nz-1)
      kN = j+1 + ny*(1-1)   + nx*ny*(nz-1)
      kS = j-1 + ny*(1-1)   + nx*ny*(nz-1)
      kE = j   + ny*(1-1+1) + nx*ny*(nz-1)
      kB = j   + ny*(1-1)   + nx*ny*(nz-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) -   Temp(kP)           ) &
                         +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(         -   Temp(kP) + Temp(kB)) &
                         )
    enddo


    ! k=1, i=1, j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(1-1)   + nx*ny*(1-1)
      kN = j+1 + ny*(1-1)   + nx*ny*(1-1)
      kS = j-1 + ny*(1-1)   + nx*ny*(1-1)
      kE = j   + ny*(1-1+1) + nx*ny*(1-1)
      kT = j   + ny*(1-1)   + nx*ny*(1-1+1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) -   Temp(kP)           ) &
                         +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(Temp(kT) -   Temp(kP)           ) &
                         )
    enddo

    ! k=1, i=nx, j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(nx-1)   + nx*ny*(1-1)
      kN = j+1 + ny*(nx-1)   + nx*ny*(1-1)
      kS = j-1 + ny*(nx-1)   + nx*ny*(1-1)
      kW = j   + ny*(nx-1-1) + nx*ny*(1-1)
      kT = j   + ny*(nx-1)   + nx*ny*(1-1+1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(         -   Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(Temp(kT) -   Temp(kP)           ) &
                         )
    enddo

    ! k=nz, i=nx, j=2..ny-1
    do j=2, ny-1
      kP = j   + ny*(nx-1)   + nx*ny*(nz-1)
      kN = j+1 + ny*(nx-1)   + nx*ny*(nz-1)
      kS = j-1 + ny*(nx-1)   + nx*ny*(nz-1)
      kW = j   + ny*(nx-1-1) + nx*ny*(nz-1)
      kB = j   + ny*(nx-1)   + nx*ny*(nz-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(         -   Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(         -   Temp(kP) + Temp(kB)) &
                         )
    enddo

    ! k=1, i=2..nx-1, j=1
    do i=2, nx-1
      kP = 1   + ny*(i-1)   + nx*ny*(1-1)
      kN = 1+1 + ny*(i-1)   + nx*ny*(1-1)
      kE = 1   + ny*(i-1+1) + nx*ny*(1-1)
      kW = 1   + ny*(i-1-1) + nx*ny*(1-1)
      kT = 1   + ny*(i-1)   + nx*ny*(1-1+1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(Temp(kN) -   Temp(kP)           ) &
                         +1/(delta_z**2)*(Temp(kT) -   Temp(kP)           ) &
                         )
    enddo

    ! k=nz, i=2..nx-1, j=1
    do i=2, nx-1
      kP = 1   + ny*(i-1)   + nx*ny*(nz-1)
      kN = 1+1 + ny*(i-1)   + nx*ny*(nz-1)
      kE = 1   + ny*(i-1+1) + nx*ny*(nz-1)
      kW = 1   + ny*(i-1-1) + nx*ny*(nz-1)
      kB = 1   + ny*(i-1)   + nx*ny*(nz-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(Temp(kN) -   Temp(kP)           ) &
                         +1/(delta_z**2)*(         -   Temp(kP) + Temp(kB)) &
                         )
    enddo

    ! k=nz, i=2..nx-1, j=ny
    do i=2, nx-1
      kP = ny   + ny*(i-1)   + nx*ny*(nz-1)
      kS = ny-1 + ny*(i-1)   + nx*ny*(nz-1)
      kE = ny   + ny*(i-1+1) + nx*ny*(nz-1)
      kW = ny   + ny*(i-1-1) + nx*ny*(nz-1)
      kB = ny   + ny*(i-1)   + nx*ny*(nz-1-1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(         -   Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(         -   Temp(kP) + Temp(kB)) &
                         )
    enddo

    ! k=1, i=2..nx-1, j=ny
    do i=2, nx-1
      kP = ny   + ny*(i-1)   + nx*ny*(1-1)
      kS = ny-1 + ny*(i-1)   + nx*ny*(1-1)
      kE = ny   + ny*(i-1+1) + nx*ny*(1-1)
      kW = ny   + ny*(i-1-1) + nx*ny*(1-1)
      kT = ny   + ny*(i-1)   + nx*ny*(1-1+1)

      Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                         +1/(delta_y**2)*(         -   Temp(kP) + Temp(kS)) &
                         +1/(delta_z**2)*(Temp(kT) -   Temp(kP)           ) &
                         )
    enddo

    !}}}

    ! /// faces /// !{{{

    ! k=2..nz-1, i=1, j=2..ny-1
    do k=2, nz-1
      do j=2, ny-1
        kP = j   + ny*(1-1)   + nx*ny*(k-1)
        kN = j+1 + ny*(1-1)   + nx*ny*(k-1)
        kS = j-1 + ny*(1-1)   + nx*ny*(k-1)
        kE = j   + ny*(1-1+1) + nx*ny*(k-1)
        kT = j   + ny*(1-1)   + nx*ny*(k-1+1)
        kB = j   + ny*(1-1)   + nx*ny*(k-1-1)

        Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) -   Temp(kP)           ) &
                           +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                           +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                           )
      enddo
    enddo

    ! k=2..nz-1, i=nx, j=2..ny-1
    do k=2, nz-1
      do j=2, ny-1
        kP = j   + ny*(nx-1)   + nx*ny*(k-1)
        kN = j+1 + ny*(nx-1)   + nx*ny*(k-1)
        kS = j-1 + ny*(nx-1)   + nx*ny*(k-1)
        kW = j   + ny*(nx-1-1) + nx*ny*(k-1)
        kT = j   + ny*(nx-1)   + nx*ny*(k-1+1)
        kB = j   + ny*(nx-1)   + nx*ny*(k-1-1)

        Temp_dt(kP) = alpha*(1/(delta_x**2)*(         -   Temp(kP) + Temp(kW)) &
                           +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                           +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                           )
      enddo
    enddo

    ! k=2..nz-1, i=2..nx-1, j=1
    do k=2, nz-1
      do i=2, nx-1
        kP = 1   + ny*(i-1)   + nx*ny*(k-1)
        kN = 1+1 + ny*(i-1)   + nx*ny*(k-1)
        kE = 1   + ny*(i-1+1) + nx*ny*(k-1)
        kW = 1   + ny*(i-1-1) + nx*ny*(k-1)
        kT = 1   + ny*(i-1)   + nx*ny*(k-1+1)
        kB = 1   + ny*(i-1)   + nx*ny*(k-1-1)

        Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                           +1/(delta_y**2)*(Temp(kN) -   Temp(kP)           ) &
                           +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                           )
      enddo
    enddo

    ! k=2..nz-1, i=2..nx-1, j=ny
    do k=2, nz-1
      do i=2, nx-1
        kP = ny   + ny*(i-1)   + nx*ny*(k-1)
        kS = ny-1 + ny*(i-1)   + nx*ny*(k-1)
        kE = ny   + ny*(i-1+1) + nx*ny*(k-1)
        kW = ny   + ny*(i-1-1) + nx*ny*(k-1)
        kT = ny   + ny*(i-1)   + nx*ny*(k-1+1)
        kB = ny   + ny*(i-1)   + nx*ny*(k-1-1)

        Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                           +1/(delta_y**2)*(         -   Temp(kP) + Temp(kS)) &
                           +1/(delta_z**2)*(Temp(kT) - 2*Temp(kP) + Temp(kB)) &
                           )
      enddo
    enddo

    ! k=1, i=2..nx-1, j=2..ny-1
    do j=2, ny-1
      do i=2, nx-1
        kP = j   + ny*(i-1)   + nx*ny*(1-1)
        kN = j+1 + ny*(i-1)   + nx*ny*(1-1)
        kS = j-1 + ny*(i-1)   + nx*ny*(1-1)
        kE = j   + ny*(i-1+1) + nx*ny*(1-1)
        kW = j   + ny*(i-1-1) + nx*ny*(1-1)
        kT = j   + ny*(i-1)   + nx*ny*(1-1+1)

        Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                           +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                           +1/(delta_z**2)*(Temp(kT) -   Temp(kP)           ) &
                           )
      enddo
    enddo

    ! k=nz, i=2..nx-1, j=2..ny-1
    do j=2, ny-1
      do i=2, nx-1
        kP = j   + ny*(i-1)   + nx*ny*(nz-1)
        kN = j+1 + ny*(i-1)   + nx*ny*(nz-1)
        kS = j-1 + ny*(i-1)   + nx*ny*(nz-1)
        kE = j   + ny*(i-1+1) + nx*ny*(nz-1)
        kW = j   + ny*(i-1-1) + nx*ny*(nz-1)
        kB = j   + ny*(i-1)   + nx*ny*(nz-1-1)

        Temp_dt(kP) = alpha*(1/(delta_x**2)*(Temp(kE) - 2*Temp(kP) + Temp(kW)) &
                           +1/(delta_y**2)*(Temp(kN) - 2*Temp(kP) + Temp(kS)) &
                           +1/(delta_z**2)*(         -   Temp(kP) + Temp(kB)) &
                           )
      enddo
    enddo

    !}}}

  end subroutine explicit_Temp_dt_compute_boundaries_isolated !}}}

#endif

end module Finite_Element_Resolution_mod
