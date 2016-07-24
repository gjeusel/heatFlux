module TVR_mod

  use Parameters_mod

  implicit none
  private

  public :: Explicit_Euler

  contains

! Explicit Euler time integration solver
  subroutine Explicit_Euler(Temp, Temp_dt, dt)
    implicit none

    real(RP), dimension(:), intent(inout) :: Temp
    real(RP), dimension(:), intent(inout) :: Temp_dt
    real(RP), intent(in) :: dt

    Temp = Temp + dt*Temp_dt

  end subroutine Explicit_Euler

end module TVR_mod
