program heatFlux
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  use IO_mod
  use Initialization_mod

  use Boundary_Conditions_mod
  use Finite_Element_Resolution_mod
  use TVR_mod

  implicit none

  character(len=STRING_LENGTH) :: namelist_path
  integer :: ntime_step
  integer :: ncounter_write_file
  real(RP) :: t

  if (command_argument_count() <= 0) then
    print *,'You forgot the path to the namelist'
    stop
  else
    call getarg(1,namelist_path)
  endif

  call initialize_default_values()
  call read_nml(namelist_path)

  ! Initial conditions :
  Temp(:) = 293.d0

  write(*,'(a70)') '----------------------------------------------------------------------'
  write(*,*) 'Starting main loop in time ...'

  ncounter_write_file = 0
  ntime_step = 0
  t = 0.d0
  do while (t<tMax .and. ntime_step < ntime_step_max)

    write(*,'(a70)') '----------------------------------------------------------------------'
    write(*,'(a12,I7,a7,E13.6,a2)')' ntime_step=',ntime_step,'     t=',t,' s'

    call update_boundary_conditions()

    call write_physical_datas(Temp, ncounter_write_file, t)

    call explicit_Temp_dt_compute(Temp, Temp_dt)

    call Explicit_Euler(Temp, Temp_dt, dt)

    t=t+dt
    ntime_step = ntime_step + 1

  enddo !while loop in time

  deallocate(Temp, Temp_dt)
  deallocate(heat_surface_gen)
  close(iunit_heatGen)

end program heatFlux
