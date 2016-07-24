program heatFlux
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  use IO_mod
  use Initialization_mod

  use update_Temp_mod
  use TVR_mod

  implicit none

  character(STRING_LENGTH) :: namelist_path
  integer :: iunit
  integer :: ntime_step
  real(RP) :: t

  real(RP), dimension(:), allocatable :: Temp, Temp_dt

  namelist /material_prop/ therm_conduct, cp, rho
  namelist /mesh_prop/ length, width, nx, ny
  namelist /temporal_integration/ tMax, dt, ntime_step_max

  if (command_argument_count() <= 0) then
    print *,'You forgot the path to the namelist'
    stop
  else
    call getarg(1,namelist_path)
    call opening_file(namelist_path, iunit)
  endif

  read(iunit, nml=material_prop)
  read(iunit, nml=mesh_prop)
  read(iunit, nml=temporal_integration)

  delta_x = length / real(nx)
  delta_y = width / real(ny)

  n_mesh = nx*ny !global variable
  allocate(Temp(n_mesh), Temp_dt(n_mesh))

  call Initialize_Temperature(Temp)
  Temp_dt(:) = 0

  write(*,'(a70)') '----------------------------------------------------------------------'
  write(*,*) 'Starting main loop in time ...'

  ntime_step = 0
  t = 0.d0
  do while (t<tMax .and. ntime_step < ntime_step_max)

    write(*,'(a70)') '----------------------------------------------------------------------'
    write(*,'(a12,I7,a7,E13.6,a2)')' ntime_step=',ntime_step,'     t=',t,' s'

    call update_Temp_dt(Temp, Temp_dt)

    call Explicit_Euler(Temp, Temp_dt, dt)

    t=t+dt
    ntime_step = ntime_step + 1

  enddo !while loop in time

  open (unit=7, file='Temperature.dat', status='replace')
  do i=1, nx
    do j=1, ny
      write(7,*) i*delta_x -0.5*delta_x, j*delta_y -0.5*delta_y, Temp(j + ny*(i-1)), Temp_dt(j + ny*(i-1))
    enddo
      write(7,*)
  enddo
  close(7)

  deallocate(Temp, Temp_dt)

end program heatFlux
