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
  integer :: ncounter_write_file
  real(RP) :: t

  real(RP), dimension(:), allocatable :: Temp, Temp_dt

  namelist /material_prop/ therm_conduct, cp, rho
IN_2D  namelist /mesh_prop/ length, width, nx, ny
IN_3D  namelist /mesh_prop/ length, width, thickness, nx, ny, nz
  namelist /temporal_integration/ tMax, dt, ntime_step_max, noutput_file

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
IN_3D delta_z = thickness / real(nz)

IN_2D  n_mesh = nx*ny                     ! global variable
IN_3D  n_mesh = nx*ny*nz                  ! global variable
  allocate(Temp(n_mesh), Temp_dt(n_mesh))
! Tijk = T(j + ny*(i-1) + nz(k-1))

  call Initialize_Temperature(Temp)
  Temp_dt(:) = 0

  write(*,'(a70)') '----------------------------------------------------------------------'
  write(*,*) 'Starting main loop in time ...'

  ncounter_write_file = 0
  ntime_step = 0
  t = 0.d0
  do while (t<tMax .and. ntime_step < ntime_step_max)

    write(*,'(a70)') '----------------------------------------------------------------------'
    write(*,'(a12,I7,a7,E13.6,a2)')' ntime_step=',ntime_step,'     t=',t,' s'

    call write_file_Temp_area(Temp, ncounter_write_file, t)

    call update_Temp_dt(Temp, Temp_dt)

    call Explicit_Euler(Temp, Temp_dt, dt)

    t=t+dt
    ntime_step = ntime_step + 1

  enddo !while loop in time

  deallocate(Temp, Temp_dt)

end program heatFlux
