module IO_mod
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  implicit none
  private

  public :: opening_file
  public :: write_file_Temp
  public :: write_file_Temp_area

  contains

  subroutine opening_file(filename, iunit) !{{{
    character(len=*), intent(in) :: filename
    integer, intent(inout) :: iunit

    integer :: error
    logical :: file_exists, file_opened

    inquire(file=filename, opened=file_opened)

    if(file_opened) then
      return
    else
      inquire(file=filename, exist=file_exists)
      if(file_exists) then
        open(newunit=iunit, file=filename, iostat=error)
        if(error /= 0) then
          print *, 'Opening ', trim(filename), ' failed ...'
          print *, 'iostat = ', error
          stop
        endif
      else
        print *, 'File named : \" ', trim(filename), '\" not found.'
        stop
      endif
    endif

  end subroutine opening_file !}}}

  subroutine write_file_Temp_area(Temp, ncounter_write_file, t) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    integer, intent(inout) :: ncounter_write_file
    real(RP), intent(in) :: t

    integer :: iunit, error
    character(len=STRING_LENGTH) :: filename, num_file
    character(len=STRING_LENGTH) :: real_format, row_format
    real(RP) :: pos_x, pos_y, pos_z

    if (ncounter_write_file*tMax/noutput_file <= t) then

      write(num_file, *) ncounter_write_file
      write(filename, *) "Temperature_area_", trim(adjustl(num_file)), ".csv"

      open(newunit=iunit, file=trim(adjustl(filename)), iostat=error)

      call write_banner_file(iunit, t)
      !write(iunit, *), "x coord, y coord, Temperature"

      real_format = "F15.10"
      row_format = "(" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // ")"

IN_3D k = nz

      do j=1, ny
        pos_y = (j-0.5)*delta_y

        do i=1, nx
          pos_x = (i-0.5)*delta_x

IN_2D     write(iunit,row_format) pos_x, ", ", pos_y, ", ", Temp(j + ny*(i-1))
IN_3D     write(iunit,row_format) pos_x, ", ", pos_y, ", ", Temp(j + ny*(i-1) + nx*ny*(k-1))

        enddo
        write(iunit,*)
      enddo

      close(iunit)
      ncounter_write_file = ncounter_write_file + 1
    endif

  end subroutine write_file_Temp_area !}}}

  subroutine write_file_Temp(Temp, ncounter_write_file, t) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    integer, intent(inout) :: ncounter_write_file
    real(RP), intent(in) :: t

    integer :: iunit, error
    character(len=STRING_LENGTH) :: filename, num_file
    character(len=STRING_LENGTH) :: real_format, row_format
    real(RP) :: pos_x, pos_y, pos_z

    if (ncounter_write_file*tMax/noutput_file <= t) then

      write(num_file, *) ncounter_write_file
      write(filename, *) "Temperature_", trim(adjustl(num_file)), ".csv"

      open(newunit=iunit, file=trim(adjustl(filename)), iostat=error)

      call write_banner_file(iunit, t)
      !write(iunit, *), "x coord, y coord, z coord, Temperature"

      real_format = "F15.10"
IN_2D row_format = "(" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // ")"
IN_3D row_format = "(" // trim(adjustl(real_format)) // "A2" //trim(adjustl(real_format)) &
IN_3D              // "A2" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // ")"

IN_3D do k=1, nz
IN_3D   pos_z = (k-0.5)*delta_z

        do j=1, ny
          pos_y = (j-0.5)*delta_y

          do i=1, nx
            pos_x = (i-0.5)*delta_x

IN_2D       write(iunit,row_format) pos_x, ", ", pos_y, ", ", Temp(j + ny*(i-1))
IN_3D       write(iunit,row_format) pos_x, ", ", pos_y, ", ", pos_z, ", ", Temp(j + ny*(i-1) + nx*ny*(k-1))

          enddo
        enddo
IN_3D enddo

      close(iunit)
      ncounter_write_file = ncounter_write_file + 1
    endif

  end subroutine write_file_Temp !}}}

  subroutine write_banner_file(iunit, t) !{{{
    integer, intent(in) :: iunit
    real(RP), intent(in) :: t
    character(len=STRING_LENGTH) :: tmp_str, tmp_str_1, tmp_str_2, tmp_str_3


    write(iunit, *) "##############################################"

    write(tmp_str, "(F10.4)") t
    write(iunit, *) "# t = ", trim(adjustl(tmp_str)), " s"

    write(tmp_str, "(F10.3)") therm_conduct
    write(iunit, *) "# k = ", trim(adjustl(tmp_str)), " W/m².K"

    write(tmp_str, "(F10.3)") cp
    write(iunit, *) "# cp = ", trim(adjustl(tmp_str)), " J/K.kg"

    write(tmp_str, "(F10.3)") rho
    write(iunit, *) "# rho = ", trim(adjustl(tmp_str)), " kg/m³"

    write(tmp_str_1, *) nx
    write(tmp_str_2, *) ny
IN_3D  write(tmp_str_3, *) nz
#ifdef heatFlux_CALC_2D
    write(tmp_str, *) "# nx = ", trim(adjustl(tmp_str_1)) &
                     ,", ny = ", trim(adjustl(tmp_str_2))
#elif defined heatFlux_CALC_3D
    write(tmp_str, *) "# nx = ", trim(adjustl(tmp_str_1)) &
                     ,", ny = ", trim(adjustl(tmp_str_2)) &
                     ,", nz = ", trim(adjustl(tmp_str_3))
#endif
    write(iunit, *) trim(adjustl(tmp_str))

    write(tmp_str_1, "(F10.3)") length
    write(tmp_str_2, "(F10.3)") width
IN_3D    write(tmp_str_3, "(F10.3)") thickness
#ifdef heatFlux_CALC_2D
    write(tmp_str, *) "# length = ", trim(adjustl(tmp_str_1)), "m" &
                     ,", width = ", trim(adjustl(tmp_str_2)), "m"
#elif defined heatFlux_CALC_3D
    write(tmp_str, *) "# length = ", trim(adjustl(tmp_str_1)), "m" &
                     ,", width = ", trim(adjustl(tmp_str_2)), "m" &
                     ,", thickness = ", trim(adjustl(tmp_str_3)), "m"
#endif
    write(iunit, *) trim(adjustl(tmp_str))

    write(iunit, *) "##############################################"
    write(iunit, *)

  end subroutine write_banner_file !}}}

end module IO_mod
