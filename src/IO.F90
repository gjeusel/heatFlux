module IO_mod
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  implicit none
  private

  public :: opening_file
  public :: write_physical_datas
  public :: read_heatGenFile_perTimeStep

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

  subroutine write_physical_datas(Temp, ncounter_write_file, t) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    integer, intent(inout) :: ncounter_write_file
    real(RP), intent(in) :: t

    character(len=STRING_LENGTH) :: o_filename, num_file
    integer :: iunit, error

    if (ncounter_write_file*tMax/noutput_file <= t) then

      write(num_file, *) ncounter_write_file
      write(o_filename, *) trim(adjustl(o_base_name)),"_",trim(adjustl(num_file)),".", trim(adjustl(o_format))
      open(newunit=iunit, file=trim(adjustl(o_filename)), iostat=error)

      call write_banner_file(iunit, t)
      call write_Full_Temp_field(Temp, iunit)
      close(iunit)

IN_3D if (o_field_area) then
IN_3D   write(o_filename, *) "plan_",trim(adjustl(o_base_name)),"_",trim(adjustl(num_file)),".", trim(adjustl(o_format))
IN_3D   open(newunit=iunit, file=trim(adjustl(o_filename)), iostat=error)
IN_3D   call write_banner_file(iunit, t)
IN_3D   call write_Plan_Temp_field(Temp, iunit)
IN_3D   close(iunit)
IN_3D endif

      ncounter_write_file = ncounter_write_file + 1
    endif

  end subroutine write_physical_datas !}}}

  subroutine write_Full_Temp_field(Temp, iunit) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    integer, intent(in) :: iunit

    character(len=STRING_LENGTH) :: real_format, row_format
    real(RP) :: pos_x, pos_y, pos_z

      ! Regarding the filetype output format desired :
      if (trim(adjustl(o_format)) == "csv") then
        real_format = "F15.10"
IN_2D   row_format = "(" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // ")"
IN_3D   row_format = "(" // trim(adjustl(real_format)) // "A2" //trim(adjustl(real_format)) &
IN_3D                // "A2" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // ")"

IN_2D   write(iunit, *) "X,               Y,               T,"
IN_3D   write(iunit, *) "X,               Y,               Z,               T,"

IN_3D   do k=1, nz
IN_3D     pos_z = (k-0.5)*delta_z

          do j=1, ny
            pos_y = (j-0.5)*delta_y

            do i=1, nx
              pos_x = (i-0.5)*delta_x

IN_2D         write(iunit,row_format) pos_x, ", ", pos_y, ", ", Temp(j + ny*(i-1))
IN_3D         write(iunit,row_format) pos_x, ", ", pos_y, ", ", pos_z, ", ", Temp(j + ny*(i-1) + nx*ny*(k-1))

            enddo
          write(iunit,*)
          enddo
IN_3D   write(iunit,*)
IN_3D   enddo
      endif

  end subroutine write_Full_Temp_field !}}}

#ifdef heatFlux_CALC_3D
  subroutine write_Plan_Temp_field(Temp, iunit) !{{{
    real(RP), dimension(:), intent(in) :: Temp
    integer, intent(in) :: iunit

    character(len=STRING_LENGTH) :: real_format, row_format
    real(RP) :: pos_x, pos_y, pos_z

    integer :: k_z

    pos_z = o_field_area_z ! from namelist

    if (pos_z == 0) then
      k_z = 1
    else
      k_z = nint(pos_z/delta_z) ! integer closer to the real
    endif

    ! Regarding the filetype output format desired :
    if (trim(adjustl(o_format)) == "csv") then
      real_format = "F15.10"
      row_format = "(" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // "A2" // trim(adjustl(real_format)) // ")"

      write(iunit, *) "# Plan Temperature Field : Z = ", pos_z
      write(iunit, *) "X,               Y,               T,"

      do j=1, ny
        pos_y = (j-0.5)*delta_y

        do i=1, nx
          pos_x = (i-0.5)*delta_x

          write(iunit,row_format) pos_x, ", ", pos_y, ", ", Temp(j + ny*(i-1) + nx*ny*(k_z-1))

        enddo
      write(iunit,*)
      enddo
    endif

  end subroutine write_Plan_Temp_field !}}}
#endif

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

  subroutine read_heatGenFile_perTimeStep(iunit, heat_surface_gen) !{{{
    integer, intent(in) :: iunit
    real(RP), dimension(:), allocatable, intent(inout) :: heat_surface_gen

    character(len=(STRING_LENGTH)) :: tmp
    integer :: nx_read, ny_read
    real(RP) :: time_read

    call read_header_heatGenFile(iunit, nx_read, ny_read, time_read)

    do i = 1, nx
      do j = 1, ny
        ! reading x, y
        read(iunit, "(A13)", advance="no") tmp
        ! reading q and skipping the end of the line :
        read(iunit, *) heat_surface_gen(j + ny*(i-1))
      enddo
    enddo

  end subroutine read_heatGenFile_perTimeStep !}}}

  subroutine read_header_heatGenFile(iunit, nx_read, ny_read, time_read) !{{{
    integer, intent(in) :: iunit
    integer, intent(inout) :: nx_read, ny_read
    real(RP), intent(inout) :: time_read

    character(len=(STRING_LENGTH)) :: tmp_str
    integer :: tmp_int

    ! Reading first two lines
    read(iunit, *) tmp_str
    read(iunit, *) tmp_str

    ! Read third line :
    read(iunit, "(A8)", advance="no") tmp_str
    read(iunit, "(I10)", advance="no") nx_read
    read(iunit, "(A3)", advance="no") tmp_str
    read(iunit, "(I10)", advance="no") ny_read
    read(iunit, "(A19)", advance="no") tmp_str
    read(iunit, "(I10)", advance="no") tmp_int
    read(iunit, "(A13)", advance="no") tmp_str
    read(iunit, *) time_read

    if(nx_read/=nx) then
      write(*,*) "Error : nx from heat generation file different from nx choosed in namelist"
      stop
    endif

    if(ny_read/=ny) then
      write(*,*) "Error : ny from heat generation file different from ny choosed in namelist"
      stop
    endif

  end subroutine read_header_heatGenFile !}}}

end module IO_mod
