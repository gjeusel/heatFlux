module IO_mod
#include "config.h"

  use Parameters_mod
  use Common_Variables_mod

  implicit none
  private

  public :: opening_file
  !public :: write_Temp

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

  !subroutine write_Temp(Temp, i0, i1, j0, j1)
  !  real(RP), dimension(:), intent(in) :: Temp
  !  integer, intent(in) :: i0, i1, j0, j1

  !  !do i=i0, i1
  !  !  do j=j0, j1
  !  !    write(*,"")
  !  !  enddo
  !  !enddo

  !end subroutine write_Temp


end module IO_mod
