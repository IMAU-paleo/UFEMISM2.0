module string_module

  use mpi_basic, only: par

  implicit none

  private

  public :: separate_strings_by_double_vertical_bars

contains

  subroutine separate_strings_by_double_vertical_bars( str_list, strs)
    !< Take a list of items separated by double vertical bars ("||"),
    !< and return them as separate strings

    ! In/output variables:
    character(len=*),                            intent(in   ) :: str_list
    character(len=*), dimension(:), allocatable, intent(  out) :: strs

    ! Local variables:
    integer                           :: i, n_doublebars,ii
    character(len=len_trim(str_list)) :: str_list_redux

    ! Count number of instances of "||"
    n_doublebars = 0
    do i = 1, len_trim(str_list)
      if (str_list(i:i+1) == '||') n_doublebars = n_doublebars + 1
    end do

    allocate(strs(n_doublebars+1))

    str_list_redux = str_list
    do i = 1, n_doublebars
      ii = index( str_list_redux,'||')
      strs(i) = str_list_redux(1:ii-1)
      str_list_redux = str_list_redux(ii+2:len_trim(str_list_redux))
    end do
    strs(n_doublebars+1) = str_list_redux

  end subroutine separate_strings_by_double_vertical_bars

end module string_module