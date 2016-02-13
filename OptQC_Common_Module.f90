module common_module

implicit none

type arr_dp_1
    double precision, allocatable :: arr(:)
    integer :: d1 = 0
contains
    procedure :: constructor => arr_dp_1_constructor
    procedure :: destructor => arr_dp_1_destructor
end type arr_dp_1

type l_arr_dp_1
    type(arr_dp_1), allocatable :: l(:)
    integer :: d = 0
contains
    procedure :: constructor => l_arr_dp_1_constructor
    procedure :: destructor => l_arr_dp_1_destructor
end type l_arr_dp_1

type arr_int_1
    integer, allocatable :: arr(:)
    integer :: d1 = 0
contains
    procedure :: constructor => arr_int_1_constructor
    procedure :: destructor => arr_int_1_destructor
end type arr_int_1

type l_arr_int_1
    type(arr_int_1), allocatable :: l(:)
    integer :: d = 0
contains
    procedure :: constructor => l_arr_int_1_constructor
    procedure :: destructor => l_arr_int_1_destructor
end type l_arr_int_1

contains

subroutine arr_dp_1_constructor(this,d1)

implicit none
class(arr_dp_1) :: this
integer :: d1

allocate(this%arr(d1))
this%d1 = d1

end subroutine arr_dp_1_constructor

subroutine arr_dp_1_destructor(this)

implicit none
class(arr_dp_1) :: this

deallocate(this%arr)
this%d1 = 0

end subroutine arr_dp_1_destructor

subroutine l_arr_dp_1_constructor(this,d)

implicit none
class(l_arr_dp_1) :: this
integer :: d

allocate(this%l(d))
this%d = d

end subroutine l_arr_dp_1_constructor

subroutine l_arr_dp_1_destructor(this)

implicit none
class(l_arr_dp_1) :: this
integer :: i

do i = 1, this%d
    call this%l(i)%destructor()
end do
deallocate(this%l)
this%d = 0

end subroutine l_arr_dp_1_destructor

subroutine arr_int_1_constructor(this,d1)

implicit none
class(arr_int_1) :: this
integer :: d1

allocate(this%arr(d1))
this%d1 = d1

end subroutine arr_int_1_constructor

subroutine arr_int_1_destructor(this)

implicit none
class(arr_int_1) :: this

deallocate(this%arr)
this%d1 = 0

end subroutine arr_int_1_destructor

subroutine l_arr_int_1_constructor(this,d)

implicit none
class(l_arr_int_1) :: this
integer :: d

allocate(this%l(d))
this%d = d

end subroutine l_arr_int_1_constructor

subroutine l_arr_int_1_destructor(this)

implicit none
class(l_arr_int_1) :: this
integer :: i

do i = 1, this%d
    call this%l(i)%destructor()
end do
deallocate(this%l)
this%d = 0

end subroutine l_arr_int_1_destructor

end module common_module
