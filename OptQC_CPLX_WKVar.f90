module arrays_cplx

implicit none

type arr_dp_4_cplx
    double complex, allocatable :: arr(:,:,:,:)
    integer :: d1, d2, d3, d4 = 0
contains
    procedure :: constructor => arr_dp_4_cplx_constructor
    procedure :: destructor => arr_dp_4_cplx_destructor
end type arr_dp_4_cplx

type l_arr_dp_4_cplx
    type(arr_dp_4_cplx), allocatable :: l(:)
    integer :: d = 0
contains
    procedure :: constructor => l_arr_dp_4_cplx_constructor
    procedure :: destructor => l_arr_dp_4_cplx_destructor
end type l_arr_dp_4_cplx

type arr_dp_2_cplx
    double complex, allocatable :: arr(:,:)
    integer :: d1, d2 = 0
contains
    procedure :: constructor => arr_dp_2_cplx_constructor
    procedure :: destructor => arr_dp_2_cplx_destructor
end type arr_dp_2_cplx

type l_arr_dp_2_cplx
    type(arr_dp_2_cplx), allocatable :: l(:)
    integer :: d = 0
contains
    procedure :: constructor => l_arr_dp_2_cplx_constructor
    procedure :: destructor => l_arr_dp_2_cplx_destructor
end type l_arr_dp_2_cplx

contains

subroutine arr_dp_4_cplx_constructor(this,d1,d2,d3,d4)

implicit none
class(arr_dp_4_cplx) :: this
integer :: d1, d2, d3, d4

allocate(this%arr(d1,d2,d3,d4))
this%d1 = d1
this%d2 = d2
this%d3 = d3
this%d4 = d4

end subroutine arr_dp_4_cplx_constructor

subroutine arr_dp_4_cplx_destructor(this)

implicit none
class(arr_dp_4_cplx) :: this

deallocate(this%arr)
this%d1 = 0
this%d2 = 0
this%d3 = 0
this%d4 = 0

end subroutine arr_dp_4_cplx_destructor

subroutine l_arr_dp_4_cplx_constructor(this,d)

implicit none
class(l_arr_dp_4_cplx) :: this
integer :: d

allocate(this%l(d))
this%d = d

end subroutine l_arr_dp_4_cplx_constructor

subroutine l_arr_dp_4_cplx_destructor(this)

implicit none
class(l_arr_dp_4_cplx) :: this
integer :: i

do i = 1, this%d
    call this%l(i)%destructor()
end do
deallocate(this%l)
this%d = 0

end subroutine l_arr_dp_4_cplx_destructor

subroutine arr_dp_2_cplx_constructor(this,d1,d2)

implicit none
class(arr_dp_2_cplx) :: this
integer :: d1, d2

allocate(this%arr(d1,d2))
this%d1 = d1
this%d2 = d2

end subroutine arr_dp_2_cplx_constructor

subroutine arr_dp_2_cplx_destructor(this)

implicit none
class(arr_dp_2_cplx) :: this

deallocate(this%arr)
this%d1 = 0
this%d2 = 0

end subroutine arr_dp_2_cplx_destructor

subroutine l_arr_dp_2_cplx_constructor(this,d)

implicit none
class(l_arr_dp_2_cplx) :: this
integer :: d

allocate(this%l(d))
this%d = d

end subroutine l_arr_dp_2_cplx_constructor

subroutine l_arr_dp_2_cplx_destructor(this)

implicit none
class(l_arr_dp_2_cplx) :: this
integer :: i

do i = 1, this%d
    call this%l(i)%destructor()
end do
deallocate(this%l)
this%d = 0

end subroutine l_arr_dp_2_cplx_destructor

end module arrays_cplx
