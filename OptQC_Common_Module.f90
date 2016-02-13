module common_module

implicit none

type cstr
    character, allocatable :: str(:)
    integer :: len
contains
    procedure :: constructor => cstr_constructor
    procedure :: destructor => cstr_destructor
    procedure :: copy => cstr_copy
end type cstr

type, extends (cstr) :: binstr
contains
	procedure :: getbinrep => binstr_getbinrep
    procedure :: qsimilar => binstr_qsimilar
end type binstr

type prog_args
    character(len=128) :: fbase
    integer :: flength
    integer :: PROG_TYPE
    integer :: ITER_LIM
    integer :: PERM_ITER_LIM
    double precision :: TOL_COEFF
contains
end type prog_args

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

subroutine cstr_constructor(this,n)

implicit none
class(cstr) :: this
integer :: n

allocate(this%str(n))
this%len = n

end subroutine cstr_constructor

subroutine cstr_destructor(this)

implicit none
class(cstr) :: this

deallocate(this%str)
this%len = 0

end subroutine cstr_destructor

subroutine cstr_copy(this,targ)

implicit none
class(cstr) :: this, targ
integer :: i, mlen

mlen = min(this%len,targ%len)
this%str(1:mlen) = targ%str(1:mlen)
if(mlen < this%len) then
    this%str(mlen+1:this%len) = ''
end if

end subroutine cstr_copy

subroutine binstr_getbinrep(this,num)

implicit none
class(binstr) :: this
integer :: num
integer :: i

this%str = '0'
do i = 1, this%len
	if(btest(num,i-1) == .true.) then
		this%str(this%len-i+1) = '1'
	end if
end do

end subroutine binstr_getbinrep

function binstr_qsimilar(this,targ,e_pos)

implicit none
class(binstr) :: this, targ
integer :: e_pos
integer :: binstr_qsimilar

integer :: i, p

! Automatically reject if bit stirngs are not of the same length
if(this%len /= targ%len) then
    binstr_qsimilar = 0
    return
end if
binstr_qsimilar = 1
! Must differ at the specified e_pos bit
p = e_pos
if(this%str(p) == targ%str(p)) then
    binstr_qsimilar = 0
    return
end if
! Must not differ at elsewhere
if(e_pos > 1) then
    do i = 1, e_pos-1
        p = i
        if(this%str(p) == targ%str(p)) then
            binstr_qsimilar = 0
            return
        end if
    end do
end if
if(e_pos < this%len-1) then
    do i = e_pos+1, this%len-1
        p = i
        if(this%str(p) == targ%str(p)) then
            binstr_qsimilar = 0
            return
        end if
    end do
end if
return

end function binstr_qsimilar

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
