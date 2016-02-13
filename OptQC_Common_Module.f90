module common_module

implicit none

type prog_args
    character(len=128) :: fbase
    integer :: flength
    integer :: PROG_TYPE
    integer :: ITER_LIM
    integer :: PERM_ITER_LIM
    double precision :: TOL_COEFF
contains
end type prog_args

type csd_write_handle
    ! Fixed global parameters
    integer :: NVAR
    integer :: ROWMAX, COLMAX, SECMAX
    ! File stream details
    character(len=128) :: fname
    integer :: FN
    ! Pointer to circuit
    character(len=32), pointer :: Circuit(:,:)
    integer :: nct                                          ! Length of circuit
    ! Output variables
    character(len=256), allocatable :: buffer(:)
    integer :: cpos
    integer :: buffer_num, buffer_col
    integer :: page_num
    logical :: env_open
    logical :: neg
    logical :: seperator
contains
    procedure :: constructor => csd_write_handle_constructor
    procedure :: destructor => csd_write_handle_destructor
    procedure :: clean => csd_write_handle_clean
    procedure :: set_file => csd_write_handle_set_file
    procedure :: assign_target => csd_write_handle_assign_target
    procedure :: nullify_target => csd_write_handle_nullify_target
    procedure :: preamble => csd_write_handle_preamble
    procedure :: postamble => csd_write_handle_postamble
    procedure :: insert_neg => csd_write_handle_insert_neg
    procedure :: insert_seperator => csd_write_handle_insert_seperator
    procedure :: push_back => csd_write_handle_push_back
    procedure :: write_buffer => csd_write_handle_write_buffer
    procedure :: write_circuit => csd_write_handle_write_circuit
end type csd_write_handle

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

subroutine csd_write_handle_constructor(this,N)

implicit none
class(csd_write_handle) :: this
integer :: N

this%NVAR = N+1
this%fname = ""
this%FN = -1
this%ROWMAX = 36
this%COLMAX = 10
this%SECMAX = 36 / this%NVAR
allocate(this%buffer(this%NVAR))
call this%clean()

end subroutine csd_write_handle_constructor

subroutine csd_write_handle_destructor(this)

implicit none
class(csd_write_handle) :: this

deallocate(this%buffer)

end subroutine csd_write_handle_destructor

subroutine csd_write_handle_clean(this)

implicit none
class(csd_write_handle) :: this

! Does not reset the file assigned
this%buffer = ""
this%cpos = 0
this%buffer_num = 0
this%buffer_col = 0
this%page_num = 0
this%env_open = .false.
this%neg = .false.
this%seperator = .false.

end subroutine csd_write_handle_clean

subroutine csd_write_handle_set_file(this,fname)

implicit none
class(csd_write_handle) :: this
character(len=128) :: fname

integer :: RINT

this%fname = fname
! Select a random unit number (10 or above)
this%FN = RINT(1000)
do while (this%FN < 10)
    this%FN = RINT(1000)
end do

end subroutine csd_write_handle_set_file

subroutine csd_write_handle_assign_target(this,Circuit,csdr_ct)

implicit none
class(csd_write_handle) :: this
character(len=32), target :: Circuit(:,:)
integer :: csdr_ct

this%Circuit => Circuit
this%nct = csdr_ct

end subroutine csd_write_handle_assign_target

subroutine csd_write_handle_nullify_target(this)

implicit none
class(csd_write_handle) :: this

nullify(this%Circuit)
this%nct = 0

end subroutine csd_write_handle_nullify_target

subroutine csd_write_handle_preamble(this)

implicit none
class(csd_write_handle) :: this

open(unit=this%FN,file=this%fname,status='replace')
if(this%FN /= -1) then
    write(this%FN,'(a)')"\documentclass{article}"
    write(this%FN,'(a)')"\usepackage[cm]{fullpage}"
    write(this%FN,'(a)')"\input{Qcircuit}"
    write(this%FN,'(a)')""
    write(this%FN,'(a)')"\begin{document}"
    write(this%FN,'(a)')""
    write(this%FN,'(a)')"\Qcircuit @C=2.0em @R=0.1em @!R{"
    this%env_open = .true.
else
    write(*,'(a)')"Attempted to write preamble to unopened .tex file. Exiting program."
    stop
end if

end subroutine csd_write_handle_preamble

subroutine csd_write_handle_postamble(this)

implicit none
class(csd_write_handle) :: this

if(this%FN /= -1) then
    if(this%buffer_col > 0) then
        call this%write_buffer()
    end if
    if(this%env_open == .true.) then
        write(this%FN,'(a)')"}"
    end if
    write(this%FN,'(a)')""
    write(this%FN,'(a)')"\end{document}"
else
    write(*,'(a)')"Attempted to write postamble to unopened .tex file. Exiting program."
    stop
end if
close(this%FN)

end subroutine csd_write_handle_postamble

subroutine csd_write_handle_insert_neg(this)

implicit none
class(csd_write_handle) :: this

! UGH inelegant fudge into push_back subroutine.
this%neg = .true.
call this%push_back()

end subroutine csd_write_handle_insert_neg

subroutine csd_write_handle_insert_seperator(this)

implicit none
class(csd_write_handle) :: this

! Let's try classifying seperators as not counting towards the buffer_col count.
! It shouldn't take up that much (figure) space right....?
! WRONG - it does; introduced a new logical variable seperator to take care of this...
this%seperator = .true.
call this%push_back()

end subroutine csd_write_handle_insert_seperator

subroutine csd_write_handle_push_back(this)

implicit none
class(csd_write_handle) :: this

integer :: row

! Check if buffer is full first
if(this%buffer_col == this%COLMAX) then
    call this%write_buffer()
    this%buffer_col = 0
    this%buffer_num = this%buffer_num + 1
    ! Close current page if maximum number of sections has been reached
    if(this%buffer_num == this%SECMAX) then
        write(this%FN,'(a)')'}'
        this%env_open = .false.
        this%buffer_num = 0
        this%page_num = this%page_num + 1
    end if
end if
if(this%seperator == .true.) then
    do row = 1, this%NVAR
        if(row /= this%NVAR) then
            this%buffer(row) = trim(adjustl(this%buffer(row))) // " & | \qw"
        else
            this%buffer(row) = trim(adjustl(this%buffer(row))) // " & |"
        end if
    end do
    this%seperator = .false.
else if(this%neg == .true.) then
    this%buffer(1) = trim(adjustl(this%buffer(1))) // " & \gate{-I_2}"
    do row = 2, this%NVAR
        if(row /= this%NVAR) then
            this%buffer(row) = trim(adjustl(this%buffer(row))) // " & \qw"
        else
            this%buffer(row) = trim(adjustl(this%buffer(row))) // " & "
        end if
    end do
    this%neg = .false.
else
    do row = 1, this%NVAR
        this%buffer(row) = trim(adjustl(this%buffer(row))) // " " // trim(adjustl(this%Circuit(row,this%cpos)))
    end do
end if
this%buffer_col = this%buffer_col + 1

end subroutine csd_write_handle_push_back

subroutine csd_write_handle_write_buffer(this)

implicit none
class(csd_write_handle) :: this

integer :: row

! Check to see if the circuit envionment is presently open - opens one if not
if(this%env_open == .false.) then
    ! Forces page break after each Qcircuit - assume that preamble has been run correctly
    write(this%FN,'(a)')""
    write(this%FN,'(a)')"\clearpage"
    write(this%FN,'(a)')""
    write(this%FN,'(a)')'\Qcircuit @C=2.0em @R=0.1em @!R{'
    this%env_open = .true.
end if
! Writes termination characters to end of buffer
do row = 1, this%NVAR
    if(row /= this%NVAR) then
        this%buffer(row) = trim(adjustl(this%buffer(row))) // ' & \qw \\'
    else
        this%buffer(row) = trim(adjustl(this%buffer(row))) // ' & \\'
    end if
end do
! Write buffer to file
do row = 1, this%NVAR
    write(this%FN,'(a)')trim(adjustl(this%buffer(row)))
end do
! Empty buffer
this%buffer = ""

end subroutine csd_write_handle_write_buffer

subroutine csd_write_handle_write_circuit(this)

implicit none
class(csd_write_handle) :: this

integer :: idx

do idx = 1, this%nct
    this%cpos = idx
    call this%push_back()
end do

end subroutine csd_write_handle_write_circuit

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
