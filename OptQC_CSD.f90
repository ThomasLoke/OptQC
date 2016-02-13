module csd_real

use arrays_real
use common_module

implicit none

double precision, parameter :: CUTOFF = 0.0000000001d0
double precision, parameter :: PI = 3.1415926535897932384626433832795028841971693993751d0

type csd_solution
    integer :: N, M
    double precision, allocatable :: X(:,:)
    character(len=15), allocatable :: Circuit(:,:)          ! Note: By design, this is empty unless run_csdr is called
    integer :: csd_ct, csdr_ct
contains
    procedure :: constructor => csd_solution_constructor
    procedure :: destructor => csd_solution_destructor
    procedure :: clean => csd_solution_clean
    procedure :: copy => csd_solution_copy
    procedure :: run_csd => csd_solution_run_csd
    procedure :: run_csdr => csd_solution_run_csdr
    procedure :: write_circuit => csd_solution_write_circuit
end type csd_solution

type csd_solution_set
    integer :: nset
    integer :: N, M
    type(csd_solution), allocatable :: arr(:)
    integer :: csd_ss_ct, csdr_ss_ct
contains
    procedure :: constructor => csd_solution_set_constructor
    procedure :: destructor => csd_solution_set_destructor
    procedure :: clean => csd_solution_set_clean
    procedure :: copy => csd_solution_set_copy
    procedure :: run_csd => csd_solution_set_run_csd
    procedure :: run_csdr => csd_solution_set_run_csdr
    procedure :: write_circuit => csd_solution_set_write_circuit
end type csd_solution_set

type csd_generator
    ! Includes standard CSD and the reduction procedures
    ! Global variables throughout all subroutines
    integer :: N, M, Mh                                     ! Reminder: M = 2**N, Mh = M/2
    integer, pointer :: index_level(:), index_pair(:,:,:)   ! Remains constant for matrices of the same dimensions
    double precision, allocatable :: GATEY(:,:), GATEPI(:,:)! Intermediate variables carrying CSD results - moved from csd_solution object
    ! Note: Could just store a pointer to the csd_solution object instead - but done this way for less indirection (i.e. more clarity)
    double precision, pointer :: X(:,:)                     ! Intent: In
    character(len=15), pointer :: Circuit(:,:)              ! Intent: Out
    integer, pointer :: csd_ct, csdr_ct                     ! Intent: Out
    ! Workspace arrays
    ! CYGR_CUTGATE:
    integer, allocatable :: Z_array(:,:), PI_array(:,:), GATEY_sign(:)
    ! CYGR_CSD:
    type(l_arr_dp_4) :: Z0, Z1
    ! CYGR_BLKCSD:
    type(l_arr_dp_2) :: X_blk, X11, X12, X21, X22, U1, U2, V1T, V2T
    type(l_arr_dp_1) :: GATEY_blk
    type(l_arr_int_1) :: IWORK
    ! ReduceSolution:
    character(len=20), allocatable :: C_Num_Bin(:,:)
    double precision, allocatable :: Type_Param(:)
    integer, allocatable :: N_Per_Type(:)
    ! Scratch variables between CYGR_CSD and CYGR_BLKCSD
    integer :: size0, size1, num0, num1
    ! Scratch variables between ReduceSolution and (GroupGates,ReduceGroups)
    integer :: N_Type, N_Total
contains
    procedure :: constructor => csd_generator_constructor
    procedure :: destructor => csd_generator_destructor
    procedure :: assign_target => csd_generator_assign_target
    procedure :: nullify_target => csd_generator_nullify_target
    procedure :: run_csd => csd_generator_run_csd
    procedure :: run_blkcsd => csd_generator_run_blkcsd
    procedure :: run_cutgate => csd_generator_run_cutgate
    procedure :: GateCount => csd_generator_GateCount
    procedure :: ReduceSolution => csd_generator_ReduceSolution
    procedure :: GroupGates => csd_generator_GroupGates
    procedure :: ReduceGroups => csd_generator_ReduceGroups
end type csd_generator

contains

subroutine csd_solution_constructor(this,N,M)

implicit none
class(csd_solution) :: this
integer :: N, M

this%N = N
this%M = M
this%csd_ct = 0
this%csdr_ct = 0
allocate(this%X(M,M))
allocate(this%Circuit(N+1,M*M/2+M))
this%X = 0.0d0
this%Circuit = ""

end subroutine csd_solution_constructor

subroutine csd_solution_destructor(this)

implicit none
class(csd_solution) :: this

deallocate(this%X)
deallocate(this%Circuit)

end subroutine csd_solution_destructor

subroutine csd_solution_clean(this)

implicit none
class(csd_solution) :: this

this%Circuit = ""
this%csd_ct = 0
this%csdr_ct = 0

end subroutine csd_solution_clean

subroutine csd_solution_copy(this,source)

implicit none
class(csd_solution) :: this
type(csd_solution) :: source

this%X = source%X
this%Circuit = source%Circuit
this%csd_ct = source%csd_ct
this%csdr_ct = source%csdr_ct

end subroutine csd_solution_copy

subroutine csd_solution_run_csd(this,generator)

implicit none
class(csd_solution) :: this
type(csd_generator) :: generator

call this%clean()
call generator%assign_target(this)
call generator%run_csd()                         ! Note: No circuit result
call generator%nullify_target()

end subroutine csd_solution_run_csd

subroutine csd_solution_run_csdr(this,generator)

implicit none
class(csd_solution) :: this
type(csd_generator) :: generator

call this%clean()
call generator%assign_target(this)
call generator%run_csd()
call generator%ReduceSolution()                  ! Note: Circuit result produced here
call generator%nullify_target()

end subroutine csd_solution_run_csdr

subroutine csd_solution_write_circuit(this,fout)

implicit none
class(csd_solution) :: this
character(len=128) :: fout

integer :: nvar, nlfull, colnum, colsum, partmax, partct, part, row, col

if(this%csdr_ct == 0) then
    write(*,'(a)')"Attempted to output empty circuit to .tex file - request ignored."
    return
end if
nvar = this%N + 1                               ! Requires N+1 rows, due to the last line being reserved for parameter values
partmax = 36 / nvar                             ! Maximum number of rows on one page set as 36, division gives the maximum number of circuit lines
open(unit=1,file=fout,status='replace')
write(1,'(a)') '\documentclass{amsart}'
write(1,'(a)') '\usepackage[matrix,frame,arrow]{xypic}'
write(1,'(a)') '\usepackage[cm]{fullpage}'
write(1,'(a)') '\input{Qcircuit}'
write(1,*)
write(1,'(a)') '\begin{document}'
colnum = 10                                     ! This determines the maximum number of columns in a line
colsum = this%csdr_ct                           ! Total number of columns - obtained from decomposition results
partct = 0
nlfull = colsum / colnum                        ! Number of full lines
do part = 1, nlfull                             ! Main looping section - colsum/colnum rounds down
    ! Start of circuit (or new page)
    if(partct == 0) then
        write(1,*)
        write(1,'(a)') '\['
        write(1,'(a)') '\Qcircuit @C=2.0em @R=0.1em @!R{'
    end if
	do row = 1,nvar
		do col = (part-1)*colnum+1,part*colnum
			write(1,303,advance='no') this%Circuit(row,col)
		end do
		if (row /= nvar) then
			write(1,303,advance='no') '& \qw \\'
		else
			write(1,303,advance='no') '&     \\'
		end if
		write(1,*)
	end do
	partct = partct + 1
	! If maximum number of lines is reached, skip to next page
	if(partct == partmax) then
        write(1,'(a)') '}'
        write(1,'(a)') '\]'
        partct = 0
    end if
end do
part = nlfull
if(part*colnum+1 <= colsum) then
    if(partct == 0) then
        write(1,*)
        write(1,'(a)') '\['
        write(1,'(a)') '\Qcircuit @C=2.0em @R=0.1em @!R{'
    end if
    do row = 1,nvar
        do col = part*colnum+1,colsum
            write(1,303,advance='no') this%Circuit(row,col)
        end do
        if (row /= nvar) then
            write(1,303,advance='no') '& \qw \\'
        else
            write(1,303,advance='no') '&     \\'
        end if
        write(1,*)
    end do
end if
write(1,'(a)') '}'
write(1,'(a)') '\]'
write(1,*)
write(1,'(a)') '\end{document}'
CLOSE(1)
303 format(A15)

end subroutine csd_solution_write_circuit

subroutine csd_solution_set_constructor(this,nset,N,M)

implicit none
class(csd_solution_set) :: this
integer :: nset, N, M

integer :: i

this%nset = nset
this%N = N
this%M = M
allocate(this%arr(nset))
do i = 1, nset
    call this%arr(i)%constructor(N,M)
end do
this%csd_ss_ct = 0
this%csdr_ss_ct = 0

end subroutine csd_solution_set_constructor

subroutine csd_solution_set_destructor(this)

implicit none
class(csd_solution_set) :: this

integer :: i

do i = 1, this%nset
    call this%arr(i)%destructor()
end do
deallocate(this%arr)

end subroutine csd_solution_set_destructor

subroutine csd_solution_set_clean(this)

implicit none
class(csd_solution_set) :: this

! Note: Cleaning of csd_solution objects is done individually in their call to run_csd or run_csdr
this%csd_ss_ct = 0
this%csdr_ss_ct = 0

end subroutine csd_solution_set_clean

subroutine csd_solution_set_copy(this,source)

implicit none
class(csd_solution_set) :: this
type(csd_solution_set) :: source

integer :: i

do i = 1, this%nset
    call this%arr(i)%copy(source%arr(i))
end do
this%csd_ss_ct = source%csd_ss_ct
this%csdr_ss_ct = source%csdr_ss_ct

end subroutine csd_solution_set_copy

subroutine csd_solution_set_run_csd(this,generator)

implicit none
class(csd_solution_set) :: this
type(csd_generator) :: generator

integer :: i

call this%clean()
do i = 1, this%nset
    call this%arr(i)%run_csd(generator)
    this%csd_ss_ct = this%csd_ss_ct + this%arr(i)%csd_ct
end do

end subroutine csd_solution_set_run_csd

subroutine csd_solution_set_run_csdr(this,generator)

implicit none
class(csd_solution_set) :: this
type(csd_generator) :: generator

integer :: i

call this%clean()
do i = 1, this%nset
    call this%arr(i)%run_csdr(generator)
    this%csd_ss_ct = this%csd_ss_ct + this%arr(i)%csd_ct
    this%csdr_ss_ct = this%csdr_ss_ct + this%arr(i)%csdr_ct
end do

end subroutine csd_solution_set_run_csdr

subroutine csd_solution_set_write_circuit(this,fout)

implicit none
class(csd_solution_set) :: this
character(len=128) :: fout

! TODO: OVERHAUL EVERYTHING RAWR
if(this%nset == 3) then
    call this%arr(2)%write_circuit(fout)
end if

end subroutine csd_solution_set_write_circuit

subroutine csd_generator_constructor(this,N,M,index_level,index_pair)

implicit none
class(csd_generator) :: this
integer :: N, M
integer, target :: index_level(M-1), index_pair(M/2,2,N)
integer :: i, Mh, size0, size1, num0, num1

Mh = M / 2
this%N = N
this%M = M
this%Mh = Mh
this%index_level => index_level
this%index_pair => index_pair
allocate(this%GATEY(M/2,M-1))
allocate(this%GATEPI(M/2,N))
call this%nullify_target()
allocate(this%Z_array(M,M))
allocate(this%PI_array(M,N))
allocate(this%GATEY_sign(Mh))
call this%Z0%constructor(N)
call this%Z1%constructor(N)
call this%X_blk%constructor(N)
call this%X11%constructor(N)
call this%X12%constructor(N)
call this%X21%constructor(N)
call this%X22%constructor(N)
call this%U1%constructor(N)
call this%U2%constructor(N)
call this%V1T%constructor(N)
call this%V2T%constructor(N)
call this%GATEY_blk%constructor(N)
call this%IWORK%constructor(N)
do i = 1, N
    size1 = 2**(N-i)
    size0 = size1 * 2
    num0 = 2**(i-1)
    num1 = num0 * 2
    call this%Z0%l(i)%constructor(size1,size1,num1,num1)
    call this%Z1%l(i)%constructor(size1,size1,num1,num1)
    call this%X_blk%l(i)%constructor(size0,size0)
    call this%X11%l(i)%constructor(size1,size1)
    call this%X12%l(i)%constructor(size1,size1)
    call this%X21%l(i)%constructor(size1,size1)
    call this%X22%l(i)%constructor(size1,size1)
    call this%U1%l(i)%constructor(size1,size1)
    call this%U2%l(i)%constructor(size1,size1)
    call this%V1T%l(i)%constructor(size1,size1)
    call this%V2T%l(i)%constructor(size1,size1)
    call this%GATEY_blk%l(i)%constructor(size1)
    call this%IWORK%l(i)%constructor(size1)
end do
allocate(this%C_Num_Bin(Mh,Mh))
allocate(this%Type_Param(Mh))
allocate(this%N_Per_Type(Mh))
this%size0 = 0
this%size1 = 0
this%num0 = 0
this%num1 = 0

end subroutine csd_generator_constructor

subroutine csd_generator_destructor(this)

implicit none
class(csd_generator) :: this

nullify(this%index_level)
nullify(this%index_pair)
deallocate(this%GATEY)
deallocate(this%GATEPI)
call this%nullify_target()
deallocate(this%Z_array,this%PI_array,this%GATEY_sign)
call this%Z0%destructor()
call this%Z1%destructor()
call this%X_blk%destructor()
call this%X11%destructor()
call this%X12%destructor()
call this%X21%destructor()
call this%X22%destructor()
call this%U1%destructor()
call this%U2%destructor()
call this%V1T%destructor()
call this%V2T%destructor()
call this%GATEY_blk%destructor()
call this%IWORK%destructor()
deallocate(this%C_Num_Bin)
deallocate(this%Type_Param)
deallocate(this%N_Per_Type)

end subroutine csd_generator_destructor

subroutine csd_generator_assign_target(this,csd_ss)

implicit none
class(csd_generator) :: this
type(csd_solution), target :: csd_ss

! Bind pointers to the csd_solution elements
this%X => csd_ss%X
this%Circuit => csd_ss%Circuit
this%Circuit = ''
this%csd_ct => csd_ss%csd_ct
this%csd_ct = 0
this%csdr_ct => csd_ss%csdr_ct
this%csdr_ct = 0

end subroutine csd_generator_assign_target

subroutine csd_generator_nullify_target(this)

implicit none
class(csd_generator) :: this

nullify(this%X)
nullify(this%Circuit)
nullify(this%csd_ct)
nullify(this%csdr_ct)

end subroutine csd_generator_nullify_target

subroutine csd_generator_run_csd(this)

implicit none
class(csd_generator) :: this

integer :: N, i, j

N = this%N
! Clean up intermediate variables first
this%GATEY = 0.0d0
this%GATEPI = 0.0d0

!	Cosine Sine Decomposition (CSD) of X recursively until the Nth level
!	The 1st level
this%size0 = this%M
this%size1 = this%Mh
this%num0 = 1
this%num1 = 2
this%Z1%l(1)%arr = 0.0d0
!index_Y = this%Mh
call this%run_blkcsd(1,1)
this%Z0%l(1)%arr = this%Z1%l(1)%arr

!	The ist level (i = 2,...,N)
do i = 2, N
    this%size1 = 2**(N-i)
    this%size0 = this%size1 * 2
    this%num0 = 2**(i-1)
    this%num1 = this%num0 * 2
    this%Z1%l(i)%arr = 0.0d0
    do j = 1, this%num0
        call this%run_blkcsd(i,j)
    end do
    this%Z0%l(i)%arr = this%Z1%l(i)%arr
end do

!	Step 3: CYGC_CSD processes the CSD results
!	        In this step, we obtain all the Ry Gates and the PI Gates for X: GATEY, GATEPI
call this%run_cutgate()

!   Step 4 (added): Obtain count of gates
call this%GateCount()

end subroutine csd_generator_run_csd

subroutine csd_generator_run_blkcsd(this,i,j)

implicit none
class(csd_generator), target :: this
integer :: i, j

integer :: index_Y, size0, size1, num0, num1
double precision, pointer :: X(:,:,:), U(:,:,:), V(:,:,:), GATEY(:), X_blk(:,:), X11(:,:), X12(:,:), X21(:,:), X22(:,:), U1(:,:), U2(:,:), V1T(:,:), V2T(:,:), GATEY_blk(:)
integer, pointer :: IWORK(:)
double precision, allocatable :: WORK(:)
integer :: LWORK, INFO
integer :: k, k1, k2

if(i==1) then
    index_Y = this%Mh
else
    index_Y = 2**(this%N-i)*(2*j-1)
end if
size0 = this%size0
size1 = this%size1
num0 = this%num0
num1 = this%num1

! Initialize pointers to elements of the csd_generator object
if(i==1) then
    nullify(X)
else
    X => this%Z0%l(i-1)%arr(:,:,:,j)
end if
U => this%Z1%l(i)%arr(:,:,:,2*j-1)
V => this%Z1%l(i)%arr(:,:,:,2*j)
GATEY => this%GATEY(:,index_Y)
X_blk => this%X_blk%l(i)%arr
X11 => this%X11%l(i)%arr
X12 => this%X12%l(i)%arr
X21 => this%X21%l(i)%arr
X22 => this%X22%l(i)%arr
U1 => this%U1%l(i)%arr
U2 => this%U2%l(i)%arr
V1T => this%V1T%l(i)%arr
V2T => this%V2T%l(i)%arr
GATEY_blk => this%GATEY_blk%l(i)%arr
IWORK => this%IWORK%l(i)%arr

X_blk = 0.0d0
X11 = 0.0d0
X12 = 0.0d0
X21 = 0.0d0
X22 = 0.0d0
U1 = 0.0d0
U2 = 0.0d0
V1T = 0.0d0
V2T = 0.0d0
GATEY_blk = 0.0d0
IWORK = 0

do k = 1, num0
	if(i==1) then
        X_blk = this%X(:,:)
    else
        X_blk = X(:,:,k)
    end if
	do k1 = 1,size0
		do k2 = 1,size0
			if(abs(X_blk(k1,k2)) <= CUTOFF) then
				X_blk(k1,k2) = 0.0
			end if
		end do
	end do
	X11 = X_blk(1:size1,1:size1)
	X12 = X_blk(1:size1,size1+1:size0)
	X21 = X_blk(size1+1:size0,1:size1)
	X22 = X_blk(size1+1:size0,size1+1:size0)
	LWORK = -1
	allocate(WORK(1))
	WORK = 0
    !	Set LWORK = -1, then a workspace query is assumed; the subroutine only calculates the optimal size of the WORK array, returns this value as the first entry of the work array.
	call DORCSD('Y', 'Y', 'Y', 'Y', 'N', 'O', size0, size1, size1, X11,&
size1, X12, size1, X21, size1, X22, size1, GATEY_blk, U1, size1, U2, &
size1, V1T, size1, V2T, size1, WORK, LWORK, IWORK, INFO)
	!	Now set LWORK = WORK(1)
	LWORK = WORK(1)
	deallocate(WORK)
	allocate(WORK(LWORK))
	WORK = 0
	call DORCSD('Y', 'Y', 'Y', 'Y', 'N', 'O', size0, size1, size1, X11, &
size1, X12, size1, X21, size1, X22, size1, GATEY_blk, U1, size1, U2, &
size1, V1T, size1, V2T, size1, WORK, LWORK, IWORK, INFO)
	deallocate(WORK)
	U(:,:,2*k-1) = U1
	U(:,:,2*k) = U2
	V(:,:,2*k-1) = V1T
	V(:,:,2*k) = V2T
	GATEY(size1*(k-1)+1:size1*k) = GATEY_blk
end do

GATEY = GATEY/PI
return

end subroutine csd_generator_run_blkcsd

subroutine csd_generator_run_cutgate(this)

implicit none
class(csd_generator), target :: this

integer :: N, M
double precision, pointer :: Z(:,:)
integer :: i, j, leng

N = this%N
M = this%M

! Initialize pointers to elements of the csd_generator object
Z => this%Z0%l(N)%arr(1,1,:,:)

!	Z(j) = diag(Z_array(j)) (j = 1,2,...,2**N), and Z(j) or Z_array(j) contains only 1 or -1.
!	In this SUBROUTINE, to make it short, the input Z(j) = Z_array(j), namely, Z(j) is a M-length array, and Z is a 2D table with the dimension
!	(M,M+N)
!	e.g. Z(j) = [1,-1,-1,1,-1,1,-1,-1]
this%Z_array = NINT(Z)
do i = 2,M
	this%Z_array(:,i) = this%Z_array(:,i-1)*this%Z_array(:,i)
end do

!	Decide the Ry Gates "GATEY"
this%GATEY_sign = 0
do i = 1,M-1
    j = this%index_level(i)
    this%GATEY_sign = this%Z_array(this%index_pair(:,1,j),i)*this%Z_array(this%index_pair(:,2,j),i)
	this%GATEY(:,i) = this%GATEY_sign*this%GATEY(:,i)
end do

!	Decide the PI Gates "GATEPI"
this%PI_array(:,:) = 0
this%PI_array(:,N) = this%Z_array(:,M)
do i = 1,N-1
	leng = 2**i
	do j = 1, 2**(N-i)
		if (this%PI_array(leng*(j-1)+1,N+1-i) == -1) then
			this%PI_array(leng*(j-1)+1:leng*j,N+1-i) = -this%PI_array(leng*(j-1)+1:leng*j,N+1-i)
			this%PI_array(leng*(j-1)+1:leng*j,N-i) = -1
		else
			this%PI_array(leng*(j-1)+1:leng*j,N-i) = 1
		end if
	end do
end do
if (this%PI_array(1,1) == -1) then
	this%PI_array(:,1) = -this%PI_array(:,1)
end if

do i = 1,N
	leng = 2**(N+1-i)
	do j = 1,2**(i-1)
		this%GATEPI(j,i) = this%PI_array(leng*j,i)
	end do
end do
return

end subroutine csd_generator_run_cutgate

subroutine csd_generator_GateCount(this)

implicit none
class(csd_generator) :: this

integer :: i, j, res

res = 0
do i = 1, this%N
    do j = 1, 2**(i-1)
		if (abs(this%GATEPI(j,i)+1.0d0) <= CUTOFF) then
			res = res + 1
		end if
	end do
end do
do i = this%M-1, 1, -1
    do j = this%Mh, 1, -1
		if (abs(this%GATEY(j,i)) > CUTOFF) then
			res = res + 1
		end if
	end do
end do
this%csd_ct = res
return

end subroutine csd_generator_GateCount


subroutine csd_generator_ReduceSolution(this)

implicit none
class(csd_generator) :: this
integer :: N, ct
integer :: i, j, k, p, ref, row, lb, rb
character(len=20) :: control, workstr
integer :: IsEmpty

! Includes update of the count of the reduced number of gates at the end
N = this%N
ct = 0
this%Circuit = ''
! Form the circuit for the reduced solution.
! GATEPI
ct = 0
do i = 1, N
    lb = 22-i
    call this%GroupGates(2**(i-1),this%GATEPI(:,i),-1.0d0)
    call this%ReduceGroups()
    ! There should only be one type, that is, Type_Param(1) = -1.0d0
    if(this%N_Type > 1) then
        write(*,*)"Error! GATEPI value different to -1.0d0 encountered. GATEPI section:"
        write(*,*)this%GATEPI(:,i)
        write(*,*)"Terminating program."
        call exit(1)
    end if
    do k = 1, this%N_Per_Type(1)
        ct = ct+1
        if(i > 1) then
            workstr = this%C_Num_Bin(1,k)
            do row = 1, i-1
                p = 21-i+row
                if(workstr(p:p)=='0') then
                    this%Circuit(row,ct) = '& \ctrlo{1}'
                else if(workstr(p:p)=='1') then
                    this%Circuit(row,ct) = '& \ctrl{1}'
                else
                    ! Case where workstr(p:p) = '*'
                    if(row > 1 .and. IsEmpty(workstr,lb,p-1)==0) then
                        this%Circuit(row,ct) = '& \qw \qwx[1]'
                    else
                        this%Circuit(row,ct) = '& \qw'
                    end if
                end if
            end do
        end if
        this%Circuit(i,ct) = '& \gate{\pi}'
        if(i < N) then
            do row = i+1,N
                this%Circuit(row,ct) = '& \qw'
            end do
        end if
        this%Circuit(N+1,ct) = '&'
    end do
end do
!GATEY
lb = 22-N
rb = 20
do i = this%M-1, 1, -1
    ref = this%index_level(i)
    call this%GroupGates(this%Mh,this%GATEY(:,i),0.0d0)
    call this%ReduceGroups()
    do j = 1, this%N_Type
        do k = 1, this%N_Per_Type(j)
            ct = ct+1
            workstr = this%C_Num_Bin(j,k)
            this%Circuit(ref,ct) = '& \gate{R_y}'
            do row = 1,N-1
                p = 21-N+row
				if(row < ref) then
					if (workstr(p:p)=='0') then
						this%Circuit(row,ct) = '& \ctrlo{1}'
					else if(workstr(p:p)=='1') then
						this%Circuit(row,ct) = '& \ctrl{1}'
                    else
                        ! Case where workstr(p:p) = '*'
                        if(row > 1 .and. IsEmpty(workstr,lb,p-1)==0) then
                            this%Circuit(row,ct) = '& \qw \qwx[1]'
                        else
                            this%Circuit(row,ct) = '& \qw'
                        end if
					end if
				else
					if(workstr(p:p)=='0') then
						this%Circuit(row+1,ct) = '& \ctrlo{-1}'
					else if(workstr(p:p)=='1') then
						this%Circuit(row+1,ct) = '& \ctrl{-1}'
                    else
                        ! Case where workstr(p:p) = '*'
                        if(row < N-1 .and. IsEmpty(workstr,p+1,rb)==0) then
                            this%Circuit(row+1,ct) = '& \qw \qwx[-1]'
                        else
                            this%Circuit(row+1,ct) = '& \qw'
                        end if
					end if
				end if
			end do
			write(this%Circuit(N+1,ct),302)this%Type_Param(j)
        end do
    end do
end do
302 format('& ',F8.4)
this%csdr_ct = ct

end subroutine csd_generator_ReduceSolution

subroutine csd_generator_GroupGates(this,extent,Gate_Param,Sig_Offset)

implicit none
class(csd_generator) :: this
integer :: extent ! Length of Gate_Param section - variable for GATEPI
double precision :: Gate_Param(:)
double precision :: Sig_Offset ! -1 for GatePI, 0 for GATEY

double precision :: ref
integer :: i, j, match_idx

this%N_Per_Type = 0
this%N_Type = 0
this%N_Total = 0
do i = 1, extent
    ref = Gate_Param(i)
    if(abs(ref+Sig_Offset) > CUTOFF) then
        match_idx = 0
        do j = 1, this%N_Type
            if(abs(ref-this%Type_Param(j)) <= CUTOFF) then
                match_idx = j
                exit
            end if
        end do
        if(match_idx == 0) then
            this%N_Type = this%N_Type + 1
            this%Type_Param(this%N_Type) = ref
            match_idx = this%N_Type
        end if
        this%N_Per_Type(match_idx) = this%N_Per_Type(match_idx) + 1
        write(this%C_Num_bin(match_idx,this%N_Per_Type(match_idx)),"(B20.20)")i-1
        this%N_Total = this%N_Total + 1
    end if
end do
return

end subroutine csd_generator_GroupGates

subroutine csd_generator_ReduceGroups(this)

implicit none
class(csd_generator) :: this

integer :: N
character(len=20) :: workstr
integer :: i, j, k, row, p, limit
integer :: QSimilar, BitPos

N = this%N
do i = 1, this%N_Type
    if(this%N_Per_Type(i) == 1) cycle
    do row = 1, N-1
        j = 1
        do while(j < this%N_Per_Type(i))
            limit = this%N_Per_Type(i)
            do k = j+1, limit
                if(QSimilar(N,this%C_Num_Bin(i,j),this%C_Num_Bin(i,k),row) == 1) then
                    p = BitPos(N,row)
                    this%C_Num_Bin(i,j)(p:p) = '*'
                    if(k /= limit) this%C_Num_Bin(i,k) = this%C_Num_Bin(i,limit)
                    this%N_Per_Type(i) = limit - 1
                    this%N_Total = this%N_Total - 1
                    exit
                end if
            end do
            j = j+1
        end do
    end do
end do
return

end subroutine csd_generator_ReduceGroups

end module csd_real
