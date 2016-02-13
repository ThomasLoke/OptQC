module csd_cplx

use arrays_cplx
use common_module

implicit none

double precision, parameter :: CUTOFF = 0.0000000001d0
double precision, parameter :: PI = 3.1415926535897932384626433832795028841971693993751d0
double complex, parameter :: C_ZERO = cmplx(0.0d0,0.0d0)

type csd_solution_cplx
    integer :: N, M
    double complex, allocatable :: X(:,:)
    character(len=32), allocatable :: Circuit(:,:)          ! Note: By design, this is empty unless run_csdr is called
    integer :: csd_ct, csdr_ct
    logical :: neg                                          ! .true. if circuit represents -X, .false. otherwise
    logical :: toggle_csd                                   ! Runs csd/csdr procedure if .true., otherwise disallow any changes to variables (except via external subroutines)
contains
    procedure :: constructor => csd_solution_cplx_constructor
    procedure :: destructor => csd_solution_cplx_destructor
    procedure :: clean => csd_solution_cplx_clean
    procedure :: copy => csd_solution_cplx_copy
    procedure :: run_csd => csd_solution_cplx_run_csd
    procedure :: run_csdr => csd_solution_cplx_run_csdr
    procedure :: write_circuit => csd_solution_cplx_write_circuit
end type csd_solution_cplx

type csd_solution_set_cplx
    integer :: nset
    integer :: N, M
    type(csd_solution_cplx), allocatable :: arr(:)
    integer :: csd_ss_ct, csdr_ss_ct
    logical :: neg
contains
    procedure :: constructor => csd_solution_set_cplx_constructor
    procedure :: destructor => csd_solution_set_cplx_destructor
    procedure :: clean => csd_solution_set_cplx_clean
    procedure :: copy => csd_solution_set_cplx_copy
    procedure :: run_csd => csd_solution_set_cplx_run_csd
    procedure :: run_csdr => csd_solution_set_cplx_run_csdr
    procedure :: write_circuit => csd_solution_set_cplx_write_circuit
end type csd_solution_set_cplx

type csd_generator_cplx
    ! Includes standard CSD and the reduction procedures
    ! Global variables throughout all subroutines
    integer :: N, M, Mh                                     ! Reminder: M = 2**N, Mh = M/2
    integer, pointer :: index_level(:), index_pair(:,:,:)   ! Remains constant for matrices of the same dimensions
    double precision, pointer :: COEFF(:,:)                 ! Remains constant for matrices of the same dimensions
    double precision :: GATEPHASE
    double precision, allocatable :: GATEY(:,:), GATEZ(:,:)
    ! Note: Could just store a pointer to the csd_solution object instead - but done this way for less indirection (i.e. more clarity)
    double complex, pointer :: X(:,:)                           ! Intent: In
    character(len=32), pointer :: Circuit(:,:)              ! Intent: Out
    integer, pointer :: csd_ct, csdr_ct                     ! Intent: Out
    logical, pointer :: neg                                 ! Intent: Out
    ! Workspace arrays
    ! CYGC_CSD:
    type(l_arr_dp_4_cplx) :: Z0, Z1
    ! CYGC_BLKCSD:
    type(l_arr_dp_2_cplx) :: X_blk, X11, X12, X21, X22, U1, U2, V1T, V2T
    type(l_arr_dp_1) :: GATEY_blk
    type(l_arr_int_1) :: IWORK
    ! CYGC_CSDPHASE:
    double precision, allocatable :: PHASEZ(:,:), PARR(:), COEFF_WORK(:,:), IPIV(:)
    ! ReduceSolution:
    character(len=20), allocatable :: C_Num_Bin(:,:)
    double precision, allocatable :: Type_Param(:)
    integer, allocatable :: N_Per_Type(:)
    ! Scratch variables between CYGR_CSD and CYGR_BLKCSD
    integer :: size0, size1, num0, num1
    ! Scratch variables between ReduceSolution and (GroupGates,ReduceGroups)
    integer :: N_Type, N_Total
contains
    procedure :: constructor => csd_generator_cplx_constructor
    procedure :: destructor => csd_generator_cplx_destructor
    procedure :: assign_target => csd_generator_cplx_assign_target
    procedure :: nullify_target => csd_generator_cplx_nullify_target
    procedure :: run_csd => csd_generator_cplx_run_csd
    procedure :: run_blkcsd => csd_generator_cplx_run_blkcsd
    procedure :: run_csdphase => csd_generator_cplx_run_csdphase
    procedure :: GateCount => csd_generator_cplx_GateCount
    procedure :: ReduceSolution => csd_generator_cplx_ReduceSolution
    procedure :: GroupGates => csd_generator_cplx_GroupGates
    procedure :: ReduceGroups => csd_generator_cplx_ReduceGroups
end type csd_generator_cplx

contains

subroutine csd_solution_cplx_constructor(this,N,M)

implicit none
class(csd_solution_cplx) :: this
integer :: N, M

this%N = N
this%M = M
this%csd_ct = 0
this%csdr_ct = 0
allocate(this%X(M,M))
allocate(this%Circuit(N+1,M*M/2+M))
this%X = C_ZERO
this%Circuit = ""
this%neg = .false.
this%toggle_csd = .true.

end subroutine csd_solution_cplx_constructor

subroutine csd_solution_cplx_destructor(this)

implicit none
class(csd_solution_cplx) :: this

deallocate(this%X)
deallocate(this%Circuit)

end subroutine csd_solution_cplx_destructor

subroutine csd_solution_cplx_clean(this)

implicit none
class(csd_solution_cplx) :: this

this%Circuit = ""
this%csd_ct = 0
this%csdr_ct = 0
this%neg = .false.

end subroutine csd_solution_cplx_clean

subroutine csd_solution_cplx_copy(this,source)

implicit none
class(csd_solution_cplx) :: this
type(csd_solution_cplx) :: source

this%X = source%X
this%Circuit = source%Circuit
this%csd_ct = source%csd_ct
this%csdr_ct = source%csdr_ct
this%neg = source%neg
this%toggle_csd = source%toggle_csd

end subroutine csd_solution_cplx_copy

subroutine csd_solution_cplx_run_csd(this,generator)

implicit none
class(csd_solution_cplx) :: this
type(csd_generator_cplx) :: generator

if(this%toggle_csd == .true.) then
    call this%clean()
    call generator%assign_target(this)
    call generator%run_csd()                         ! Note: No circuit result
    call generator%nullify_target()
end if

end subroutine csd_solution_cplx_run_csd

subroutine csd_solution_cplx_run_csdr(this,generator)

implicit none
class(csd_solution_cplx) :: this
type(csd_generator_cplx) :: generator

if(this%toggle_csd == .true.) then
    call this%clean()
    call generator%assign_target(this)
    call generator%run_csd()
    call generator%ReduceSolution()                  ! Note: Circuit result produced here
    call generator%nullify_target()
end if

end subroutine csd_solution_cplx_run_csdr

subroutine csd_solution_cplx_write_circuit(this,wh)

implicit none
class(csd_solution_cplx) :: this
type(csd_write_handle) :: wh

if(this%csdr_ct == 0) then
    write(*,'(a)')"Attempted to output empty circuit to .tex file - request ignored."
    return
end if
call wh%clean()
call wh%preamble()
if(this%neg == .true.) call wh%insert_neg()
call wh%assign_target(this%Circuit,this%csdr_ct)
call wh%write_circuit()
call wh%nullify_target()
call wh%postamble()

end subroutine csd_solution_cplx_write_circuit

subroutine csd_solution_set_cplx_constructor(this,nset,N,M)

implicit none
class(csd_solution_set_cplx) :: this
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

end subroutine csd_solution_set_cplx_constructor

subroutine csd_solution_set_cplx_destructor(this)

implicit none
class(csd_solution_set_cplx) :: this

integer :: i

do i = 1, this%nset
    call this%arr(i)%destructor()
end do
deallocate(this%arr)

end subroutine csd_solution_set_cplx_destructor

subroutine csd_solution_set_cplx_clean(this)

implicit none
class(csd_solution_set_cplx) :: this

! Note: Cleaning of csd_solution objects is done individually in their call to run_csd or run_csdr
this%csd_ss_ct = 0
this%csdr_ss_ct = 0
this%neg = .false.

end subroutine csd_solution_set_cplx_clean

subroutine csd_solution_set_cplx_copy(this,source)

implicit none
class(csd_solution_set_cplx) :: this
type(csd_solution_set_cplx) :: source

integer :: i

do i = 1, this%nset
    call this%arr(i)%copy(source%arr(i))
end do
this%csd_ss_ct = source%csd_ss_ct
this%csdr_ss_ct = source%csdr_ss_ct
this%neg = source%neg

end subroutine csd_solution_set_cplx_copy

subroutine csd_solution_set_cplx_run_csd(this,generator)

implicit none
class(csd_solution_set_cplx) :: this
type(csd_generator_cplx) :: generator

integer :: i

call this%clean()
do i = 1, this%nset
    call this%arr(i)%run_csd(generator)
    this%csd_ss_ct = this%csd_ss_ct + this%arr(i)%csd_ct
    if(this%arr(i)%neg == .true.) this%neg = .not.this%neg
end do

end subroutine csd_solution_set_cplx_run_csd

subroutine csd_solution_set_cplx_run_csdr(this,generator)

implicit none
class(csd_solution_set_cplx) :: this
type(csd_generator_cplx) :: generator

integer :: i

call this%clean()
do i = 1, this%nset
    call this%arr(i)%run_csdr(generator)
    this%csd_ss_ct = this%csd_ss_ct + this%arr(i)%csd_ct
    this%csdr_ss_ct = this%csdr_ss_ct + this%arr(i)%csdr_ct
    if(this%arr(i)%neg == .true.) this%neg = .not.this%neg
end do

end subroutine csd_solution_set_cplx_run_csdr

subroutine csd_solution_set_cplx_write_circuit(this,wh)

implicit none
class(csd_solution_set_cplx) :: this
type(csd_write_handle) :: wh

integer :: i

call wh%clean()
call wh%preamble()
if(this%neg == .true.) call wh%insert_neg()
! Note: Order reversal in circuit, since U = U_1 U_2 .. U_N means that U_N is applied first
do i = this%nset, 1, -1
    if(this%arr(i)%csdr_ct > 0) then
        call wh%assign_target(this%arr(i)%Circuit,this%arr(i)%csdr_ct)
        call wh%write_circuit()
        if(i /= 1) then
            call wh%insert_seperator()
        end if
        call wh%nullify_target()
    end if
end do
call wh%postamble()

end subroutine csd_solution_set_cplx_write_circuit

subroutine csd_generator_cplx_constructor(this,N,M,index_level,index_pair,COEFF)

implicit none
class(csd_generator_cplx) :: this
integer :: N, M
integer, target :: index_level(M-1), index_pair(M/2,2,N)
double precision, target :: COEFF(M,M)
integer :: i, Mh, size0, size1, num0, num1

Mh = M / 2
this%N = N
this%M = M
this%Mh = Mh
this%index_level => index_level
this%index_pair => index_pair
this%COEFF => COEFF
allocate(this%GATEY(M/2,M-1))
allocate(this%GATEZ(M/2,M+N-1))
call this%nullify_target()
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
allocate(this%PHASEZ(M,M+N))
allocate(this%PARR(M))
allocate(this%COEFF_WORK(M,M))
allocate(this%IPIV(M))
allocate(this%C_Num_Bin(Mh,Mh))
allocate(this%Type_Param(Mh))
allocate(this%N_Per_Type(Mh))
this%size0 = 0
this%size1 = 0
this%num0 = 0
this%num1 = 0

end subroutine csd_generator_cplx_constructor

subroutine csd_generator_cplx_destructor(this)

implicit none
class(csd_generator_cplx) :: this

nullify(this%index_level)
nullify(this%index_pair)
nullify(this%COEFF)
deallocate(this%GATEY)
deallocate(this%GATEZ)
call this%nullify_target()
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
deallocate(this%PHASEZ)
deallocate(this%PARR)
deallocate(this%COEFF_WORK)
deallocate(this%IPIV)
deallocate(this%C_Num_Bin)
deallocate(this%Type_Param)
deallocate(this%N_Per_Type)

end subroutine csd_generator_cplx_destructor

subroutine csd_generator_cplx_assign_target(this,csd_ss)

implicit none
class(csd_generator_cplx) :: this
type(csd_solution_cplx), target :: csd_ss

! Bind pointers to the csd_solution elements
this%X => csd_ss%X
this%Circuit => csd_ss%Circuit
this%Circuit = ''
this%csd_ct => csd_ss%csd_ct
this%csd_ct = 0
this%csdr_ct => csd_ss%csdr_ct
this%csdr_ct = 0
this%neg => csd_ss%neg
this%neg = .false.

end subroutine csd_generator_cplx_assign_target

subroutine csd_generator_cplx_nullify_target(this)

implicit none
class(csd_generator_cplx) :: this

nullify(this%X)
nullify(this%Circuit)
nullify(this%csd_ct)
nullify(this%csdr_ct)
nullify(this%neg)

end subroutine csd_generator_cplx_nullify_target

subroutine csd_generator_cplx_run_csd(this)

implicit none
class(csd_generator_cplx) :: this

integer :: N, i, j

N = this%N
! Clean up intermediate variables first
this%GATEPHASE = 0.0d0
this%GATEY = 0.0d0
this%GATEZ = 0.0d0

!	Cosine Sine Decomposition (CSD) of X recursively until the Nth level
!	The 1st level
this%size0 = this%M
this%size1 = this%Mh
this%num0 = 1
this%num1 = 2
this%Z1%l(1)%arr = C_ZERO
call this%run_blkcsd(1,1)
this%Z0%l(1)%arr = this%Z1%l(1)%arr

!	The ist level (i = 2,...,N)
do i = 2, N
    this%size1 = 2**(N-i)
    this%size0 = this%size1 * 2
    this%num0 = 2**(i-1)
    this%num1 = this%num0 * 2
    this%Z1%l(i)%arr = C_ZERO
    do j = 1, this%num0
        call this%run_blkcsd(i,j)
    end do
    this%Z0%l(i)%arr = this%Z1%l(i)%arr
end do

!	Step 3: CYGC_CSD processes the CSD results
!	In this step, we obtain all the Rz Gates and the final Phase Gate for X: GATEZ, GATEPHASEcall this%run_cutgate()
this%COEFF_WORK = this%COEFF
call this%run_csdphase()

!   Step 4 (added): Obtain count of gates
call this%GateCount()

end subroutine csd_generator_cplx_run_csd

subroutine csd_generator_cplx_run_blkcsd(this,i,j)

implicit none
class(csd_generator_cplx), target :: this
integer :: i, j

integer :: index_Y, size0, size1, num0, num1
double complex, pointer :: X(:,:,:), U(:,:,:), V(:,:,:), X_blk(:,:), X11(:,:), X12(:,:), X21(:,:), X22(:,:), U1(:,:), U2(:,:), V1T(:,:), V2T(:,:)
double precision, pointer :: GATEY(:), GATEY_blk(:)
integer, pointer :: IWORK(:)
double complex, allocatable :: WORK(:)
double precision, allocatable :: RWORK(:)
integer :: LWORK, LRWORK, INFO
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

X_blk = C_ZERO
X11 = C_ZERO
X12 = C_ZERO
X21 = C_ZERO
X22 = C_ZERO
U1 = C_ZERO
U2 = C_ZERO
V1T = C_ZERO
V2T = C_ZERO
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
				X_blk(k1,k2) = C_ZERO
			end if
		end do
	end do
	X11 = X_blk(1:size1,1:size1)
	X12 = X_blk(1:size1,size1+1:size0)
	X21 = X_blk(size1+1:size0,1:size1)
	X22 = X_blk(size1+1:size0,size1+1:size0)
	LWORK = -1
	LRWORK = -1
	allocate(WORK(1))
	WORK = C_ZERO
	allocate(RWORK(1))
	RWORK = 0.0d0
    !	Set LWORK = -1, then a workspace query is assumed; the subroutine only calculates the optimal size of the WORK array, returns this value as the first entry of the work array.
    !	Set LRWORK = -1, then a workspace query is assumed; the subroutine only calculates the optimal size of the RWORK array, returns this value as the first entry of the work array.
	call ZUNCSD('Y', 'Y', 'Y', 'Y', 'N', 'O', size0, size1, size1, X11, &
size1, X12, size1, X21, size1, X22, size1, GATEY_blk, U1, size1, U2,&
size1, V1T, size1, V2T, size1, WORK, LWORK, RWORK, LRWORK, IWORK, &
INFO)
	!	Now set LWORK = WORK(1)
    !	Now set LRWORK = RWORK(1)
	LWORK = int(WORK(1))
	LRWORK = int(RWORK(1))
	deallocate(WORK)
	deallocate(RWORK)
	allocate(WORK(LWORK))
	WORK = C_ZERO
	allocate(RWORK(LRWORK))
	RWORK = 0.0d0
	call ZUNCSD('Y', 'Y', 'Y', 'Y', 'N', 'O', size0, size1, size1, X11, &
size1, X12, size1, X21, size1, X22, size1, GATEY_blk, U1, size1, U2, &
size1, V1T, size1, V2T, size1, WORK, LWORK, RWORK, LRWORK, IWORK, &
INFO)
	deallocate(WORK)
	deallocate(RWORK)
	U(:,:,2*k-1) = U1
	U(:,:,2*k) = U2
	V(:,:,2*k-1) = V1T
	V(:,:,2*k) = V2T
	GATEY(size1*(k-1)+1:size1*k) = GATEY_blk
end do

GATEY = GATEY/PI
return

end subroutine csd_generator_cplx_run_blkcsd

subroutine csd_generator_cplx_run_csdphase(this)

implicit none
class(csd_generator_cplx), target :: this

integer :: N, M
double complex, pointer :: Z(:,:)
integer :: i, j, idx, INFO

N = this%N
M = this%M

! Initialize pointers to elements of the csd_generator object
Z => this%Z0%l(N)%arr(1,1,:,:)

!	Z(j) = diag(exp(i*Pi*PHASEZ(j,1:M))) (j = 1,2,...,2**N)
!	In this SUBROUTINE, to make it short, the input Z(j) = exp(i*Pi*PHASEZ(j,1:M)), namely, Z(j) is a M-length array, and Z is a 2D table with the dimension
!	(M,M+N)
!	Here we transform Z/ into PHASEZ
this%PHASEZ = atan2(aimag(Z),real(Z))/PI

!	Decide the Rz Gates "GATEZ" (The first (M-1) Rz Gates)
do i = 1, M-1
	idx = this%index_level(i)
	this%PARR = 0.0d0
	do j = 1, M/2
		this%PARR(this%index_pair(j,:,idx)) = -0.5d0*sum(this%PHASEZ(this%index_pair(j,:,idx),i))
	end do
	this%PHASEZ(:,i) = this%PHASEZ(:,i) + this%PARR
	this%GATEZ(:,i) = this%PHASEZ(this%index_pair(:,1,idx),i)
	this%PHASEZ(:,i+1) = this%PHASEZ(:,i+1) - this%PARR
end do

!	Decide the Rz Gates "GATEZ" (The last N Rz Gates) and the Final Phase Gate "GATEPHASE"
this%IPIV = 0.0d0
call DGESV(M,1,this%COEFF_WORK,M,this%IPIV,this%PHASEZ(:,M),M,INFO)
do i = 1, N
	this%GATEZ(1:2**(N-i),M+i-1) = this%PHASEZ(2**N-2**(N+1-i)+1:2**N-2**(N-i),M)
end do
this%GATEPHASE = this%PHASEZ(M,M)

!	Set "GATEZ" and "GATEPHASE" within the range of (-1,1)
this%GATEZ = this%GATEZ-2.0d0*nint(this%GATEZ*0.5d0)
this%GATEPHASE = this%GATEPHASE-2.0d0*nint(this%GATEPHASE*0.5d0)
return

end subroutine csd_generator_cplx_run_csdphase

subroutine csd_generator_cplx_GateCount(this)

implicit none
class(csd_generator_cplx) :: this

integer :: i, j, idx, res

if(abs(this%GATEPHASE) > CUTOFF) then
    res = 1
else
    res = 0
end if
do i = 1, this%N
    idx = this%M+this%N-i
    do j = 1, 2**(i-1)
		if (abs(this%GATEZ(j,idx)) > CUTOFF) then
			res = res + 1
		end if
	end do
end do
do i = this%M-1, 1, -1
    do j = 1, this%M/2
		if (abs(this%GATEY(j,i)) > CUTOFF) then
			res = res + 1
		end if
		if (abs(this%GATEZ(j,i)) > CUTOFF) then
			res = res + 1
		end if
	end do
end do
this%csd_ct = res
return

end subroutine csd_generator_cplx_GateCount

subroutine csd_generator_cplx_ReduceSolution(this)

implicit none
class(csd_generator_cplx) :: this
integer :: N, M, ct
integer :: i, j, k, p, ref, row, lb, rb
character(len=20) :: workstr
integer :: IsEmpty

! Includes update of the count of the reduced number of gates at the end
N = this%N
M = this%M
ct = 0
this%Circuit = ""
! Form the circuit for the reduced solution.
! GATEPHASE
if(abs(this%GATEPHASE) > CUTOFF) then
    this%Circuit(1,1) = '& \gate{\Phi}'
    this%Circuit(2:N,1) = '& \qw'
    write(this%Circuit(N+1,1),302) this%GATEPHASE
    ct = 1
end if
! GATEZ
do i = 1,N
    lb = 22-i
    call this%GroupGates(2**(i-1),this%GATEZ(:,M+N-i),0.0d0)
    call this%ReduceGroups()
    do j = 1, this%N_Type
        do k = 1, this%N_Per_Type(j)
            ct = ct+1
            if(i > 1) then
                workstr = this%C_Num_Bin(j,k)
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
            this%Circuit(i,ct) = '& \gate{R_z}'
            if(i < N) then
                do row = i+1,N
                    this%Circuit(row,ct) = '& \qw'
                end do
            end if
            write(this%Circuit(N+1,ct),302)this%Type_Param(j)
        end do
    end do
end do
! GATEY,GATEZ
lb = 22-N
rb = 20
do i = M-1, 1, -1
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
                    ! Qubits before target
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
                    ! Qubits after target
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
    call this%GroupGates(this%Mh,this%GATEZ(:,i),0.0d0)
    call this%ReduceGroups()
    do j = 1, this%N_Type
        do k = 1, this%N_Per_Type(j)
            ct = ct+1
            workstr = this%C_Num_Bin(j,k)
            this%Circuit(ref,ct) = '& \gate{R_z}'
            do row = 1,N-1
                p = 21-N+row
				if(row < ref) then
                    ! Qubits before target
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
                    ! Qubits after target
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
302 format('& ',F7.4)
this%csdr_ct = ct

end subroutine csd_generator_cplx_ReduceSolution

subroutine csd_generator_cplx_GroupGates(this,extent,Gate_Param,Sig_Offset)

implicit none
class(csd_generator_cplx) :: this
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

end subroutine csd_generator_cplx_GroupGates

subroutine csd_generator_cplx_ReduceGroups(this)

implicit none
class(csd_generator_cplx) :: this

integer :: N
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

end subroutine csd_generator_cplx_ReduceGroups

end module csd_cplx
