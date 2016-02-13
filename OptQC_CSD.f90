module csd_tools

use arrays_real
use arrays_cplx
use common_module
use rng

implicit none

double precision, parameter :: CUTOFF = 0.0000000001d0
double precision, parameter :: PI = 3.1415926535897932384626433832795028841971693993751d0
double complex, parameter :: C_ZERO = cmplx(0.0d0,0.0d0)

type csd_solution
    integer :: N, M
    integer :: obj_type                                     ! 0 if real, 1 if complex
    double precision, allocatable :: X(:,:)
    double complex, allocatable :: Xc(:,:)
    character(len=32), allocatable :: Circuit(:,:)          ! Note: By design, this is empty unless run_csdr is called
    integer :: csd_ct, csdr_ct
    logical :: neg                                          ! .true. if circuit represents -X, .false. otherwise - NOT USED IF CPLX
    logical :: toggle_csd                                   ! Runs csd/csdr procedure if .true., otherwise disallow any changes to variables (except via external subroutines)
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
    logical :: neg
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
    integer :: obj_type                                     ! 0 if used only for real matrices, 1 if otherwise
    integer, pointer :: index_level(:), index_pair(:,:,:)   ! Remains constant for matrices of the same dimensions
    double precision, pointer :: COEFF(:,:)                 ! Remains constant for matrices of the same dimensions
    ! Intermediate result variables
    !!!!!!!!!!!!!!!
    !!! COMMON  !!!
    double precision, allocatable :: GATEY(:,:)
    !!! IF REAL !!!
    double precision, allocatable :: GATEPI(:,:)
    !!! IF CPLX !!!
    double precision :: GATEPHASE
    double precision, allocatable :: GATEZ(:,:)
    !!!!!!!!!!!!!!!
    ! Note: Could just store a pointer to the csd_solution object instead - but done this way for less indirection (i.e. more clarity)
    integer :: targ_type                                    ! Copied
    double precision, pointer :: X(:,:)                     ! Intent: In
    double complex, pointer :: Xc(:,:)                      ! Intent: In
    character(len=32), pointer :: Circuit(:,:)              ! Intent: Out
    integer, pointer :: csd_ct, csdr_ct                     ! Intent: Out
    logical, pointer :: neg                                 ! Intent: Out
    ! Workspace arrays
    ! CSD:
    !!!!!!!!!!!!!!!
    !!! IF REAL !!!
    type(l_arr_dp_4) :: Z0, Z1
    !!! IF CPLX !!!
    type(l_arr_dp_4_cplx) :: Z0c, Z1c
    !!!!!!!!!!!!!!!
    ! BLKCSD:
    !!!!!!!!!!!!!!!
    !!! COMMON  !!!
    type(l_arr_dp_1) :: GATEY_blk
    type(l_arr_int_1) :: IWORK
    !!! IF REAL !!!
    type(l_arr_dp_2) :: X_blk, X11, X12, X21, X22, U1, U2, V1T, V2T
    !!! IF CPLX !!!
    type(l_arr_dp_2_cplx) :: X_blkc, X11c, X12c, X21c, X22c, U1c, U2c, V1Tc, V2Tc
    !!!!!!!!!!!!!!!
    ! CUTGATE (REAL ONLY):
    integer, allocatable :: Z_array(:,:), PI_array(:,:), GATEY_sign(:)
    ! CSDPHASE (CPLX ONLY):
    double complex, allocatable :: Z_temp(:,:)
    double precision, allocatable :: PHASEZ(:,:), PARR(:), COEFF_WORK(:,:), IPIV(:)
    ! ReduceSolution:
    type(binstr), allocatable :: C_Num_Bin(:,:)
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
    procedure :: run_cplx_blkcsd => csd_generator_run_cplx_blkcsd
    procedure :: run_cutgate => csd_generator_run_cutgate
    procedure :: run_csdphase => csd_generator_run_csdphase
    procedure :: GateCount => csd_generator_GateCount
    procedure :: ReduceSolution => csd_generator_ReduceSolution
    procedure :: GroupGates => csd_generator_GroupGates
    procedure :: ReduceGroups => csd_generator_ReduceGroups
end type csd_generator

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

contains

subroutine csd_solution_constructor(this,N,M,obj_type)

implicit none
class(csd_solution) :: this
integer :: N, M, obj_type

this%N = N
this%M = M
this%obj_type = obj_type
this%csd_ct = 0
this%csdr_ct = 0
! Note: Edited circuit allocation size to be equal for both types
! Warning: Hardcoded - assumed in OptQC_Perm.f90 as well
if(obj_type == 0) then
    allocate(this%X(M,M))
    this%X = 0.0d0
    !allocate(this%Circuit(N+1,M*M/2+M))
else
    allocate(this%Xc(M,M))
    this%Xc = C_ZERO
    !allocate(this%Circuit(N+1,M*M+1))
end if
allocate(this%Circuit(N+1,M*M+1))
this%Circuit = ""
this%neg = .false.
this%toggle_csd = .true.

end subroutine csd_solution_constructor

subroutine csd_solution_destructor(this)

implicit none
class(csd_solution) :: this

if(this%obj_type == 0) then
    deallocate(this%X)
else
    deallocate(this%Xc)
end if
deallocate(this%Circuit)

end subroutine csd_solution_destructor

subroutine csd_solution_clean(this)

implicit none
class(csd_solution) :: this

this%Circuit = ""
this%csd_ct = 0
this%csdr_ct = 0
this%neg = .false.

end subroutine csd_solution_clean

subroutine csd_solution_copy(this,source)

implicit none
class(csd_solution) :: this, source

if(this%obj_type /= source%obj_type) then
    write(*,'(a)')"Attempted to copy between csd_solution objects of different type. Terminating program."
    call exit(1)
end if
if(this%obj_type == 0) then
    this%X = source%X
else
    this%Xc = source%Xc
end if
this%Circuit = source%Circuit
this%csd_ct = source%csd_ct
this%csdr_ct = source%csdr_ct
this%neg = source%neg
this%toggle_csd = source%toggle_csd

end subroutine csd_solution_copy

subroutine csd_solution_run_csd(this,generator)

implicit none
class(csd_solution) :: this
class(csd_generator) :: generator

if(this%toggle_csd == .true.) then
    call this%clean()
    call generator%assign_target(this)
    call generator%run_csd()                         ! Note: No circuit result
    call generator%nullify_target()
end if

end subroutine csd_solution_run_csd

subroutine csd_solution_run_csdr(this,generator)

implicit none
class(csd_solution) :: this
class(csd_generator) :: generator

if(this%toggle_csd == .true.) then
    call this%clean()
    call generator%assign_target(this)
    call generator%run_csd()
    call generator%ReduceSolution()                  ! Note: Circuit result produced here
    call generator%nullify_target()
end if

end subroutine csd_solution_run_csdr

subroutine csd_solution_write_circuit(this,wh)

implicit none
class(csd_solution) :: this
class(csd_write_handle) :: wh

if(this%csdr_ct == 0) then
    write(*,'(a)')"Attempted to output empty circuit to .tex file - request ignored."
    return
end if
call wh%clean()
call wh%preamble()
if(this%neg == .true.) call wh%insert_neg()
call wh%assign_target(this)
call wh%write_circuit()
call wh%nullify_target()
call wh%postamble()

end subroutine csd_solution_write_circuit

subroutine csd_solution_set_constructor(this,N,M,nset,type_spec)

implicit none
class(csd_solution_set) :: this
integer :: N, M, nset
integer :: type_spec(nset)

integer :: i

this%N = N
this%M = M
this%nset = nset
allocate(this%arr(nset))
do i = 1, nset
    call this%arr(i)%constructor(N,M,type_spec(i))
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
this%neg = .false.

end subroutine csd_solution_set_clean

subroutine csd_solution_set_copy(this,source)

implicit none
class(csd_solution_set) :: this, source
integer :: i

do i = 1, this%nset
    call this%arr(i)%copy(source%arr(i))
end do
this%csd_ss_ct = source%csd_ss_ct
this%csdr_ss_ct = source%csdr_ss_ct
this%neg = source%neg

end subroutine csd_solution_set_copy

subroutine csd_solution_set_run_csd(this,generator)

implicit none
class(csd_solution_set) :: this
class(csd_generator) :: generator

integer :: i

call this%clean()
do i = 1, this%nset
    call this%arr(i)%run_csd(generator)
    this%csd_ss_ct = this%csd_ss_ct + this%arr(i)%csd_ct
    if(this%arr(i)%neg == .true.) this%neg = .not.this%neg
end do

end subroutine csd_solution_set_run_csd

subroutine csd_solution_set_run_csdr(this,generator)

implicit none
class(csd_solution_set) :: this
class(csd_generator) :: generator

integer :: i

call this%clean()
do i = 1, this%nset
    call this%arr(i)%run_csdr(generator)
    this%csd_ss_ct = this%csd_ss_ct + this%arr(i)%csd_ct
    this%csdr_ss_ct = this%csdr_ss_ct + this%arr(i)%csdr_ct
    if(this%arr(i)%neg == .true.) this%neg = .not.this%neg
end do

end subroutine csd_solution_set_run_csdr

subroutine csd_solution_set_write_circuit(this,wh)

implicit none
class(csd_solution_set) :: this
class(csd_write_handle) :: wh

integer :: i, sepidx, sepct

call wh%clean()
call wh%preamble()
sepidx = 0
sepct = -1
if(this%neg == .true.) call wh%insert_neg()
! Precount number of seperators required
do i = this%nset, 1, -1
    if(this%arr(i)%csdr_ct > 0) then
        sepct = sepct + 1
    end if
end do
! Note: Order reversal in circuit, since U = U_1 U_2 .. U_N means that U_N is applied first
do i = this%nset, 1, -1
    if(this%arr(i)%csdr_ct > 0) then
        call wh%assign_target(this%arr(i))
        call wh%write_circuit()
        if(sepidx < sepct) then
            call wh%insert_seperator()
            sepidx = sepidx + 1
        end if
        call wh%nullify_target()
    end if
end do
call wh%postamble()

end subroutine csd_solution_set_write_circuit

subroutine csd_generator_constructor(this,N,M,obj_type,index_level,index_pair,COEFF)

implicit none
class(csd_generator) :: this
integer :: N, M, obj_type
integer, target :: index_level(M-1), index_pair(M/2,2,N)
double precision, target :: COEFF(M,M)
integer :: i, j, Mh, size0, size1, num0, num1

Mh = M / 2
this%N = N
this%M = M
this%Mh = Mh
this%obj_type = obj_type
this%index_level => index_level
this%index_pair => index_pair
this%COEFF => COEFF
allocate(this%GATEY(Mh,M-1))
allocate(this%GATEPI(Mh,N))
if(obj_type /= 0) then
    allocate(this%GATEZ(Mh,M+N-1))
end if
call this%nullify_target()
call this%Z0%constructor(N)
call this%Z1%constructor(N)
if(obj_type /= 0) then
    call this%Z0c%constructor(N)
    call this%Z1c%constructor(N)
end if
call this%GATEY_blk%constructor(N)
call this%IWORK%constructor(N)
call this%X_blk%constructor(N)
call this%X11%constructor(N)
call this%X12%constructor(N)
call this%X21%constructor(N)
call this%X22%constructor(N)
call this%U1%constructor(N)
call this%U2%constructor(N)
call this%V1T%constructor(N)
call this%V2T%constructor(N)
if(obj_type /= 0) then
    call this%X_blkc%constructor(N)
    call this%X11c%constructor(N)
    call this%X12c%constructor(N)
    call this%X21c%constructor(N)
    call this%X22c%constructor(N)
    call this%U1c%constructor(N)
    call this%U2c%constructor(N)
    call this%V1Tc%constructor(N)
    call this%V2Tc%constructor(N)
end if
do i = 1, N
    size1 = 2**(N-i)
    size0 = size1 * 2
    num0 = 2**(i-1)
    num1 = num0 * 2
    call this%GATEY_blk%l(i)%constructor(size1)
    call this%IWORK%l(i)%constructor(size1)
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
    if(obj_type /= 0) then
        call this%Z0c%l(i)%constructor(size1,size1,num1,num1)
        call this%Z1c%l(i)%constructor(size1,size1,num1,num1)
        call this%X_blkc%l(i)%constructor(size0,size0)
        call this%X11c%l(i)%constructor(size1,size1)
        call this%X12c%l(i)%constructor(size1,size1)
        call this%X21c%l(i)%constructor(size1,size1)
        call this%X22c%l(i)%constructor(size1,size1)
        call this%U1c%l(i)%constructor(size1,size1)
        call this%U2c%l(i)%constructor(size1,size1)
        call this%V1Tc%l(i)%constructor(size1,size1)
        call this%V2Tc%l(i)%constructor(size1,size1)
    end if
end do
allocate(this%Z_array(M,M))
allocate(this%PI_array(M,N))
allocate(this%GATEY_sign(Mh))
if(obj_type /= 0) then
    allocate(this%Z_temp(M,M+N))
    allocate(this%PHASEZ(M,M+N))
    allocate(this%PARR(M))
    allocate(this%COEFF_WORK(M,M))
    allocate(this%IPIV(M))
end if
allocate(this%C_Num_Bin(Mh,Mh))
do i = 1, Mh
    do j = 1, Mh
        call this%C_Num_Bin(i,j)%constructor(N-1)
    end do
end do
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
integer :: i, j

nullify(this%index_level)
nullify(this%index_pair)
nullify(this%COEFF)
deallocate(this%GATEY)
deallocate(this%GATEPI)
if(this%obj_type /= 0) then
    deallocate(this%GATEZ)
end if
call this%nullify_target()
call this%Z0%destructor()
call this%Z1%destructor()
if(this%obj_type /= 0) then
    call this%Z0c%destructor()
    call this%Z1c%destructor()
end if
call this%GATEY_blk%destructor()
call this%IWORK%destructor()
call this%X_blk%destructor()
call this%X11%destructor()
call this%X12%destructor()
call this%X21%destructor()
call this%X22%destructor()
call this%U1%destructor()
call this%U2%destructor()
call this%V1T%destructor()
call this%V2T%destructor()
if(this%obj_type /= 0) then
    call this%X_blkc%destructor()
    call this%X11c%destructor()
    call this%X12c%destructor()
    call this%X21c%destructor()
    call this%X22c%destructor()
    call this%U1c%destructor()
    call this%U2c%destructor()
    call this%V1Tc%destructor()
    call this%V2Tc%destructor()
end if
deallocate(this%Z_array)
deallocate(this%PI_array)
deallocate(this%GATEY_sign)
if(this%obj_type /= 0) then
    deallocate(this%Z_temp)
    deallocate(this%PHASEZ)
    deallocate(this%PARR)
    deallocate(this%COEFF_WORK)
    deallocate(this%IPIV)
end if
do i = 1, this%Mh
    do j = 1, this%Mh
        call this%C_Num_Bin(i,j)%destructor()
    end do
end do
deallocate(this%C_Num_Bin)
deallocate(this%Type_Param)
deallocate(this%N_Per_Type)

end subroutine csd_generator_destructor

subroutine csd_generator_assign_target(this,csd_ss)

implicit none
class(csd_generator) :: this
class(csd_solution), target :: csd_ss

! Bind pointers to the csd_solution elements
this%targ_type = csd_ss%obj_type
if(this%targ_type == 1 .and. this%obj_type == 0) then
    write(*,*)"Error! Illegal usage of a real csd_generator object to decompose a complex matrix. Terminating program."
    call exit(1)
end if
if(this%targ_type == 0) then
    this%X => csd_ss%X
else
    this%Xc => csd_ss%Xc
end if
this%Circuit => csd_ss%Circuit
this%Circuit = ''
this%csd_ct => csd_ss%csd_ct
this%csd_ct = 0
this%csdr_ct => csd_ss%csdr_ct
this%csdr_ct = 0
this%neg => csd_ss%neg
this%neg = .false.

end subroutine csd_generator_assign_target

subroutine csd_generator_nullify_target(this)

implicit none
class(csd_generator) :: this

if(this%targ_type == 0) then
    nullify(this%X)
else
    nullify(this%Xc)
end if
nullify(this%Circuit)
nullify(this%csd_ct)
nullify(this%csdr_ct)
nullify(this%neg)
this%targ_type = 0

end subroutine csd_generator_nullify_target

subroutine csd_generator_run_csd(this)

implicit none
class(csd_generator) :: this

integer :: N, i, j

N = this%N
! Clean up intermediate variables first
this%GATEY = 0.0d0

if(this%targ_type == 0) then
    this%GATEPI = 0.0d0
    !	Cosine Sine Decomposition (CSD) of X recursively until the Nth level
    !	The 1st level
    this%size0 = this%M
    this%size1 = this%Mh
    this%num0 = 1
    this%num1 = 2
    this%Z1%l(1)%arr = 0.0d0
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
else
    this%GATEPHASE = 0.0d0
    this%GATEZ = 0.0d0
    !	Cosine Sine Decomposition (CSD) of X recursively until the Nth level
    !	The 1st level
    this%size0 = this%M
    this%size1 = this%Mh
    this%num0 = 1
    this%num1 = 2
    this%Z1c%l(1)%arr = C_ZERO
    call this%run_cplx_blkcsd(1,1)
    this%Z0c%l(1)%arr = this%Z1c%l(1)%arr
    !	The ist level (i = 2,...,N)
    do i = 2, N
        this%size1 = 2**(N-i)
        this%size0 = this%size1 * 2
        this%num0 = 2**(i-1)
        this%num1 = this%num0 * 2
        this%Z1c%l(i)%arr = C_ZERO
        do j = 1, this%num0
            call this%run_cplx_blkcsd(i,j)
        end do
        this%Z0c%l(i)%arr = this%Z1c%l(i)%arr
    end do
    !	Step 3: CYGC_CSD processes the CSD results
    !	In this step, we obtain all the Rz Gates and the final Phase Gate for X: GATEZ, GATEPHASE
    this%COEFF_WORK = this%COEFF
    call this%run_csdphase()
end if

!   Last step (added): Obtain count of gates
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
				X_blk(k1,k2) = 0.0d0
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
	LWORK = int(WORK(1))
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

subroutine csd_generator_run_cplx_blkcsd(this,i,j)

implicit none
class(csd_generator), target :: this
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
    X => this%Z0c%l(i-1)%arr(:,:,:,j)
end if
U => this%Z1c%l(i)%arr(:,:,:,2*j-1)
V => this%Z1c%l(i)%arr(:,:,:,2*j)
GATEY => this%GATEY(:,index_Y)
X_blk => this%X_blkc%l(i)%arr
X11 => this%X11c%l(i)%arr
X12 => this%X12c%l(i)%arr
X21 => this%X21c%l(i)%arr
X22 => this%X22c%l(i)%arr
U1 => this%U1c%l(i)%arr
U2 => this%U2c%l(i)%arr
V1T => this%V1Tc%l(i)%arr
V2T => this%V2Tc%l(i)%arr
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
        X_blk = this%Xc(:,:)
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

end subroutine csd_generator_run_cplx_blkcsd

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
	this%neg = .true.
else
    this%neg = .false.
end if

do i = 1,N
	leng = 2**(N+1-i)
	do j = 1,2**(i-1)
		this%GATEPI(j,i) = this%PI_array(leng*j,i)
	end do
end do
return

end subroutine csd_generator_run_cutgate

subroutine csd_generator_run_csdphase(this)

implicit none
class(csd_generator) :: this

integer :: N, M
integer :: i, j, idx, INFO

N = this%N
M = this%M
INFO = 0

! Initialize Z_temp to prescribed block
this%Z_temp = C_ZERO
this%Z_temp(:,1:M) = this%Z0c%l(N)%arr(1,1,:,:)

!	Z(j) = diag(exp(i*Pi*PHASEZ(j,1:M))) (j = 1,2,...,2**N)
!	In this SUBROUTINE, to make it short, the input Z(j) = exp(i*Pi*PHASEZ(j,1:M)), namely, Z(j) is a M-length array, and Z is a 2D table with the dimension
!	(M,M+N)
!	Here we transform Z/ into PHASEZ
this%PHASEZ = atan2(aimag(this%Z_temp),real(this%Z_temp))/PI

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
	this%GATEZ(1:2**(N-i),M+i-1) = this%PHASEZ(M-2**(N+1-i)+1:2**N-2**(N-i),M)
end do
this%GATEPHASE = this%PHASEZ(M,M)

!	Set "GATEZ" and "GATEPHASE" within the range of (-1,1)
this%GATEZ = this%GATEZ-2.0d0*nint(this%GATEZ*0.5d0)
this%GATEPHASE = this%GATEPHASE-2.0d0*nint(this%GATEPHASE*0.5d0)
return

end subroutine csd_generator_run_csdphase

subroutine csd_generator_GateCount(this)

implicit none
class(csd_generator) :: this

integer :: i, j, idx, res

res = 0
if(this%targ_type == 0) then
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
else
    if(abs(this%GATEPHASE) > CUTOFF) then
        res = 1
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
end if
this%csd_ct = res
return

end subroutine csd_generator_GateCount

subroutine csd_generator_ReduceSolution(this)

implicit none
class(csd_generator) :: this
integer :: N, M, ct
integer :: i, j, k, p, ref, row, lb, rb
type(binstr) :: workstr

! Includes update of the count of the reduced number of gates at the end
N = this%N
M = this%M
ct = 0
this%Circuit = ""
call workstr%constructor(N-1)
if(this%targ_type == 0) then
    ! Form the circuit for the reduced solution.
    ! GATEPI
    ct = 0
    do i = 1, N
        lb = N+1-i ! Equivalent to workstr%len+2-i
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
                call workstr%copy(this%C_Num_Bin(1,k))
                do row = 1, i-1
                    p = lb-1+row ! Equivalent to N-i+row
                    if(workstr%str(p)=='0') then
                        this%Circuit(row,ct) = '& \ctrlo{1}'
                    else if(workstr%str(p)=='1') then
                        this%Circuit(row,ct) = '& \ctrl{1}'
                    else
                        ! Case where workstr%str(p) = '*'
                        if(row > 1 .and. workstr%isempty(lb,p-1)==0) then
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
    lb = 1
    rb = N-1
    do i = M-1, 1, -1
        ref = this%index_level(i)
        call this%GroupGates(this%Mh,this%GATEY(:,i),0.0d0)
        call this%ReduceGroups()
        do j = 1, this%N_Type
            do k = 1, this%N_Per_Type(j)
                ct = ct+1
                call workstr%copy(this%C_Num_Bin(j,k))
                this%Circuit(ref,ct) = '& \gate{R_y}'
                do row = 1,N-1
                    p = row
                    if(row < ref) then
                        if (workstr%str(p)=='0') then
                            this%Circuit(row,ct) = '& \ctrlo{1}'
                        else if(workstr%str(p)=='1') then
                            this%Circuit(row,ct) = '& \ctrl{1}'
                        else
                            ! Case where workstr%str(p) = '*'
                            if(row > 1 .and. workstr%isempty(lb,p-1)==0) then
                                this%Circuit(row,ct) = '& \qw \qwx[1]'
                            else
                                this%Circuit(row,ct) = '& \qw'
                            end if
                        end if
                    else
                        if(workstr%str(p)=='0') then
                            this%Circuit(row+1,ct) = '& \ctrlo{-1}'
                        else if(workstr%str(p)=='1') then
                            this%Circuit(row+1,ct) = '& \ctrl{-1}'
                        else
                            ! Case where workstr%str(p) = '*'
                            if(row < N-1 .and. workstr%isempty(p+1,rb)==0) then
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
else
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
        lb = N+1-i
        call this%GroupGates(2**(i-1),this%GATEZ(:,M+N-i),0.0d0)
        call this%ReduceGroups()
        do j = 1, this%N_Type
            do k = 1, this%N_Per_Type(j)
                ct = ct+1
                if(i > 1) then
                    call workstr%copy(this%C_Num_Bin(j,k))
                    do row = 1, i-1
                        p = row
                        if(workstr%str(p)=='0') then
                            this%Circuit(row,ct) = '& \ctrlo{1}'
                        else if(workstr%str(p)=='1') then
                            this%Circuit(row,ct) = '& \ctrl{1}'
                        else
                            ! Case where workstr%str(p) = '*'
                            if(row > 1 .and. workstr%isempty(lb,p-1)==0) then
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
    lb = 1
    rb = N-1
    do i = M-1, 1, -1
        ref = this%index_level(i)
        call this%GroupGates(this%Mh,this%GATEY(:,i),0.0d0)
        call this%ReduceGroups()
        do j = 1, this%N_Type
            do k = 1, this%N_Per_Type(j)
                ct = ct+1
                call workstr%copy(this%C_Num_Bin(j,k))
                this%Circuit(ref,ct) = '& \gate{R_y}'
                do row = 1,N-1
                    p = row
                    if(row < ref) then
                        ! Qubits before target
                        if (workstr%str(p)=='0') then
                            this%Circuit(row,ct) = '& \ctrlo{1}'
                        else if(workstr%str(p)=='1') then
                            this%Circuit(row,ct) = '& \ctrl{1}'
                        else
                            ! Case where workstr%str(p) = '*'
                            if(row > 1 .and. workstr%isempty(lb,p-1)==0) then
                                this%Circuit(row,ct) = '& \qw \qwx[1]'
                            else
                                this%Circuit(row,ct) = '& \qw'
                            end if
                        end if
                    else
                        ! Qubits after target
                        if(workstr%str(p)=='0') then
                            this%Circuit(row+1,ct) = '& \ctrlo{-1}'
                        else if(workstr%str(p)=='1') then
                            this%Circuit(row+1,ct) = '& \ctrl{-1}'
                        else
                            ! Case where workstr%str(p) = '*'
                            if(row < N-1 .and. workstr%isempty(p+1,rb)==0) then
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
                call workstr%copy(this%C_Num_Bin(j,k))
                this%Circuit(ref,ct) = '& \gate{R_z}'
                do row = 1,N-1
                    p = row
                    if(row < ref) then
                        ! Qubits before target
                        if (workstr%str(p)=='0') then
                            this%Circuit(row,ct) = '& \ctrlo{1}'
                        else if(workstr%str(p)=='1') then
                            this%Circuit(row,ct) = '& \ctrl{1}'
                        else
                            ! Case where workstr%str(p) = '*'
                            if(row > 1 .and. workstr%isempty(lb,p-1)==0) then
                                this%Circuit(row,ct) = '& \qw \qwx[1]'
                            else
                                this%Circuit(row,ct) = '& \qw'
                            end if
                        end if
                    else
                        ! Qubits after target
                        if(workstr%str(p)=='0') then
                            this%Circuit(row+1,ct) = '& \ctrlo{-1}'
                        else if(workstr%str(p)=='1') then
                            this%Circuit(row+1,ct) = '& \ctrl{-1}'
                        else
                            ! Case where workstr%str(p) = '*'
                            if(row < N-1 .and. workstr%isempty(p+1,rb)==0) then
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
end if
! Just enough precision to look nice in the circuit plot
302 format('& ',F7.4)
! Excessive precision for result checking - looks awful in circuit pllot
!302 format('& ',F12.9)
this%csdr_ct = ct
call workstr%destructor()

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
        call this%C_Num_bin(match_idx,this%N_Per_Type(match_idx))%getbinrep(i-1)
        this%N_Total = this%N_Total + 1
    end if
end do
return

end subroutine csd_generator_GroupGates

subroutine csd_generator_ReduceGroups(this)

implicit none
class(csd_generator) :: this

integer :: N
integer :: i, j, k, row, p, limit

N = this%N
do i = 1, this%N_Type
    if(this%N_Per_Type(i) == 1) cycle
    do row = 1, N-1
        j = 1
        do while(j < this%N_Per_Type(i))
            limit = this%N_Per_Type(i)
            do k = j+1, limit
                if(this%C_Num_Bin(i,j)%qsimilar(this%C_Num_Bin(i,k),row) == 1) then
                    p = row
                    this%C_Num_Bin(i,j)%str(p) = '*'
                    if(k /= limit) call this%C_Num_Bin(i,k)%copy(this%C_Num_Bin(i,limit))
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

this%fname = fname
! Select a random unit number (10 or above)
this%FN = rng_inst%rint(1000)
do while (this%FN < 10)
    this%FN = rng_inst%rint(1000)
end do

end subroutine csd_write_handle_set_file

subroutine csd_write_handle_assign_target(this,targ)

implicit none
class(csd_write_handle) :: this
class(csd_solution), target :: targ

this%Circuit => targ%Circuit
this%nct = targ%csdr_ct

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

end module csd_tools
