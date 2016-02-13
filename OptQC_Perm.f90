module csd_perm

use csd_tools
use rng

implicit none
private :: GGset, GGsetdim, perm_temp, perm_temp2, binstr_temp
save :: GGset, GGsetdim, perm_temp, perm_temp2, binstr_temp
character, allocatable :: GGset(:)
integer :: GGsetdim
integer, allocatable :: perm_temp(:), perm_temp2(:)
type(binstr) :: binstr_temp

contains

subroutine initialize_csd_perm(N,M)

implicit none
integer :: N, M

GGsetdim = 4
allocate(GGset(GGsetdim))
GGset(1:GGsetdim) = (/ '0', '1', '*', 't' /)
allocate(perm_temp(M))
allocate(perm_temp2(M))
perm_temp = 0
perm_temp2 = 0
call binstr_temp%constructor(N)

end subroutine initialize_csd_perm

subroutine finalize_csd_perm()

implicit none

GGsetdim = 0
deallocate(GGset)
deallocate(perm_temp)
deallocate(perm_temp2)
call binstr_temp%destructor()

end subroutine finalize_csd_perm

function gatecounttarg(gstr)

implicit none
type(binstr) :: gstr
integer :: gatecounttarg

integer :: i

gatecounttarg = 0
do i = 1, gstr%len
    if(gstr%str(i) == 't') then
        gatecounttarg = gatecounttarg + 1
    end if
end do

end function gatecounttarg

subroutine generategate(gstr)

implicit none
type(binstr) :: gstr

integer :: i, temp

! Generates a gate with at least one NOT gate
gstr%str = ''
! Generate with the constraint that the number of t's is >= 1
do while(gatecounttarg(gstr) < 1)
    do i = 1, gstr%len
        temp = rng_inst%rint(GGsetdim)
        gstr%str(i) = GGset(temp)
    end do
end do

end subroutine generategate

! subroutine generategate(gstr)

! implicit none
! type(binstr) :: gstr

! integer :: i, temp, pos

! ! Generates a gate with only one NOT gate
! gstr%str = ''
! pos = rng_inst%rint(gstr%len)
! gstr%str(pos) = 't'
! ! Generate the rest of the gate
! do i = 1, gstr%len
!     if(i /= pos) then
!         temp = rng_inst%rint(GGsetdim-1)
!         gstr%str(i) = GGset(temp)
!     end if
! end do

! end subroutine generategate

subroutine getpermofgstr(N,M,gstr,Perm)

implicit none
integer :: N, M
type(binstr) :: gstr
integer :: Perm(M)

integer :: i, j
character :: gsec, ssec
logical :: flag

Perm = 0
do i = 1, M
    Perm(i) = i
    call binstr_temp%getbinrep(i-1)
    ! If flag is .true., then use the new bit string value as the target; otherwise, if flag is .false., then the action is the identity (presumably due to not matching conditional)
    flag = .true.
    do j = 1, N
        gsec = gstr%str(j)
        ssec = binstr_temp%str(j)
        ! If conditional is not matched, then the gate action is the identity
        if(gsec == '0' .or. gsec == '1') then
            if(gsec /= ssec) then
                flag = .false.
                exit
            end if
        ! If target, flip bit
        else if(gsec == 't') then
            if(ssec == '0') then
                binstr_temp%str(j) = '1'
            else
                binstr_temp%str(j) = '0'
            end if
        ! The only other case is gsec = '*', which is skipped
        else
            cycle
        end if
    end do
    ! WARNING: Potentially assumed the wrong convention for permutations (that is, could be transposed or not) - problem should only arise for non-self-invertible gates, which shoudn't be used anyway
    if(flag == .true.) then
        Perm(i) = binstr_temp%getdecrep()+1
    end if
end do

end subroutine getpermofgstr

subroutine findextentgstr(gstr,l,r)

implicit none
type(binstr) :: gstr
integer :: l, r

integer :: i, n
logical :: flag

l = 0
r = 0
n = gstr%len
flag = .false.
! Note that if l = 0, then r = 0 (and vice versa)
! Check if there are any conditonals - if there aren't any, then the (conditional) gate extent is defined to be zero
do i = 1, n
    if(gstr%str(i) == '0' .or. gstr%str(i) == '1') then
        flag = .true.
        exit
    end if
end do
! Calculate the extent
if(flag == .true.) then
    do i = 1, n
        if(gstr%str(i) == '0' .or. gstr%str(i) == '1' .or. gstr%str(i) == 't') then
            if(l == 0) then
                l = i
                exit
            end if
        end if
    end do
    do i = n, 1, -1
         if(gstr%str(i) == '0' .or. gstr%str(i) == '1' .or. gstr%str(i) == 't') then
            if(r == 0) then
                r = i
                exit
            end if
        end if
    end do
end if

end subroutine findextentgstr

subroutine NeighbourhoodOpt(N,M,csdss_Xinit,csdss_source,csdss_targ,Perm)

implicit none
type(csd_solution_set) :: csdss_Xinit, csdss_source, csdss_targ
integer :: N, M
integer :: Perm(M)

type(binstr) :: gstr
integer :: i, col, l, r, idx

call gstr%constructor(N)
! Check first if circuit is close to hardcoded limit of M*M+1, if so, perform a reset - if not just transfer the circuit and gate count as normal
if(csdss_source%arr(2)%csdr_ct == M*M) then
    do i = 1, M
        Perm(i) = i
    end do
    call csdss_targ%arr(2)%clean()
    call csdss_targ%arr(4)%clean()
else
    ! Note: Only need to copy the circuit for P - P' is derived from P later
    csdss_targ%arr(4)%Circuit = csdss_source%arr(4)%Circuit
    csdss_targ%arr(2)%csd_ct = csdss_source%arr(2)%csd_ct
    csdss_targ%arr(2)%csdr_ct = csdss_source%arr(2)%csdr_ct
    csdss_targ%arr(4)%csd_ct = csdss_source%arr(4)%csd_ct
    csdss_targ%arr(4)%csdr_ct = csdss_source%arr(4)%csdr_ct
end if
! Generate a random permutation gate
call generategate(gstr)
call findextentgstr(gstr,l,r)
perm_temp = 0
call getpermofgstr(N,M,gstr,perm_temp)
perm_temp2 = Perm
! Compute the new permutation
! perm_temp1 = pg, perm_temp2 = pn
! Then matrixform[pnew] = matrixform[pg].matrixform[pn] is equivalent to pnew(i) = pn(pg(i))
do i = 1, M
    Perm(i) = perm_temp2(perm_temp(i))
end do
! Add the random permutation gate to the solution set
! Update the permutation matrices
call permlisttomatrixtr(M,Perm,csdss_targ%arr(2)%X)
call permlisttomatrix(M,Perm,csdss_targ%arr(4)%X)
! Update the U' matrix
if(csdss_targ%arr(3)%obj_type == 0) then
    call ApplyPerm(M,csdss_Xinit%arr(3)%X,csdss_targ%arr(3)%X,Perm)
else
    call ApplyPerm_CPLX(M,csdss_Xinit%arr(3)%Xc,csdss_targ%arr(3)%Xc,Perm)
end if
! Update the gate count
csdss_targ%arr(2)%csd_ct = csdss_targ%arr(2)%csd_ct + 1
csdss_targ%arr(2)%csdr_ct = csdss_targ%arr(2)%csdr_ct + 1
csdss_targ%arr(4)%csd_ct = csdss_targ%arr(4)%csd_ct + 1
csdss_targ%arr(4)%csdr_ct = csdss_targ%arr(4)%csdr_ct + 1
! Update the circuit for P
col = csdss_targ%arr(4)%csdr_ct
do i = 1, N
    select case (gstr%str(i))
        case('0')
            if(i < r) then
                csdss_targ%arr(4)%Circuit(i,col) = '& \ctrlo{1}'
            else
                csdss_targ%arr(4)%Circuit(i,col) = '& \ctrlo{-1}'
            end if
        case('1')
            if(i < r) then
                csdss_targ%arr(4)%Circuit(i,col) = '& \ctrl{1}'
            else
                csdss_targ%arr(4)%Circuit(i,col) = '& \ctrl{-1}'
            end if
        case('t')
            if(l == 0) then
                csdss_targ%arr(4)%Circuit(i,col) = '& \qswap'
            else
                if(i < r) then
                    csdss_targ%arr(4)%Circuit(i,col) = '& \targ \qw \qwx[1]'
                else
                    csdss_targ%arr(4)%Circuit(i,col) = '& \targ \qw'
                end if
            end if
        case default
            ! Should be just '*'
            if(i > l .and. i < r) then
                csdss_targ%arr(4)%Circuit(i,col) = '& \qw \qwx[1]'
            else
                csdss_targ%arr(4)%Circuit(i,col) = '& \qw'
            end if
    end select
end do
csdss_targ%arr(4)%Circuit(N+1,col) = '&'
! Copy reverse cicuit to P'
idx = col
do i = 1, col
    csdss_targ%arr(2)%Circuit(:,idx) = csdss_targ%arr(4)%Circuit(:,i)
    idx = idx - 1
end do
call gstr%destructor()

end subroutine NeighbourhoodOpt

subroutine qperm_compute(N,M,csd_obj,qperm)

implicit none
integer :: N, M
type(csd_solution) :: csd_obj
integer :: qperm(N)

integer, allocatable :: qperm_temp(:), perm(:)
integer :: i, j, pt1, pt2, temp, ct

! Allocate temporary variables - not terribly efficient but meh could be worse
allocate(qperm_temp(N))
allocate(perm(M))
! Initialize qperm_temp to the identity qubit permutation
do i = 1, N
    qperm_temp(i) = i
end do
call csd_obj%clean()
! Determine associated state permutation
call qpermtoperm(N,M,qperm,perm)
! Construct matrix from the state permutation
! Note: Transpose of permutation - MAGIC DON'T TOUCH PLOX
! Assumed to be of real type
call permlisttomatrixtr(M,perm,csd_obj%X)
! Determine swap gate ordering
ct = 0
do i = 1, N-1
    pt1 = i
    pt2 = -1
    temp = qperm(i)
    do j = i, N
        if(temp == qperm_temp(j)) then
            pt2 = j
            exit
        end if
    end do
    ! Reminder: pt2 >= pt1 assuming a valid qubit permutation
    if(pt2 /= -1) then
        qperm_temp(pt2) = qperm_temp(pt1)
        qperm_temp(pt1) = temp
        if(pt1 /= pt2) then
            ct = ct + 1
            do j = 1, N
                if(j == pt1) then
                    csd_obj%Circuit(j,ct) = "& \qswap \qwx[1]"
                else if(j > pt1 .and. j < pt2) then
                    csd_obj%Circuit(j,ct) = "& \qw \qwx[1]"
                else if(j == pt2) then
                    csd_obj%Circuit(j,ct) = "& \qswap"
                else
                    csd_obj%Circuit(j,ct) = "& \qw"
                end if
            end do
            csd_obj%Circuit(N+1,ct) = "&"
        end if
    else
        write(*,*)"Invalid qubit permutation: ",qperm
        stop
    end if
end do
! Update gate count in csd_obj
csd_obj%csd_ct = ct
csd_obj%csdr_ct = ct
! Deallocate temporary variables
deallocate(qperm_temp)
deallocate(perm)

end subroutine qperm_compute

subroutine qperm_process(N,M,csdss_obj,csdgen_obj,QPerm,X,ecur)

implicit none
integer :: N, M
type(csd_solution_set) :: csdss_obj
type(csd_generator) :: csdgen_obj
integer :: QPerm(N)
double precision :: X(M,M)
integer :: ecur

! Construct state permutation matrix from qubit permutation
call qperm_compute(N,M,csdss_obj%arr(5),QPerm)                              ! Q
call qperm_reverse(N,csdss_obj%arr(5),csdss_obj%arr(1))                     ! Q^T
! Note: Q U Q^T = P^T U' P, so we treat Q U Q^T as the matrix to be decomposed
csdss_obj%arr(3)%X = matmul(csdss_obj%arr(5)%X,X)                           ! Q U
csdss_obj%arr(3)%X = matmul(csdss_obj%arr(3)%X,csdss_obj%arr(1)%X)          ! Q U Q^T
! Count initial number of gates (including reduction)
call csdss_obj%arr(3)%run_csdr(csdgen_obj)
csdss_obj%csd_ss_ct = csdss_obj%arr(3)%csd_ct + csdss_obj%arr(1)%csd_ct + csdss_obj%arr(5)%csd_ct
csdss_obj%csdr_ss_ct = csdss_obj%arr(3)%csdr_ct + csdss_obj%arr(1)%csdr_ct + csdss_obj%arr(5)%csdr_ct
ecur = csdss_obj%csdr_ss_ct

end subroutine qperm_process

subroutine qperm_process_CPLX(N,M,csdss_obj,csdgen_obj,QPerm,X,ecur)

implicit none
integer :: N, M
type(csd_solution_set) :: csdss_obj
type(csd_generator) :: csdgen_obj
integer :: QPerm(N)
double complex :: X(M,M)
integer :: ecur

! Construct state permutation matrix from qubit permutation
call qperm_compute(N,M,csdss_obj%arr(5),QPerm)                              ! Q
call qperm_reverse(N,csdss_obj%arr(5),csdss_obj%arr(1))                     ! Q^T
! Note: Q U Q^T = P^T U' P, so we treat Q U Q^T as the matrix to be decomposed
csdss_obj%arr(3)%Xc = matmul(csdss_obj%arr(5)%X,X)                          ! Q U
csdss_obj%arr(3)%Xc = matmul(csdss_obj%arr(3)%Xc,csdss_obj%arr(1)%X)        ! Q U Q^T
! Count initial number of gates (including reduction)
call csdss_obj%arr(3)%run_csdr(csdgen_obj)
csdss_obj%csd_ss_ct = csdss_obj%arr(3)%csd_ct + csdss_obj%arr(1)%csd_ct + csdss_obj%arr(5)%csd_ct
csdss_obj%csdr_ss_ct = csdss_obj%arr(3)%csdr_ct + csdss_obj%arr(1)%csdr_ct + csdss_obj%arr(5)%csdr_ct
ecur = csdss_obj%csdr_ss_ct

end subroutine qperm_process_CPLX

subroutine qperm_reverse(N,csd_obj_source,csd_obj_targ)

implicit none
integer :: N
type(csd_solution) :: csd_obj_source, csd_obj_targ

integer :: i, j, idx, ct

call csd_obj_targ%clean()
ct = csd_obj_source%csd_ct
csd_obj_targ%X = transpose(csd_obj_source%X)
do i = 1, ct
    idx = ct+1-i
    csd_obj_targ%Circuit(:,idx) = csd_obj_source%Circuit(:,i)
end do
csd_obj_targ%csd_ct = ct
csd_obj_targ%csdr_ct = ct

end subroutine qperm_reverse

end module csd_perm