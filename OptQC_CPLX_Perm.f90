subroutine NeighbourhoodOpt_CPLX(M,csdss_source,csdss_targ,Perm)

use csd_cplx

implicit none
type(csd_solution_set_cplx) :: csdss_source, csdss_targ
integer :: M
integer :: Perm(M)

integer :: s1, s2, temp
integer :: RINT

M = csdss_source%M
s1 = RINT(M)
s2 = RINT(M)
do while (s2 == s1)
    s2 = RINT(M)
end do
! Perform swap between elements s1 and s2
call SwapElem_CPLX(M,csdss_source%arr(3)%X,csdss_targ%arr(3)%X,s1,s2)
! Reflect change in the permutation
temp = Perm(s1)
Perm(s1) = Perm(s2)
Perm(s2) = temp
! Update permutation in target solution set
call permlisttomatrixtr_CPLX(M,Perm,csdss_targ%arr(2)%X)
call permlisttomatrix_CPLX(M,Perm,csdss_targ%arr(4)%X)

end subroutine NeighbourhoodOpt_CPLX

subroutine SwapElem_CPLX(M,X,Z,s1,s2)

implicit none
integer :: M
double complex :: X(M,M), Z(M,M)
integer :: s1, s2

integer :: i

Z = X
do i = 1, M
    if(i /= s1 .and. i /= s2) then
        Z(i,s1) = X(i,s2)
        Z(i,s2) = X(i,s1)
        Z(s1,i) = X(s2,i)
        Z(s2,i) = X(s1,i)
    else
        if(i == s1) then
            Z(s1,s1) = X(s2,s2)
            Z(s2,s2) = X(s1,s1)
        end if
        if(i == s2) then
            Z(s2,s1) = X(s1,s2)
            Z(s1,s2) = X(s2,s1)
        end if
    end if
end do

return

end subroutine SwapElem_CPLX

subroutine ApplyPerm_CPLX(M,X,Z,Perm)

implicit none
integer :: M
double complex :: X(M,M), Z(M,M)
integer :: Perm(M)

integer :: i, j, row

do i = 1, M
    row = Perm(i)
    do j = 1, M
        Z(i,j) = X(row,Perm(j))
    end do
end do

return

end subroutine ApplyPerm_CPLX

subroutine permlisttomatrix_CPLX(M,PermList,PermMat)

implicit none
integer :: M
integer :: PermList(M)
double complex :: PermMat(M,M)

integer :: i, j, temp

do i = 1, M
    temp = PermList(i)
    do j = 1, M
        if(j == temp) then
            PermMat(i,j) = cmplx(1.0d0)
        else
            PermMat(i,j) = cmplx(0.0d0)
        end if
    end do
end do

end subroutine permlisttomatrix_CPLX

subroutine permlisttomatrixtr_CPLX(M,PermList,PermMatTr)

implicit none
integer :: M
integer :: PermList(M)
double complex :: PermMatTr(M,M)

integer :: i, j, temp

do i = 1, M
    temp = PermList(i)
    do j = 1, M
        if(j == temp) then
            PermMatTr(j,i) = cmplx(1.0d0)
        else
            PermMatTr(j,i) = cmplx(0.0d0)
        end if
    end do
end do

end subroutine permlisttomatrixtr_CPLX

subroutine qperm_compute_CPLX(N,M,csd_obj,qperm)

use csd_cplx

implicit none
integer :: N, M
type(csd_solution_cplx) :: csd_obj
integer :: qperm(N)

integer, allocatable :: qperm_temp(:), perm(:)
logical, allocatable :: bitval(:)
integer :: i, j, pt1, pt2, temp, ct

! Allocate temporary variables - this subroutine should only be run once per thread
allocate(qperm_temp(N))
allocate(perm(M))
allocate(bitval(N))
! Initialize qperm_temp to the identity qubit permutation
do i = 1, N
    qperm_temp(i) = i
end do
call csd_obj%clean()
! Determine associated state permutation
do i = 1, M
    bitval = .false.
    ! Find bit string representation
    temp = i-1
    do j = 1, N
        bitval(N+1-j) = btest(temp,j-1)
    end do
    ! Determine new value after swapping bits from old bit string directly
    temp = 0
    do j = 1, N
        if(bitval(qperm(j)) == .true.) then
            temp = temp + (2**(N-j))
        end if
    end do
    perm(i) = temp + 1
end do
! Construct matrix from the state permutation
! Note: Transpose of permutation - MAGIC DON'T TOUCH PLOX
call permlisttomatrixtr_CPLX(M,perm,csd_obj%X)
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
deallocate(bitval)

end subroutine qperm_compute_CPLX

subroutine qperm_process_CPLX(N,M,csdss_obj,csdgen_obj,QPerm,X,ecur)

use csd_cplx

implicit none
integer :: N, M
type(csd_solution_set_cplx) :: csdss_obj
type(csd_generator_cplx) :: csdgen_obj
integer :: QPerm(N)
double complex :: X(M,M)
integer :: ecur

! Construct state permutation matrix from qubit permutation
call qperm_compute_CPLX(N,M,csdss_obj%arr(5),QPerm)                              ! Q
call qperm_reverse_CPLX(N,csdss_obj%arr(5),csdss_obj%arr(1))                     ! Q^T
! Note: Q U Q^T = P^T U' P, so we treat Q U Q^T as the matrix to be decomposed
csdss_obj%arr(3)%X = matmul(csdss_obj%arr(5)%X,X)                           ! Q U
csdss_obj%arr(3)%X = matmul(csdss_obj%arr(3)%X,csdss_obj%arr(1)%X)          ! Q U Q^T
! Count initial number of gates (including reduction)
call csdss_obj%arr(3)%run_csdr(csdgen_obj)
csdss_obj%csd_ss_ct = csdss_obj%arr(3)%csd_ct + csdss_obj%arr(1)%csd_ct + csdss_obj%arr(5)%csd_ct
csdss_obj%csdr_ss_ct = csdss_obj%arr(3)%csdr_ct + csdss_obj%arr(1)%csdr_ct + csdss_obj%arr(5)%csdr_ct
ecur = csdss_obj%csdr_ss_ct

end subroutine qperm_process_CPLX

subroutine qperm_reverse_CPLX(N,csd_obj_source,csd_obj_targ)

use csd_cplx

implicit none
integer :: N
type(csd_solution_cplx) :: csd_obj_source, csd_obj_targ

integer :: i, j, idx, ct

call csd_obj_targ%clean()
ct = csd_obj_source%csd_ct
csd_obj_targ%X = transpose(csd_obj_source%X)
do i = 1, ct
    idx = ct+1-i
    do j = 1, N
        csd_obj_targ%Circuit(j,idx) = csd_obj_source%Circuit(j,i)
    end do
end do
csd_obj_targ%csd_ct = ct
csd_obj_targ%csdr_ct = ct

end subroutine qperm_reverse_CPLX
