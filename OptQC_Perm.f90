subroutine NeighbourhoodOpt(M,csdss_source,csdss_targ,Perm)

use csd_real

implicit none
type(csd_solution_set) :: csdss_source, csdss_targ
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
call SwapElem(M,csdss_source%arr(3)%X,csdss_targ%arr(3)%X,s1,s2)
! Reflect change in the permutation
temp = Perm(s1)
Perm(s1) = Perm(s2)
Perm(s2) = temp
! Update permutation in target solution set
call permlisttomatrixtr(M,Perm,csdss_targ%arr(2)%X)
call permlisttomatrix(M,Perm,csdss_targ%arr(4)%X)

end subroutine NeighbourhoodOpt

subroutine SwapElem(M,X,Z,s1,s2)

implicit none
integer :: M
double precision :: X(M,M), Z(M,M)
integer :: s1, s2

integer :: i, j

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

end subroutine SwapElem

subroutine ApplyPerm(M,X,Z,Perm)

implicit none
integer :: M
double precision :: X(M,M), Z(M,M)
integer :: Perm(M)

integer :: i, j, row

do i = 1, M
    row = Perm(i)
    do j = 1, M
        Z(i,j) = X(row,Perm(j))
    end do
end do

return

end subroutine ApplyPerm

subroutine permlisttomatrix(M,PermList,PermMat)

implicit none
integer :: M
integer :: PermList(M)
double precision :: PermMat(M,M)

integer :: i, j, temp

do i = 1, M
    temp = PermList(i)
    do j = 1, M
        if(j == temp) then
            PermMat(i,j) = 1.0d0
        else
            PermMat(i,j) = 0.0d0
        end if
    end do
end do

end subroutine permlisttomatrix

subroutine permlisttomatrixtr(M,PermList,PermMatTr)

implicit none
integer :: M
integer :: PermList(M)
double precision :: PermMatTr(M,M)

integer :: i, j, temp

do i = 1, M
    temp = PermList(i)
    do j = 1, M
        if(j == temp) then
            PermMatTr(j,i) = 1.0d0
        else
            PermMatTr(j,i) = 0.0d0
        end if
    end do
end do

end subroutine permlisttomatrixtr

subroutine qperm_compute(N,M,csd_obj,qperm,my_rank)

use csd_real

implicit none
integer :: N, M
type(csd_solution) :: csd_obj
integer :: qperm(N)
integer(4) :: my_rank

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
deallocate(bitval)

end subroutine qperm_compute

subroutine qperm_generate(N,qperm)

implicit none
integer :: N
integer :: qperm(N)

integer :: i, idx, temp
integer :: RINT

! Start with the identity qubit permutation
do i = 1, N
    qperm(i) = i
end do
! Choose each element randomly and fix the chosen ones from the left of the array
do i = 1, N-1
    idx = i-1+RINT(N-i+1)
    temp = qperm(i)
    qperm(i) = qperm(idx)
    qperm(idx) = temp
end do

end subroutine qperm_generate

subroutine qperm_reverse(N,M,csd_obj_source,csd_obj_targ)

use csd_real

implicit none
integer :: N, M
type(csd_solution) :: csd_obj_source, csd_obj_targ

integer :: i, j, idx, ct

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

end subroutine qperm_reverse
