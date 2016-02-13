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
call SwapElem(M,csdss_source%arr(2)%X,csdss_targ%arr(2)%X,s1,s2)
! Reflect change in the permutation
temp = Perm(s1)
Perm(s1) = Perm(s2)
Perm(s2) = temp
! Update permutation in target solution set
call permlisttomatrixtr(M,Perm,csdss_targ%arr(1)%X)
call permlisttomatrix(M,Perm,csdss_targ%arr(3)%X)

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
