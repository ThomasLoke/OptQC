! Copied code from CYG_INDEXTABLE

subroutine CYG_INDEXTABLE(N, M, index_level, index_pair)

!
!	FUNCTION
!	========
!
!	1. Create an array "index_level"
!		"index_level" is used to decide the recursive level of each decomposed matrix in the matrix sequence A and B.
!		index_table(2**(N-1)*(2j-1)) = i (i = 1,...,N; j = 1,...,2**(i-1))
!		e.g.
!		if N = 3, then
!			index_table = 3,2,3,1,3,2,3
!		if N = 4, then
!			index_table = 4,3,4,2,4,3,4,1,4,3,4,2,4,3,4
!	2. Create a 3-dimension table "index_level"
!		"index_pair" is used to decide the structure of each decomposed matrix in the matrix sequence A and B and helps with the process afterwards.
!		e.g.
!		if N = 3, then
!			i = 1, index_pair(:,:,1) =
!				[ 1,5 ]
!				[ 2,6 ]
!				[ 3,7 ]
!				[ 4,8 ]
!			i = 2, index_pair(:,:,2) =
!				[ 1,3 ]
!				[ 2,4 ]
!				[ 5,7 ]
!				[ 6,8 ]
!			i = 3, index_pair(:,:,3) =
!				[ 1,2 ]
!				[ 3,4 ]
!				[ 5,6 ]
!				[ 7,8 ]
!		if N = 4, then
!			i = 1, index_pair(:,:,1) =
!				[  1, 9 ]
!				[  2,10 ]
!				[  3,11 ]
!				[  4,12 ]
!				[  5,13 ]
!				[  6,14 ]
!				[  7,15 ]
!				[  8,16 ]
!			i = 2, index_pair(:,:,2) =
!				[  1, 5 ]
!				[  2, 6 ]
!				[  3, 7 ]
!				[  4, 8 ]
!				[  9,13 ]
!				[ 10,14 ]
!				[ 11,15 ]
!				[ 12,16 ]
!			i = 3, index_pair(:,:,3) =
!				[  1, 3 ]
!				[  2, 4 ]
!				[  5, 7 ]
!				[  6, 8 ]
!				[  9,11 ]
!				[ 10,12 ]
!				[ 13,15 ]
!				[ 14,16 ]
!			i = 4, index_pair(:,:,4) =
!				[  1, 2 ]
!				[  3, 4 ]
!				[  5, 6 ]
!				[  7, 8 ]
!				[  9,10 ]
!				[ 11,12 ]
!				[ 13,14 ]
!				[ 15,16 ]
!
!	ARGUMENT
!	========
!
!	N: (Input) Integer
!		The number of qubits for the quantum circuit
!	M: (Input) Integer
!		The size of the unitary matrix X. M = 2**N.
!	index_level: (Input & Output) Integer, dimension(M-1)
!		"index_level" is used to decide the recursive level of each decomposed matrix in the matrix sequence A and B.
!	index_pair: (Input & Output) Integer, dimension(M/2,2,N)
!		"index_pair" is used to decide the structure of each decomposed matrix in the matrix sequence A and B and helps with the process afterwards.
!

integer, intent(in)    :: N,M
integer, intent(inout) :: index_level(M-1)
integer, intent(inout) :: index_pair(M/2,2,N)

integer                :: i,j,k
integer                :: length, height
integer, allocatable   :: index_pair_i(:,:), index_line(:), index_layer(:,:,:)

allocate(index_pair_i(M/2,2))
index_pair_i(:,:) = 0
allocate(index_line(M))
index_line(:) = 0
!	Create a consecutive integer array "index_line" (Dimension(M)) from 1 to M
do j = 1,M
	index_line(j) = j
end do

do i = 1, N
	do j = 1, 2**(i-1)
    	!	Create the array "index_level"
		index_level(2**(N-i)*(2*j-1)) = i
	end do

	length = 2**(N-i)
	height = 2**(i-1)
	allocate(index_layer(length,2,height))
	index_layer(:,:,:) = 0
    !	For each "i", Reshape the array "index_line" (Dimension(M)) into the 3D table "index_layer"  (Dimension(length,2,height))
    !	e.g. N = 3
    !		index_line = [1,2,3,4,5,6,7,8]
    !		If i = 1, then
    !		index_layer = (Dimension(4,2,1))
    !			[1,5]
    !			[2,6]
    !			[3,7]
    !			[4,8]
    !		If i = 2, then
    !		index_layer = (Dimension(2,2,2))
    !			[1,3]	[5,7]
    !			[2,4]	[6,8]
    !		If i = 3, then
    !		index_layer = (Dimension(1,2,4))
    !			[1,2]	[3,4]	[5,6]	[7,8]
	index_layer = reshape(index_line,(/length,2,height/))
	index_pair_i(:,:) = 0
    !	Link index_layer and create index_pair_i
    !	e.g. N = 3
    !		If i = 1, then
    !			index_layer =
    !				[1,5]
    !				[2,6]
    !				[3,7]
    !				[4,8]
    !			index_pair_i =
    !				[1,5]
    !				[2,6]
    !				[3,7]
    !				[4,8]
    !		If i = 2, then
    !			index_layer =
    !				[1,3]	[5,7]
    !				[2,4]	[6,8]
    !			index_pair_i =
    !				[1,3]
    !				[2,4]
    !				[5,7]
    !				[6,8]
    !		If i = 3, then
    !			index_layer =
    !				[1,2]	[3,4]	[5,6]	[7,8]
    !			index_pair_i =
    !				[1,2]
    !				[3,4]
    !				[5,6]
    !				[7,8]
	do k = 1,height
		index_pair_i(length*(k-1)+1:length*k,:) = index_layer(:,:,k)
	end do
    !	Create "index_pair", which satisfies "index_pair(:,:,i) = index_pair_i"
	index_pair(:,:,i) = index_pair_i
	deallocate(index_layer)
end do
deallocate(index_pair_i, index_line)
return

end subroutine CYG_INDEXTABLE

! Copied code from CYGC_COEFF

subroutine CYGC_COEFF(N, M, COEFF)

!
!	FUNCTION
!	========
!
!	To obtain the equivalent quantum gates of Z(2**N), namely the last N Rz Gates and a Phase Gate, we need to solve the following linear equation:
!
!	e.g. N = 3
!		[  1  0  0  0 |  1  0 |  1 | 0 ]  [ GateZ'_{3,1} ]   [ PHASEZ_{8,1} ]
!		[ -1  0  0  0 |  1  0 |  1 | 0 ]  [ GateZ'_{3,2} ]   [ PHASEZ_{8,2} ]
!		[  0  1  0  0 | -1  0 |  1 | 0 ]  [ GateZ'_{3,3} ]   [ PHASEZ_{8,3} ]
!		[  0 -1  0  0 | -1  0 |  1 | 0 ]  [ GateZ'_{3,4} ] = [ PHASEZ_{8,4} ]
!		[  0  0  1  0 |  0  1 | -1 | 1 ]  [ GateZ'_{2,1} ]   [ PHASEZ_{8,5} ]
!		[  0  0 -1  0 |  0  1 | -1 | 1 ]  [ GateZ'_{2,2} ]   [ PHASEZ_{8,6} ]
!		[  0  0  0  1 |  0 -1 | -1 | 1 ]  [ GateZ'_{1,1} ]   [ PHASEZ_{8,7} ]
!		[  0  0  0 -1 |  0 -1 | -1 | 1 ]  [  GatePhase   ]   [ PHASEZ_{8,8} ]
!
!	e.g. N = 4
!		[  1  0  0  0  0  0  0  0 |  1  0  0  0 |  1  0 |  1 | 0 ]  [ GateZ'_{4,1} ]   [ PHASEZ_{16, 1} ]
!		[ -1  0  0  0  0  0  0  0 |  1  0  0  0 |  1  0 |  1 | 0 ]  [ GateZ'_{4,2} ]   [ PHASEZ_{16, 2} ]
!		[  0  1  0  0  0  0  0  0 | -1  0  0  0 |  1  0 |  1 | 0 ]  [ GateZ'_{4,3} ]   [ PHASEZ_{16, 3} ]
!		[  0 -1  0  0  0  0  0  0 | -1  0  0  0 |  1  0 |  1 | 0 ]  [ GateZ'_{4,4} ]   [ PHASEZ_{16, 4} ]
!		[  0  0  1  0  0  0  0  0 |  0  1  0  0 | -1  0 |  1 | 0 ]  [ GateZ'_{4,5} ]   [ PHASEZ_{16, 5} ]
!		[  0  0 -1  0  0  0  0  0 |  0  1  0  0 | -1  0 |  1 | 0 ]  [ GateZ'_{4,6} ]   [ PHASEZ_{16, 6} ]
!		[  0  0  0  1  0  0  0  0 |  0 -1  0  0 | -1  0 |  1 | 0 ]  [ GateZ'_{4,7} ]   [ PHASEZ_{16, 7} ]
!		[  0  0  0 -1  0  0  0  0 |  0 -1  0  0 | -1  0 |  1 | 0 ]  [ GateZ'_{4,8} ] = [ PHASEZ_{16, 8} ]
!		[  0  0  0  0  1  0  0  0 |  0  0  1  0 |  0  1 | -1 | 1 ]  [ GateZ'_{3,1} ]   [ PHASEZ_{16, 9} ]
!		[  0  0  0  0 -1  0  0  0 |  0  0  1  0 |  0  1 | -1 | 1 ]  [ GateZ'_{3,2} ]   [ PHASEZ_{16,10} ]
!		[  0  0  0  0  0  1  0  0 |  0  0 -1  0 |  0  1 | -1 | 1 ]  [ GateZ'_{3,3} ]   [ PHASEZ_{16,11} ]
!		[  0  0  0  0  0 -1  0  0 |  0  0 -1  0 |  0  1 | -1 | 1 ]  [ GateZ'_{3,4} ]   [ PHASEZ_{16,12} ]
!		[  0  0  0  0  0  0  1  0 |  0  0  0  1 |  0 -1 | -1 | 1 ]  [ GateZ'_{2,1} ]   [ PHASEZ_{16,13} ]
!		[  0  0  0  0  0  0 -1  0 |  0  0  0  1 |  0 -1 | -1 | 1 ]  [ GateZ'_{2,2} ]   [ PHASEZ_{16,14} ]
!		[  0  0  0  0  0  0  0  1 |  0  0  0 -1 |  0 -1 | -1 | 1 ]  [ GateZ'_{1,1} ]   [ PHASEZ_{16,15} ]
!		[  0  0  0  0  0  0  0 -1 |  0  0  0 -1 |  0 -1 | -1 | 1 ]  [  GatePhase   ]   [ PHASEZ_{16,16} ]
!
!	CYGC_COEFF create the coefficient matrix.
!
!	ARGUMENT
!	========
!
!	N: (Input) Integer
!		The number of qubits for the quantum circuit
!	M: (Input) Integer
!		The size of the unitary matrix X. M = 2**N.
!	COEFF: (Input & Output) Double Precision, dimension(M,M)
!		The coefficient matrix of the linear equation.
!

implicit none
integer :: N, M
double precision :: COEFF(M,M)

integer :: i, j

COEFF = 0.0d0
do i = 1,N
	do j = 1,2**(N-i)
		COEFF(2**i*(j-1)+1:2**i*j-2**(i-1),M-2**(N+1-i)+j) = 1.0d0
		COEFF(2**i*j-2**(i-1)+1:2**i*j,M-2**(N+1-i)+j) = -1.0d0
	end do
end do
COEFF(M/2+1:M,M) = 1.0d0
return

end subroutine CYGC_COEFF

! Functions/Subroutines for Main

function CalcTol(TOL_COEFF,ecur) result(res)

implicit none
double precision :: TOL_COEFF
integer :: ecur

integer :: res

res = ceiling(TOL_COEFF * ecur)
return

end function CalcTol

function FindMinPos(p,arr) result(res)

implicit none
integer :: p
integer :: arr(p)

integer :: i, curmin
integer :: res

res = 1
curmin = arr(1)
do i = 2, p
    if(arr(i) < curmin) then
        curmin = arr(i)
        res = i
    end if
end do
return

end function FindMinPos

function ChooseN(d) result(res)

implicit none
integer :: d
double precision :: temp
integer :: nceil, nround

integer :: res

temp = log(real(d,8))/log(2.0d0)
nceil = ceiling(temp)
nround = nint(temp)
res = nceil
! Dangerous case
if(d == 2**nround) then
    res = nround
end if
return

end function ChooseN

! Functions/Subroutines for Perm

! Find the first position (from the left) where an element of used is zero.
function FindUnused(M,used)

implicit none
integer :: M
integer :: used(M)

integer :: i, FindUnused

FindUnused = 1
do i = 1, M
    if(used(i) == 0) exit
    FindUnused = FindUnused + 1
end do

return

end function FindUnused

subroutine permlisttomatrix(M,PermList,PermMat)

implicit none
integer :: M
integer :: PermList(M)
double precision :: PermMat(M,M)

integer :: i

PermMat = 0.0d0
do i = 1, M
    PermMat(i,PermList(i)) = 1.0d0
end do

end subroutine permlisttomatrix

subroutine permlisttomatrixtr(M,PermList,PermMatTr)

implicit none
integer :: M
integer :: PermList(M)
double precision :: PermMatTr(M,M)

integer :: i

PermMatTr = 0.0d0
do i = 1, M
    PermMatTr(PermList(i),i) = 1.0d0
end do

end subroutine permlisttomatrixtr

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

subroutine qpermtoperm(N,M,qperm,perm)

implicit none
integer :: N, M
integer :: qperm(N), perm(M)

logical, allocatable :: bitval(:)
integer :: i, j, temp

! Allocate temporary variables
allocate(bitval(N))
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
deallocate(bitval)

end subroutine qpermtoperm

subroutine invertperm(len,Perm,PermTarg)

implicit none
integer :: len
integer :: Perm(len), PermTarg(len)

integer :: i

do i = 1, len
    PermTarg(Perm(i)) = i
end do

end subroutine invertperm

! Generates a random permutation of length M
subroutine perm_generate(M,perm)

use rng

implicit none
integer :: M
integer :: perm(M)

integer :: i, idx, temp

! Start with the identity permutation
do i = 1, M
    perm(i) = i
end do
! Choose each element randomly and fix the chosen ones from the left of the array
do i = 1, M-1
    idx = i-1+rng_inst%rint(M-i+1)
    if(i /= idx) then
        temp = perm(i)
        perm(i) = perm(idx)
        perm(idx) = temp
    end if
end do

end subroutine perm_generate