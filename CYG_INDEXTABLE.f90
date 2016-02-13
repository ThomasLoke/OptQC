
SUBROUTINE CYG_INDEXTABLE(N, M, index_level, index_pair)

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

INTEGER, INTENT(IN)    :: N,M
INTEGER, INTENT(INOUT) :: index_level(M-1)
INTEGER, INTENT(INOUT) :: index_pair(M/2,2,N)

INTEGER                :: i,j,k
INTEGER                :: length, height
INTEGER, ALLOCATABLE   :: index_pair_i(:,:), index_line(:), index_layer(:,:,:)

ALLOCATE(index_pair_i(M/2,2))
index_pair_i(:,:) = 0
ALLOCATE(index_line(M))
index_line(:) = 0
!	Create a consecutive integer array "index_line" (Dimension(M)) from 1 to M
DO j = 1,M
	index_line(j) = j
END DO

DO i = 1, N
	DO j = 1, 2**(i-1)
    	!	Create the array "index_level"
		index_level(2**(N-i)*(2*j-1)) = i
	END DO

	length = 2**(N-i)
	height = 2**(i-1)
	ALLOCATE(index_layer(length,2,height))
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
	DO k = 1,height
		index_pair_i(length*(k-1)+1:length*k,:) = index_layer(:,:,k)
	END DO
    !	Create "index_pair", which satisfies "index_pair(:,:,i) = index_pair_i"
	index_pair(:,:,i) = index_pair_i
	DEALLOCATE(index_layer)
END DO
DEALLOCATE(index_pair_i, index_line)

RETURN
END
