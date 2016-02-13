SUBROUTINE CYGR_CUTGATE (Z, GATEPI, GATEY, sgn, N, M, index_level, index_pair, Z_array, PI_array, GATEY_sign)

!
!	FUNCTION
!	========
!
!	The ultimate purpose of the SUBROUTINE CYGC_CSDPHASE is to decide the Ry Gates "GATEY" and the PI Gates "GATEPI" in the equivalent circuit
!	of the REAL unitary matrix X.
!
!		Decide the Ry Gates "GATEY" and the PI Gates "GATEPI"
!
!		If X is an M-by-M (M = 2**N) unitary COMPLEX matrix, we do CS decomposition recursively until the Nth level:
!
!			X = Z(1) Y(1) Z(2) Y(2) ... Z(2**N-1) Y(2**N-1) Z(2**N)
!
!		where Z(j) = diag(Z_array(j)) (j = 1,2,...,2**N), and Z(j) or Z_array(j) contains only 1 or -1.
!		In this SUBROUTINE, to make it short, the input Z(j) = Z_array(j), namely, Z(j) is a M-length array, and Z is a 2D table with the dimension
!		(M,M+N)
!		e.g. Z(j) = [1,-1,-1,1,-1,1,-1,-1]
!
!		The SUBROUTINE CYGC_CUTGATE processes this result and then:
!
!			X = Y(:,:,1) Y(:,:,2) ... Y(:,:,2**N-1) PI_MAT(:,:,1) PI_MAT(:,:,2) ... PI_MAT(:,:,N)
!
!		Meanwhile, the Ry Gates "GATEY" and the PI Gates "GATEPI" are decided.
!
!	ARGUMENT
!	========
!
!	Z: (Input) Double Precision, dimension(M,M)
!		The original matrix sequence after the decomposition.
!			X = Z(1) Y(1) Z(2) Y(2) ... Z(2**N-1) Y(2**N-1) Z(2**N)
!		or
!			X = Z(:,:,1) Y(:,:,1) Z(:,:,2) Y(:,:,2) ... Z(:,:,2**N-1) Y(:,:,2**N-1) Z(:,:,2**N)
!		In this SUBROUTINE, to make it short, the input Z(j) = Z_array(j), namely, Z(j) is a M-length array, and Z is a 2D table with the dimension
!		(M,M+N)
!	GATEPI: (Input & Output) Double Precision, dimension(M/2,N)
!		The PI Gates in the final equivalent circuit of X.
!	GATEY: (Input & Output) Double Precision, dimension(M/2,M-1)
!		The Ry Gates in the final equivalent circuit of X.
!	sgn: (Input & Output) Integer
!		If sgn = 1, then the final quantum circuit represents X
!		If sgn = -1, then the final quantum circuit represents -X
!	N: (Input) Integer
!		The number of qubits for the quantum circuit
!	M: (Input) Integer
!		The size of the unitary matrix X. M = 2**N.
!	index_level: (Input) Integer, dimension(M-1)
!		"index_level" is used to decide the recursive level of each decomposed matrix in the matrix sequence A and B.
!	index_pair: (Input) Integer, dimension(M/2,2,N)
!		"index_pair" is used to decide the structure of each decomposed matrix in the matrix sequence A and B and helps with the process afterwards.
!

implicit none
INTEGER, INTENT(IN)             :: N, M
DOUBLE PRECISION, INTENT(IN)    :: Z(M,M)
DOUBLE PRECISION, INTENT(INOUT) :: GATEPI(M/2,N)
DOUBLE PRECISION, INTENT(INOUT) :: GATEY(M/2,M-1)
INTEGER, INTENT(INOUT)          :: sgn
INTEGER, INTENT(IN)             :: index_level(M-1)
INTEGER, INTENT(IN)             :: index_pair(M/2,2,N)

INTEGER                         :: i, j, leng
INTEGER                         :: Z_array(M,M)
INTEGER                         :: PI_array(M,N)
INTEGER                         :: GATEY_sign(M/2)

!	Z(j) = diag(Z_array(j)) (j = 1,2,...,2**N), and Z(j) or Z_array(j) contains only 1 or -1.
!	In this SUBROUTINE, to make it short, the input Z(j) = Z_array(j), namely, Z(j) is a M-length array, and Z is a 2D table with the dimension
!	(M,M+N)
!	e.g. Z(j) = [1,-1,-1,1,-1,1,-1,-1]
Z_array = NINT(Z)
DO i = 2,M
	Z_array(:,i) = Z_array(:,i-1)*Z_array(:,i)
END DO

!	Decide the Ry Gates "GATEY"
GATEY_sign = 0
DO i = 1,M-1
    j = index_level(i)
    GATEY_sign = Z_array(index_pair(:,1,j),i)*Z_array(index_pair(:,2,j),i)
	GATEY(:,i) = GATEY_sign*GATEY(:,i)
END DO

!	Decide the PI Gates "GATEPI"
PI_array(:,:) = 0
PI_array(:,N) = Z_array(:,M)
DO i = 1,N-1
	leng = 2**i
	DO j = 1, 2**(N-i)
		IF (PI_array(leng*(j-1)+1,N+1-i) == -1) THEN
			PI_array(leng*(j-1)+1:leng*j,N+1-i) = -PI_array(leng*(j-1)+1:leng*j,N+1-i)
			PI_array(leng*(j-1)+1:leng*j,N-i) = -1
		ELSE
			PI_array(leng*(j-1)+1:leng*j,N-i) = 1
		END IF
	END DO
END DO
IF (PI_array(1,1) == -1) THEN
	PI_array(:,1) = -PI_array(:,1)
	sgn = -1
ELSE
	sgn = 1
END IF

DO i = 1,N
	leng = 2**(N+1-i)
	DO j = 1,2**(i-1)
		GATEPI(j,i) = PI_array(leng*j,i)
	END DO
END DO

RETURN
END
