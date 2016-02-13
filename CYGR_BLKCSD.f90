SUBROUTINE CYGR_BLKCSD (X,U,V,GATEY,N,M,size0,size1,num0,num1,X_blk,X11,X12,X21,X22,U1,U2,V1T,V2T,GATEY_blk,IWORK)

!
!	FUNCTION
!	========
!
!	DEFINITION
!
!	CS decomposition of an M-by-M (M = 2**N) unitary matrix X is:
!
!	     [ U1 |    ] [  C | S ] [ V1 |    ]**H
!	X =  [---------] [--------] [---------]
!	     [    | U2 ] [ -S | C ] [    | V2 ]
!	  =  U Y V
!
!	where theta is a M/2 array, C = diag(cos(theta)), S = diag(sin(theta))
!
!	FUNCTION
!
!	Since U1, U2, V1, V2 are still unitary matrices, CS decomposition can be applied recursively. In the ith level, the input unitary matrix X is:
!
!	    [ X(1)                       ]
!	    [      X(2)                  ]
!	X = [           X(3)             ]
!	    [                ...         ]
!	    [                    X(num0) ]
!
!	where num0 = 2**(N-1), and X(1), X(2), ... X(num0) are P0-by-P0 unitary matrices, P0 = 2**(N-i+1)
!	The SUBROUTINE CYGC_BLKCSD calculate the CS decomposition of X in the ith level.
!	First, for X(1), we can do the CS decomposition:
!
!	        [ U1 |    ] [  C | S ] [ V1 |    ]**H
!	X(1) =  [---------] [--------] [---------]
!	        [    | U2 ] [ -S | C ] [    | V2 ]
!	     =  U(1) Y(1) V(1)
!
!	where C = diag(cos(GATEY(1))), S = diag(sin(GATEY(1)))
!
!	Then, we repeat the previous step for X(k) (k = 1,2,...,num0):
!
!	        [ U1 |    ] [  C | S ] [ V1 |    ]**H
!	X(k) =  [---------] [--------] [---------]
!	        [    | U2 ] [ -S | C ] [    | V2 ]
!	     =  U(k) Y(k) V(k)
!	(k = 1,2,...,num0)
!
!	Finally the result can be summerized as:
!
!	    [ X(1)                       ]   [ U(1)                       ] [ Y(1)                       ] [ V(1)                       ]**H
!	    [      X(2)                  ]   [      U(2)                  ] [      Y(2)                  ] [      V(2)                  ]
!	X = [           X(3)             ] = [           U(3)             ] [           Y(3)             ] [           V(3)             ]
!	    [                ...         ]   [                ...         ] [                ...         ] [                ...         ]
!	    [                    X(num0) ]   [                    U(num0) ] [                    Y(num0) ] [                    V(num0) ]
!
!	                                   = U Y V
!
!	GATEY = [GATEY(1), GATEY(2), ..., GATEY(num0)]
!	Note that finally GATEY is a M/2 array.
!
!	ARGUMENT
!	========
!
!	N: (Input) Integer
!		The number of qubits for the quantum circuit
!	M: (Input) Integer
!		The size of the unitary matrix X. M = 2**N.
!	X: (Input) Double Precision, dimension(M,M)
!		The block diagonal unitary matrix that is to be decomposed.
!		    [ X(1)                       ]
!		    [      X(2)                  ]
!		X = [           X(3)             ]
!		    [                ...         ]
!		    [                    X(num0) ]
!	U: (Input & Output) Double Precision, dimension(M,M)
!		The block diagonal unitary matrix we obtain after the decomposition.
!		    [ U(1)                       ]
!		    [      U(2)                  ]
!		U = [           U(3)             ]
!		    [                ...         ]
!		    [                    U(num0) ]
!	V: (Input & Output) Double Precision, dimension(M,M)
!		The block diagonal unitary matrix we obtain after the decomposition.
!		    [ V(1)                       ]**H
!		    [      V(2)                  ]
!		V = [           V(3)             ]
!		    [                ...         ]
!		    [                    V(num0) ]
!	GATEY: (Input & Output) Double Precision, dimension(M/2)
!		The Ry Gate we obtain after the decomposition. But for the REAL X, this Ry Gate is not the one in the final equivalent circuit of X.
!	size0: (Input) Integer
!		If i is the recursive level of the CSD, then
!		size0 = 2**(N-i+1) is the size of the matrix blocks X(1), X(2), ..., X(num0)
!	size1: (Input) Integer
!		If i is the recursive level of the CSD, then
!		size1 = 2**(N-i)
!	num0: (Input) Integer
!		If i is the recursive level of the CSD, then
!		num0 = 2**(i-1) is the number of the matrix blocks X(1), X(2), ..., X(num0)
!	num1: (Input) Integer
!		If i is the recursive level of the CSD, then
!		num1 = 2**i
!

IMPLICIT NONE
INTEGER, INTENT(IN)             :: N, M
INTEGER, INTENT(IN)             :: size0, size1, num0, num1
DOUBLE PRECISION, INTENT(IN)    :: X(size0,size0,num0)
DOUBLE PRECISION, INTENT(INOUT) :: U(size1,size1,num1), V(size1,size1,num1)
DOUBLE PRECISION, INTENT(INOUT) :: GATEY(M/2)

integer                         :: k, k1, k2
DOUBLE PRECISION                :: X_blk(size0,size0)
DOUBLE PRECISION                :: GATEY_blk(size1)
DOUBLE PRECISION                :: X11(size1,size1), X12(size1,size1), X21(size1,size1), X22(size1,size1), U1(size1,size1), U2(size1,size1), V1T(size1,size1), V2T(size1,size1)
DOUBLE PRECISION, ALLOCATABLE   :: WORK(:)
INTEGER                         :: IWORK(size1)
INTEGER                         :: LWORK, INFO
DOUBLE PRECISION                :: Pi
DOUBLE PRECISION                :: threshold

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

threshold = 1.0E-20
DO k = 1, num0

	!	    [ X(1)                       ]
	!	    [      X(2)                  ]
	!	X = [           X(3)             ]
	!	    [                ...         ]
	!	    [                    X(num0) ]
	!
	!	For X(k) (k = 1,...,num0), we can do the CS decomposition:
	!
	!	X_blk = X(k) = X(P0*(k-1)+1:P0*k,P0*(k-1)+1:P0*k)
	!
	!	         [ U1 |    ] [  C | S ] [ V1 |    ]**H
	!	X_blk =  [---------] [--------] [---------]
	!	         [    | U2 ] [ -S | C ] [    | V2 ]
	!	      =  U_blk Y_blk V_blk
	!
	!	where C = diag(cos(GATEY_blk)), S = diag(sin(GATEY_blk))

	X_blk = X(:,:,k)
	DO k1 = 1,size0
		DO k2 = 1,size0
			IF(ABS(X_blk(k1,k2)) <= threshold) THEN
				X_blk(k1,k2) = 0.0
			END IF
		END DO
	END DO

	!
	!	CSD
    !	===
    !
    !	The following part: Applying Sutton's CSD package (DORCSD)
	X11 = X_blk(1:size1,1:size1)
	X12 = X_blk(1:size1,size1+1:size0)
	X21 = X_blk(size1+1:size0,1:size1)
	X22 = X_blk(size1+1:size0,size1+1:size0)

	LWORK = -1
	ALLOCATE(WORK(1))
	WORK = 0
		!	Set LWORK = -1, then a workspace query is assumed; the subroutine only calculates the optimal size of the WORK array, returns this value as the first entry of the work array.
	CALL DORCSD('Y', 'Y', 'Y', 'Y', 'N', 'O', size0, size1, size1, X11,&
size1, X12, size1, X21, size1, X22, size1, GATEY_blk, U1, size1, U2, &
size1, V1T, size1, V2T, size1, WORK, LWORK, IWORK, INFO)
		!	Now set LWORK = WORK(1)
	LWORK = WORK(1)
	DEALLOCATE(WORK)
	ALLOCATE(WORK(LWORK))
	WORK = 0
	CALL DORCSD('Y', 'Y', 'Y', 'Y', 'N', 'O', size0, size1, size1, X11, &
size1, X12, size1, X21, size1, X22, size1, GATEY_blk, U1, size1, U2, &
size1, V1T, size1, V2T, size1, WORK, LWORK, IWORK, INFO)
	DEALLOCATE(WORK)
	!	Applying Sutton's CSD package (DORCSD): The end

    !	U(P0*(k-1)+1:P0*k,P0*(k-1)+1:P0*k) = U(k) = U_blk
	!	V(P0*(k-1)+1:P0*k,P0*(k-1)+1:P0*k) = V(k) = V_blk
	!	GATEY(P1*(k-1)+1:P1*k) = GATEY(k) = GATEY_blk
    !	Finally, we have
    !
	!	    [ X(1)                       ]   [ U(1)                       ] [ Y(1)                       ] [ V(1)                       ]**H
	!	    [      X(2)                  ]   [      U(2)                  ] [      Y(2)                  ] [      V(2)                  ]
	!	X = [           X(3)             ] = [           U(3)             ] [           Y(3)             ] [           V(3)             ]
	!	    [                ...         ]   [                ...         ] [                ...         ] [                ...         ]
	!	    [                    X(num0) ]   [                    U(num0) ] [                    Y(num0) ] [                    V(num0) ]
	!
	!	                                   = U Y V
    !
    !	GATEY = [GATEY(1), GATEY(2), ..., GATEY(num0)]
    !	Note that we skip the construction of the matrix Y because we only care about the Ry Gate rather than the matrix Y.
	U(:,:,2*k-1) = U1
	U(:,:,2*k) = U2
	V(:,:,2*k-1) = V1T
	V(:,:,2*k) = V2T
	GATEY(size1*(k-1)+1:size1*k) = GATEY_blk
END DO

Pi = 3.141592653589793d0
GATEY = GATEY/Pi

RETURN
END
