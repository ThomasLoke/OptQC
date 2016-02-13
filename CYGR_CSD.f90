SUBROUTINE CYGR_CSD(N,M,X)

use memwork_real

!
!	FUNCTION
!	========
!
!	For an arbitray REAL unitary matrix X, the SUBROUTINE CYGR_CSD applies Cosine Sine Decomposition (CSD) method to transform X
!	into its equivalent quantum circuit.
!	More specifically,
!	1. The SUBROUTINE CYGR_CSD reads in the unitary matrix X for the input file and makes X into a 2**N-by-2**N unitary matrix.
!	2. CYGR_CSD creates two useful array/table "index_level" & "index_pair", and does Cosine Sine Decomposition (CSD) of X recursively
!	   until the Nth level
!	3. CYGR_CSD processes the CSD results
!	   In this step, we obtain all the Ry Gates and the PI Gates for X: GATEY, GATEPI
!	4. CYGR_CSD writes down the complete CSD result and the information of every gate in the quantum circuit.
!

IMPLICIT NONE
INTEGER                       :: d_sum, N, M
INTEGER                       :: i, j
INTEGER                       :: size0, size1, num0, num1
DOUBLE PRECISION              :: X(M,M)
INTEGER                       :: sgn
INTEGER                       :: index_Y

GATEY = 0.0d0
GATEPI = 0.0d0

!	Cosine Sine Decomposition (CSD) of X recursively until the Nth level
!	The 1st level
size0 = M
size1 = M/2
num0 = 1
num1 = 2
Z1%l(1)%arr = 0.0d0
index_Y = 2**(N-1)
CALL CYGR_BLKCSD(X,Z1%l(1)%arr(:,:,:,1),Z1%l(1)%arr(:,:,:,2),GATEY(:,index_Y),N,M,size0,size1,num0,num1,X_blk%l(1)%arr,X11%l(1)%arr,X12%l(1)%arr,X21%l(1)%arr,X22%l(1)%arr,U1%l(1)%arr,U2%l(1)%arr,V1T%l(1)%arr,V2T%l(1)%arr,GATEY_blk%l(1)%arr,IWORK%l(1)%arr)
Z0%l(1)%arr = Z1%l(1)%arr

!	The ist level (i = 2,...,N)
DO i = 2, N
    size1 = 2**(N-i)
    size0 = size1 * 2
    num0 = 2**(i-1)
    num1 = num0 * 2
    Z1%l(i)%arr = 0.0d0
    DO j = 1, num0
        index_Y = 2**(N-i)*(2*j-1)
        CALL CYGR_BLKCSD(Z0%l(i-1)%arr(:,:,:,j),Z1%l(i)%arr(:,:,:,2*j-1),Z1%l(i)%arr(:,:,:,2*j),GATEY(:,index_Y),N,M,size0,size1,num0,num1,X_blk%l(i)%arr,X11%l(i)%arr,X12%l(i)%arr,X21%l(i)%arr,X22%l(i)%arr,U1%l(i)%arr,U2%l(i)%arr,V1T%l(i)%arr,V2T%l(i)%arr,GATEY_blk%l(i)%arr,IWORK%l(i)%arr)
    END DO
    Z0%l(i)%arr = Z1%l(i)%arr
END DO

!	Step 3: CYGC_CSD processes the CSD results
!	        In this step, we obtain all the Ry Gates and the PI Gates for X: GATEY, GATEPI
CALL CYGR_CUTGATE(Z0%l(N)%arr(1,1,:,:), GATEPI, GATEY, sgn, N, M, index_level, index_pair, Z_array, PI_array, GATEY_sign)

RETURN
END
