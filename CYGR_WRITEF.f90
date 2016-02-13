SUBROUTINE CYGR_WRITEF(N, M, GATEY, GATEPI, index_level, fgate, ftex, ftexr, r_Circuit, r_col)

!
!	FUNCTION
!	========
!
!	The SUBROUTINE CYGC_CSD puts down all the information of X's quantum circuit into 2 files:
!	1. fgate: The details of all the gates in X's quantum circuit, Decimal form
!	2. ftex: The non-reduced quantum circuit of X
!   3. ftexr: The reduced quantum circuit of X
!

INTEGER, INTENT(IN)           :: N
INTEGER, INTENT(IN)           :: M
DOUBLE PRECISION, INTENT(IN)  :: GATEY(M/2,M-1)
DOUBLE PRECISION, INTENT(IN)  :: GATEPI(M/2,N)
INTEGER, INTENT(IN)           :: index_level(M-1)
character(len=128)            :: fgate, ftex, ftexr
CHARACTER(len=15)             :: r_Circuit(N+1,M*M/2+M) ! Reduced circuit from ReduceSolution call externally
integer                       :: r_col

INTEGER                       :: i, j, k, p, comma
INTEGER                       :: part, row, col, colsum, colnum
CHARACTER(len=15)             :: Circuit(N+1,M*M/2+M)   ! Non-reduced circuit, used locally
CHARACTER(len=20)             :: control, workstr
double precision              :: tol, partmax, partct

tol = 0.00000001d0
partmax = 36 / (N+1)
open(unit=1,file=fgate,status='replace')
DO i = 1,N
	WRITE(1,112)
	112 FORMAT('GATEPI')
!<ADD
	WRITE(1,1121,advance='NO') i
	1121 FORMAT(I3,';')
	IF (i .ne. 1) THEN
		DO k = 1,i-1
			IF (k .ne. i-1) THEN
				WRITE(1,11221,advance='NO') k
				11221 FORMAT(I3,',')
			ELSE
				WRITE(1,11222,advance='NO') k
				11222 FORMAT(I3)
			END IF
		END DO
	END IF
	WRITE(1,*)
!ADD>
!<ADD
	DO k = 1, 2**(i-1)
		IF (ABS(GATEPI(k,i)-1.0) <= tol) THEN
			WRITE(1,1123,advance='NO')
			1123 FORMAT('  N')
		ELSE
			WRITE(1,1124,advance='NO')
			1124 FORMAT('  Y')
		END IF
	END DO
	WRITE(1,*)
!ADD>
	WRITE(1,*)
END DO

DO j = M-1, 1, -1
	WRITE(1,113)
	113 FORMAT('GATEY')
	WRITE(1,1131,advance='NO') index_level(j)
	1131 FORMAT(I3,';')
	comma = 1
	DO k = 1,N
		IF (k .ne. index_level(j)) THEN
			IF (comma .ne. N-1) THEN
				WRITE(1,11321,advance='NO') k
				11321 FORMAT(I3,',')
			ELSE
				WRITE(1,11322,advance='NO') k
				11322 FORMAT(I3)
			END IF
			comma = comma+1
		END IF
	END DO
	WRITE(1,*)
	WRITE(1,201) GATEY(:,j)
	WRITE(1,*)
END DO
CLOSE(1)
201 FORMAT(100F8.4)

! Form the circuit for the non-reduced solution
! GATEPI
col = 0
DO i = 1,N
	DO k = 1, 2**(i-1)
		WRITE(control,301) k-1
		IF (ABS(GATEPI(k,i)+1.0) <= tol) THEN
			col = col+1
			IF (i .ne. 1) THEN
				DO row = 1,i-1
                    p = 21-i+row
					IF (control(p:p)=='0') THEN
						Circuit(row,col) = '& \ctrlo{1}'
					ELSE
						Circuit(row,col) = '& \ctrl{1}'
					END IF
				END DO
			END IF
			row = i
			Circuit(row,col) = '& \gate{\pi}'
			IF (i .ne. N) THEN
				DO row = i+1,N
					Circuit(row,col) = '& \qw'
				END DO
			END IF
			row = N+1
			Circuit(row,col) = '&'
		END IF
	END DO
END DO
!GATEY
DO j = M-1, 1, -1
	DO k = 1,M/2
		IF (ABS(GATEY(k,j)) > tol) THEN
			col = col+1
			row = index_level(j)
			Circuit(row,col) = '& \gate{R_y}'
			WRITE(control,301) k-1
			DO row = 1,N-1
                p = 21-N+row
				IF (row < index_level(j)) THEN
					IF (control(p:p)=='0') THEN
						Circuit(row,col) = '& \ctrlo{1}'
					ELSE
						Circuit(row,col) = '& \ctrl{1}'
					END IF
				ELSE
					IF (control(p:p)=='0') THEN
						Circuit(row+1,col) = '& \ctrlo{-1}'
					ELSE
						Circuit(row+1,col) = '& \ctrl{-1}'
					END IF
				END IF
			END DO
			WRITE(Circuit(N+1,col),302) GATEY(k,j)
		END IF
	END DO
END DO

! Write non-reduced solution to file
OPEN(unit=2,file=ftex,status='replace')
WRITE(2,'(a)') '\documentclass{amsart}'
WRITE(2,'(a)') '\usepackage[matrix,frame,arrow]{xypic}'
WRITE(2,'(a)') '\input{Qcircuit}'
WRITE(2,*)
WRITE(2,'(a)') '\begin{document}'
colnum = 10
colsum = col
partct = 0
DO part = 1,colsum/colnum
    if(partct == 0) then
        write(2,*)
        write(2,'(a)') '\['
        write(2,'(a)') '\Qcircuit @C=2.0em @R=0.1em @!R{'
    end if
	DO row = 1,N+1
		DO col = (part-1)*colnum+1,part*colnum
			WRITE(2,303,advance='NO') Circuit(row,col)
		END DO
		IF (row .ne. N+1) THEN
			WRITE(2,303,advance='NO') '& \qw \\'
		ELSE
			WRITE(2,303,advance='NO') '&     \\'
		END IF
		WRITE(2,*)
	END DO
	partct = partct + 1
	if(partct == partmax) then
        write(2,'(a)') '}'
        write(2,'(a)') '\]'
        partct = 0
    end if
END DO
part = colsum/colnum
if(part*colnum+1 <= colsum) then
    if(partct == 0) then
        write(2,*)
        write(2,'(a)') '\['
        write(2,'(a)') '\Qcircuit @C=2.0em @R=0.1em @!R{'
    end if
    DO row = 1,N+1
        DO col = part*colnum+1,colsum
            WRITE(2,303,advance='NO') Circuit(row,col)
        END DO
        IF (row .ne. N+1) THEN
            WRITE(2,303,advance='NO') '& \qw \\'
        ELSE
            WRITE(2,303,advance='NO') '&     \\'
        END IF
        WRITE(2,*)
    END DO
end if
WRITE(2,'(a)') '}'
WRITE(2,'(a)') '\]'
WRITE(2,*)
WRITE(2,'(a)') '\end{document}'
CLOSE(2)

open(unit=3,file=ftexr,status='replace')
WRITE(3,'(a)') '\documentclass{amsart}'
WRITE(3,'(a)') '\usepackage[matrix,frame,arrow]{xypic}'
WRITE(3,'(a)') '\input{Qcircuit}'
WRITE(3,*)
WRITE(3,'(a)') '\begin{document}'
colnum = 10
colsum = r_col
partct = 0
DO part = 1,colsum/colnum
    if(partct == 0) then
        write(3,*)
        write(3,'(a)') '\['
        write(3,'(a)') '\Qcircuit @C=2.0em @R=0.1em @!R{'
    end if
	DO row = 1,N+1
		DO col = (part-1)*colnum+1,part*colnum
			WRITE(3,303,advance='NO') r_Circuit(row,col)
		END DO
		IF (row .ne. N+1) THEN
			WRITE(3,303,advance='NO') '& \qw \\'
		ELSE
			WRITE(3,303,advance='NO') '&     \\'
		END IF
		WRITE(3,*)
	END DO
	partct = partct + 1
	if(partct == partmax) then
        write(3,'(a)') '}'
        write(3,'(a)') '\]'
        partct = 0
    end if
END DO
part = colsum/colnum
if(part*colnum+1 <= colsum) then
    if(partct == 0) then
        write(3,*)
        write(3,'(a)') '\['
        write(3,'(a)') '\Qcircuit @C=2.0em @R=0.1em @!R{'
    end if
    DO row = 1,N+1
        DO col = part*colnum+1,colsum
            WRITE(3,303,advance='NO') r_Circuit(row,col)
        END DO
        IF (row .ne. N+1) THEN
            WRITE(3,303,advance='NO') '& \qw \\'
        ELSE
            WRITE(3,303,advance='NO') '&     \\'
        END IF
        WRITE(3,*)
    END DO
end if
WRITE(3,'(a)') '}'
WRITE(3,'(a)') '\]'
WRITE(3,*)
WRITE(3,'(a)') '\end{document}'
CLOSE(3)
301 FORMAT(B20.20)
302 FORMAT('& ',F8.4)
303 FORMAT(A15)

RETURN
END
