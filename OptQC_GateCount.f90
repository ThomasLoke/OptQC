function GateCount(N, M, GATEY, GATEPI)

implicit none
integer :: N
integer :: M
double precision :: GATEY(M/2,M-1)
double precision :: GATEPI(M/2,N)

integer :: i, j, GateCount
double precision :: tol

tol = 0.00000001d0
GateCount = 0
do i = 1, N
    do j = 1, 2**(i-1)
		if (abs(GATEPI(j,i)+1.0) <= tol) then
			GateCount = GateCount + 1
		end if
	end do
end do

do i = M-1, 1, -1
    do j = M/2, 1, -1
		if (abs(GATEY(j,i)) > 10**(-8)) then
			GateCount = GateCount + 1
		end if
	end do
end do

return

end function GateCount

function DecomposeAndCount(N, M, X, GATEY, GATEPI)

implicit none
integer :: N, M
double precision :: X(M,M)
double precision :: GATEY(M/2,M-1)
double precision :: GATEPI(M/2,N)
integer :: GateCount

integer :: DecomposeAndCount

call CYGR_CSD(N,M,X)
DecomposeAndCount = GateCount(N,M,GATEY,GATEPI)
return

end function DecomposeAndCount

function DecomposeAndCountReduced(N, M, X, GATEY, GATEPI, index_level, r_Circuit)

implicit none
integer :: N, M
double precision :: X(M,M)
double precision :: GATEY(M/2,M-1)
double precision :: GATEPI(M/2,N)
integer :: index_level(M-1)
character(len=15) :: r_Circuit(N+1,M*M/2+M)

integer :: DecomposeAndCountReduced

call CYGR_CSD(N,M,X)
call ReduceSolution(N,M,GATEY,GATEPI,index_level,r_Circuit,DecomposeAndCountReduced)
return

end function DecomposeAndCountReduced

function DecomposeAndCountReducedWP(N, M, X, Perm, Pmat, PTmat, GATEY, GATEPI, index_level, r_Circuit)

implicit none
integer :: N, M
double precision :: X(M,M)
integer :: Perm(M)
double precision :: Pmat(M,M), PTmat(M,M)
double precision :: GATEY(M/2,M-1)
double precision :: GATEPI(M/2,N)
integer :: index_level(M-1)
character(len=15) :: r_Circuit(N+1,M*M/2+M)
integer :: DecomposeAndCountReduced

integer :: DecomposeAndCountReducedWP
integer :: nX, nP, nPT

call permlisttomatrix(M,Perm,Pmat)
PTmat = transpose(Pmat)
nP = DecomposeAndCountReduced(N,M,Pmat,GATEY,GATEPI,index_level,r_Circuit)
nPT = DecomposeAndCountReduced(N,M,PTmat,GATEY,GATEPI,index_level,r_Circuit)
nX = DecomposeAndCountReduced(N,M,X,GATEY,GATEPI,index_level,r_Circuit)
DecomposeAndCountReducedWP = nX + nP + nPT
return

end function DecomposeAndCountReducedWP
