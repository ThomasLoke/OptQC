subroutine OptQC_REAL(my_rank,p,fbase,flength,ITER_LIM)

use memwork_real
use mpi

implicit none
character(len=128) :: fbase, finput, fhist, fperm, fgate, ftex, ftexr, fin, fmat
integer :: flength
integer :: N, M, d, i, j, Msq
integer :: size0, size1, num0, num1
integer :: ITER_LIM
integer :: ecur, enew
integer :: nP, nPT, ntot
integer, allocatable :: history(:), history_Delta(:)
character(len=15), allocatable :: r_Circuit(:,:)
integer :: r_col, m_pos
integer :: delta, tol, t0, e0
! Input arrays
double precision, allocatable :: X(:,:)
! Functions
integer :: DecomposeAndCount, DecomposeAndCountReduced, DecomposeAndCountReducedWP, FindMinPos, CalcTol
! Workspace arrays
! Local:
integer, allocatable :: Perm(:), Perm_new(:), Perm_sol(:)
double precision, allocatable :: Xcur(:,:), Xnew(:,:)
double precision, allocatable :: Pmat(:,:), PTmat(:,:)

! Storage of best solution
integer :: esol, idx_sol
double precision, allocatable :: X_sol(:,:)

! MPI variables
integer(4) :: my_rank, p, ierr, root, dummy
double precision, allocatable :: X_MPI_temp(:)
integer, allocatable :: r_col_array(:)
integer(4) :: mpi_stat(MPI_STATUS_SIZE)
! Timer variables
double precision :: c_start, c_end, e_time

! Choose root process as 0
root = 0
! Initialize a random seed for the RNG
call init_random_seed()

if(my_rank == root) then
    write(*,*)
    write(*,'(a)')"Performing decomposition of a real orthogonal matrix."
    write(finput,'(a,a)')fbase(1:flength),".txt"
    write(*,'(a,a)')"Retrieving unitary matrix from file ",finput
    open(unit=1,file=finput,action='read')
    read(1,*)d
    N = ceiling(log(real(d))/log(2.0))
    M = 2**N
    allocate(X(M,M))
    allocate(X_MPI_temp(M*M))
    X = 0.0d0
    do i = 1, d
    	do j = 1, d
    		read(1,*)X(i,j)
    	end do
    end do
    if(M /= d) then
    	do i = d+1, M
    		X(i,i) = 1
    	end do
    end if
    close(1)
    Msq = M*M
end if

! Broadcast the matrix X from the root process to all the other processes
! Uses a flattened matrix for broadcasting, which is reshaped into a matrix.
call MPI_Bcast(N,1,MPI_INTEGER8,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    M = 2**N
    Msq = M*M
    allocate(X(M,M))            ! Currently not in use....
    allocate(X_MPI_temp(Msq))
end if

allocate(history(ITER_LIM+1))
allocate(history_Delta(ITER_LIM))
allocate(r_col_array(p))
! Allocate input and output arrays
allocate(index_level(M-1))
allocate(index_pair(M/2,2,N))
index_level = 0
index_pair = 0
! index_level and index_pair are invariant per CSD invocation
call CYG_INDEXTABLE(N,M,index_level,index_pair)

! Allocate workspace arrays
allocate(Perm(M))
allocate(Perm_new(M))
allocate(Perm_sol(M))
allocate(Xcur(M,M))
allocate(Xnew(M,M))
allocate(Pmat(M,M))
allocate(PTmat(M,M))
allocate(X_sol(M,M))

if(my_rank == root) then
    ! Find a 'well-arranged' permutation to start from
    ! Initialize the matrix of Hamming distances between all i and j
    esol = DecomposeAndCount(N,M,X,GATEY,GATEPI)
    write(*,'(a,i8)')"Original number of gates (before reduction) = ",esol
    esol = DecomposeAndCountReduced(N,M,X,GATEY,GATEPI,index_level,r_Circuit)
    write(*,'(a,i8)')"Original number of gates (after reduction) = ",esol
    Xcur = X
    ecur = esol
    do i = 1, M
        Perm(i) = i
    end do
    call Flatten(M,Xcur,X_MPI_temp)
end if

call MPI_Bcast(ecur,1,MPI_INTEGER8,root,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Perm,M,MPI_INTEGER8,root,MPI_COMM_WORLD,ierr)
call MPI_Bcast(X_MPI_temp,Msq,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    call Box(M,X_MPI_temp,Xcur)
end if

esol = ecur
idx_sol = 0
X_sol = Xcur
Perm_sol = Perm
history(1) = ecur
Xnew = Xcur
enew = ecur
Perm_new = Perm
tol = CalcTol(ecur)
t0 = tol
e0 = ecur

c_start = MPI_Wtime()
do i = 1, ITER_LIM
    Perm_new = Perm
    call NeighbourhoodOpt(N,M,Xcur,Xnew,Perm_new)
    enew = DecomposeAndCountReducedWP(N,M,Xnew,Perm_new,Pmat,PTmat,GATEY,GATEPI,index_level,r_Circuit)
    delta = enew - ecur
    !write(*,*)"New number of gates = ",ecur,"/",enew
    if(enew <= ecur) then
        Xcur = Xnew
        ecur = enew
        Perm = Perm_new
        if(enew < esol) then
            esol = ecur
            idx_sol = i
            X_sol = Xcur
            Perm_sol = Perm
        end if
    else
        ! No improvement to the cost function
        if(delta <= tol) then
            Xcur = Xnew
            ecur = enew
            Perm = Perm_new
        end if
    end if
    history(i+1) = ecur
    history_Delta(i) = delta
    ! Limit the maximum tolerance to the initial tolerance value
    if(ecur < e0) then
        tol = CalcTol(ecur)
    else
        tol = t0
    end if
end do
if(ecur < esol) then
    esol = ecur
    idx_sol = ITER_LIM
    X_sol = Xcur
    Perm_sol = Perm
end if
r_col = DecomposeAndCountReduced(N,M,X_sol,GATEY,GATEPI,index_level,r_Circuit)
c_end = MPI_Wtime()
e_time = c_end - c_start

! Write the .tex file based on the best solution
call MPI_Allgather(esol,1,MPI_INTEGER8,r_col_array,1,MPI_INTEGER8,MPI_COMM_WORLD,ierr)
m_pos = FindMinPos(int(p,8),r_col_array)    ! Typecast to integer(8) since p is of type integer(4)
! Perform I/O operations from the process with the least number of gates - eliminates need
! to communicate results back to the root.
if(my_rank == m_pos-1) then
    write(*,'(a,i4,a,i8,a)')"Process ",my_rank," has the least number of gates: ",esol," gates."
    write(*,'(a)')"Writing solution to files."
    write(*,*)
    ! Initialize filenames
    write(fhist,'(a,a)')fbase(1:flength),"_history.dat"
    write(fperm,'(a,a)')fbase(1:flength),"_perm.dat"
    write(fgate,'(a,a)')fbase(1:flength),"_gates.txt"
    write(ftex,'(a,a)')fbase(1:flength),"_plot.tex"
    write(ftexr,'(a,a)')fbase(1:flength),"_plot_reduced.tex"
    ! Write the .tex files
    call CYGR_WRITEF(N,M,GATEY,GATEPI,index_level,fgate,ftex,ftexr,r_Circuit,r_col)
    ! Write the history of the algorithm to a file
    open(unit=1,file=fhist,action='write')
    do i = 0, ITER_LIM
        write(1,'(i8,a,i8)')i," ",history(i+1)
    end do
    close(1)
    ! Write the optimal permutation and the corresponding matrix to a file
    open(unit=2,file=fperm,action='write')
    do i = 1, M
        write(2,'(i15,a)',advance='no')Perm_sol(i)," "
    end do
    write(2,*)
    do i = 1, M
        do j = 1, M
            write(2,'(f15.9,a)',advance='no')X_sol(i,j)," "
        end do
        write(2,*)
    end do
    close(2)
end if

call permlisttomatrix(M,Perm_sol,Pmat)
PTmat = transpose(Pmat)
nP = DecomposeAndCountReduced(N,M,Pmat,GATEY,GATEPI,index_level,r_Circuit)
nPT = DecomposeAndCountReduced(N,M,PTmat,GATEY,GATEPI,index_level,r_Circuit)
ntot = r_col + nP + nPT

call flush(6)
call MPI_Barrier(MPI_COMM_WORLD,ierr)
dummy = 0
! Assuming root = 0
if(my_rank /= root) then
    call MPI_Recv(dummy,1,MPI_INTEGER,my_rank-1,101,MPI_COMM_WORLD,mpi_stat,ierr)
else
    write(*,'(a,i4)')"Size of communicator: ",p
end if
write(*,'(a,i4,a,i8,a,i8)')"Process ",my_rank,": Final number of gates = ",esol," at step number ",idx_sol
write(*,'(a,i4,a,i8,a,i8,a,i8,a,i8,a)')"Process ",my_rank,": (r_col,nP,nPT,ntot) = (",r_col,",",nP,",",nPT,",",ntot,")"
write(*,'(a,i4,a,f15.9,a)')"Process ",my_rank,": Time taken = ",e_time," seconds."
call flush(6)
if(my_rank /= p-1) then
    call MPI_Send(dummy,1,MPI_INTEGER,my_rank+1,101,MPI_COMM_WORLD,ierr)
end if

deallocate(history)
deallocate(history_Delta)
deallocate(r_col_array)
deallocate(X,index_level,index_pair)
deallocate(Perm)
deallocate(Perm_new)
deallocate(Perm_sol)
deallocate(Xcur)
deallocate(Xnew)
deallocate(Pmat)
deallocate(PTmat)
deallocate(X_sol)
deallocate(X_MPI_temp)

end subroutine OptQC_REAL

subroutine Box(M,source,dest)

implicit none
integer :: M
double precision :: source(M*M), dest(M,M)

integer :: i, j, ct

ct = 1
do i = 1, M
    do j = 1, M
        dest(i,j) = source(ct)
        ct = ct + 1
    end do
end do
return

end subroutine Box

subroutine Flatten(M,source,dest)

implicit none
integer :: M
double precision :: source(M,M), dest(M*M)

integer :: i, j, ct

ct = 1
do i = 1, M
    do j = 1, M
        dest(ct) = source(i,j)
        ct = ct + 1
    end do
end do
return

end subroutine Flatten
