subroutine OptQC_REAL(my_rank,p,fbase,flength,ITER_LIM)

use csd_real
use mpi

implicit none
character(len=128) :: fbase, finput, fhist, fperm, ftex, fin, fmat
integer :: flength
integer :: N, M, d, i, j, Msq
integer :: ITER_LIM
integer :: ecur, enew, esol, idx_sol
integer, allocatable :: history(:)
integer :: m_pos, delta, tol, t0, e0
! Input arrays
double precision, allocatable :: X(:,:)
! Functions
integer :: FindMinPos, CalcTol
! Workspace arrays
integer, allocatable :: index_level(:), index_pair(:,:,:)
integer, allocatable :: Perm(:), Perm_new(:), Perm_sol(:)

! MPI variables
integer(4) :: my_rank, p, ierr, root, dummy
double precision, allocatable :: X_MPI_temp(:)
integer, allocatable :: r_col_array(:)
integer(4) :: mpi_stat(MPI_STATUS_SIZE)
! Timer variables
double precision :: c_start, c_end, e_time

! Objects
type(csd_generator) :: csdgen_obj
type(csd_solution_set) :: csdss_Xinit, csdss_Xcur, csdss_Xnew, csdss_Xsol

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
    Msq = M*M
    allocate(X(M,M))
    allocate(X_MPI_temp(Msq))
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
    call Flatten(M,X,X_MPI_temp)
end if

! Broadcast the matrix X from the root process to all the other processes
! Uses a flattened matrix for broadcasting, which is reshaped into a matrix.
call MPI_Bcast(N,1,MPI_INTEGER8,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    M = 2**N
    Msq = M*M
    allocate(X(M,M))
    allocate(X_MPI_temp(Msq))
end if
call MPI_Bcast(X_MPI_temp,Msq,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    call Box(M,X_MPI_temp,X)
end if

! Allocate record arrays
allocate(history(ITER_LIM+1))
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
! Start with the identity permutation
do i = 1, M
    Perm(i) = i
    Perm_new(i) = i
    Perm_sol(i) = i
end do
! Run object constructors
call csdgen_obj%constructor(N,M,index_level,index_pair)
call csdss_Xinit%constructor(3,N,M)
call csdss_Xcur%constructor(3,N,M)
call csdss_Xnew%constructor(3,N,M)
call csdss_Xsol%constructor(3,N,M)
! Convention: U = P^T U' P - PLEASE CHECK THIS
call permlisttomatrixtr(M,Perm,csdss_Xinit%arr(1)%X)
csdss_Xinit%arr(2)%X = X
call permlisttomatrix(M,Perm,csdss_Xinit%arr(3)%X)
! Count initial number of gates (including reduction)
call csdss_Xinit%arr(2)%run_csdr(csdgen_obj)
csdss_Xinit%csd_ss_ct = csdss_Xinit%arr(2)%csd_ct
csdss_Xinit%csdr_ss_ct = csdss_Xinit%arr(2)%csdr_ct
ecur = csdss_Xinit%csdr_ss_ct

!goto 3010

! Initialize solution to initial setup
call csdss_Xcur%copy(csdss_Xinit)
enew = ecur
call csdss_Xnew%copy(csdss_Xinit)
esol = ecur
call csdss_Xsol%copy(csdss_Xinit)
tol = CalcTol(ecur)
t0 = tol
e0 = ecur
idx_sol = 0
history(1) = ecur

c_start = MPI_Wtime()
do i = 1, ITER_LIM
    Perm_new = Perm
    call NeighbourhoodOpt(M,csdss_Xcur,csdss_Xnew,Perm_new)
    call csdss_Xnew%run_csdr(csdgen_obj)
    enew = csdss_Xnew%csdr_ss_ct
    delta = enew - ecur
    !write(*,*)"New number of gates = ",ecur,"/",enew
    if(enew <= ecur) then
        call csdss_Xcur%copy(csdss_Xnew)
        ecur = enew
        Perm = Perm_new
        if(enew < esol) then
            call csdss_Xsol%copy(csdss_Xcur)
            esol = ecur
            Perm_sol = Perm
            idx_sol = i
        end if
    else
        ! No improvement to the cost function
        if(delta <= tol) then
            call csdss_Xcur%copy(csdss_Xnew)
            ecur = enew
            Perm = Perm_new
        end if
    end if
    history(i+1) = ecur
    ! Limit the maximum tolerance to the initial tolerance value
    if(ecur < e0) then
        tol = CalcTol(ecur)
    else
        tol = t0
    end if
end do
if(ecur < esol) then
    call csdss_Xsol%copy(csdss_Xcur)
    esol = ecur
    Perm_sol = Perm
    idx_sol = ITER_LIM
end if
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
    write(ftex,'(a,a)')fbase(1:flength),"_circuit.tex"
    ! Write the history of the algorithm to a file
    open(unit=1,file=fhist,action='write')
    do i = 0, ITER_LIM
        write(1,'(i8,a,i8)')i," ",history(i+1)
    end do
    close(1)
    ! Write the optimal permutation and the corresponding matrix to a file - PLEASE KEEP IN VIEW
    open(unit=2,file=fperm,action='write')
    do i = 1, M
        write(2,'(i15,a)',advance='no')Perm_sol(i)," "
    end do
    write(2,*)
    do i = 1, M
        do j = 1, M
            write(2,'(f15.9,a)',advance='no')csdss_Xsol%arr(2)%X(i,j)," "
        end do
        write(2,*)
    end do
    close(2)
    ! Write the .tex files
    call csdss_Xsol%write_circuit(ftex)
end if

call flush(6)
call MPI_Barrier(MPI_COMM_WORLD,ierr)
dummy = 0
! Assuming root = 0
if(my_rank /= root) then
    call MPI_Recv(dummy,1,MPI_INTEGER,my_rank-1,101,MPI_COMM_WORLD,mpi_stat,ierr)
else
    write(*,'(a,i4)')"Size of communicator: ",p
end if
write(*,'(a,i4,a,i8,a,i8)')"Process ",my_rank,": Initial number of gates before/after reduction = ",csdss_Xinit%csd_ss_ct,"/",csdss_Xinit%csdr_ss_ct
write(*,'(a,i4,a,i8,a,i8,a,i8,a,i8,a,i8)')"Process ",my_rank,": Breakdown of solution = ",csdss_Xsol%arr(1)%csdr_ct,"/",csdss_Xsol%arr(2)%csdr_ct,"/",csdss_Xsol%arr(3)%csdr_ct,"/",csdss_Xsol%csdr_ss_ct," at step number ",idx_sol
write(*,'(a,i4,a,f15.9,a)')"Process ",my_rank,": Time taken = ",e_time," seconds."
call flush(6)
if(my_rank /= p-1) then
    call MPI_Send(dummy,1,MPI_INTEGER,my_rank+1,101,MPI_COMM_WORLD,ierr)
end if

! Free up memory
deallocate(X)
deallocate(X_MPI_temp)
deallocate(history)
deallocate(r_col_array)
deallocate(index_level)
deallocate(index_pair)
deallocate(Perm)
deallocate(Perm_new)
deallocate(Perm_sol)

! Run object destructors
call csdgen_obj%destructor()
call csdss_Xinit%destructor()
call csdss_Xcur%destructor()
call csdss_Xnew%destructor()
call csdss_Xsol%destructor()

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
