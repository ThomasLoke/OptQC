subroutine OptQC_REAL(root,my_rank,p,args_obj)

use common_module
use csd_tools
use mpi

implicit none
character(len=128) :: finput, fhist, fperm, ftex
integer :: N, M, d, i, j, Msq
integer :: ecur, enew, esol, idx_sol
integer, allocatable :: history(:)
integer :: m_pos, delta, tol, t0, e0
! Input arrays
double precision, allocatable :: X(:,:)
! Functions
integer :: FindMinPos, CalcTol
! Workspace arrays
integer, allocatable :: index_level(:), index_pair(:,:,:)
double precision, allocatable :: COEFF(:,:)
integer, allocatable :: Perm(:), Perm_new(:), Perm_sol(:)
integer, allocatable :: QPerm(:)

! MPI variables
integer(4) :: my_rank, p, ierr, root, dummy
double precision, allocatable :: X_MPI_temp(:)
integer, allocatable :: r_col_array(:)
integer(4) :: mpi_stat(MPI_STATUS_SIZE)
! Timer variables
double precision :: c_start, c_mid, c_end, c_time1, c_time2
! Comparison variables for initial permutation
integer :: einit_root

! Objects and specification variables
integer :: nset
integer, allocatable :: type_spec(:)
type(prog_args) :: args_obj
type(csd_generator) :: csdgen_obj
type(csd_solution_set) :: csdss_Xinit, csdss_Xcur, csdss_Xnew, csdss_Xsol
type(csd_write_handle) :: csdwh_obj

! Initialize a random seed for the RNG
call init_random_seed(my_rank)
! Set the number of matrices to be 5
nset = 5

if(my_rank == root) then
    write(*,'(a)')"Performing decomposition of a real orthogonal matrix."
    write(finput,'(a,a)')args_obj%fbase(1:args_obj%flength),".txt"
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

! Start timing here so as to include time spent choosing and decomposing the initial permutation
c_start = MPI_Wtime()
! Allocate record arrays
allocate(history(args_obj%ITER_LIM+1))
allocate(r_col_array(p))
! Allocate input and output arrays
allocate(index_level(M-1))
allocate(index_pair(M/2,2,N))
allocate(COEFF(M,M))
index_level = 0
index_pair = 0
! index_level and index_pair are invariant per CSD invocation
call CYG_INDEXTABLE(N,M,index_level,index_pair)
! COEFF is also invariant per CSD invocation
call CYGC_COEFF(N,M,COEFF)
! Allocate workspace arrays
allocate(Perm(M))
allocate(Perm_new(M))
allocate(Perm_sol(M))
allocate(QPerm(N))
! Start with the identity permutation
do i = 1, M
    Perm(i) = i
    Perm_new(i) = i
    Perm_sol(i) = i
end do
! Initialize root process to identity qubit permutation - all others are randomized
if(my_rank == root) then
    do i = 1, N
        QPerm(i) = i
    end do
else
    call qperm_generate(N,QPerm)
end if
! Run object constructors
allocate(type_spec(nset))
do i = 1, nset
    ! Set all matrices to be of real type
    type_spec(i) = 0
end do
call csdgen_obj%constructor(N,M,0,index_level,index_pair,COEFF)
call csdss_Xinit%constructor(N,M,nset,type_spec)
csdss_Xinit%arr(1)%toggle_csd = .false.
csdss_Xinit%arr(5)%toggle_csd = .false.
call csdss_Xcur%constructor(N,M,nset,type_spec)
call csdss_Xnew%constructor(N,M,nset,type_spec)
call csdss_Xsol%constructor(N,M,nset,type_spec)
call csdwh_obj%constructor(N)
! Convention: U = Q^T P^T U' P Q - CHECKED AND VERIFIED
call permlisttomatrixtr(M,Perm,csdss_Xinit%arr(2)%X)                            ! P^T
call permlisttomatrix(M,Perm,csdss_Xinit%arr(4)%X)                              ! P
call qperm_process(N,M,csdss_Xinit,csdgen_obj,QPerm,X,ecur)
! Repeatedly find random permutations until the number of gates has been lowered relative to the root process (i.e. the identity permutation)
einit_root = ecur
call MPI_Bcast(einit_root,1,MPI_INTEGER8,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    j = 0
    ! Keep testing new qubit permutations until a lower number of gates is found, or until the limit PERM_ITER_LIM has been reached,
    ! in which case the identity permutation is used instead.
    do while(einit_root <= ecur)
        if(j == args_obj%PERM_ITER_LIM) then
            do i = 1, N
                QPerm(i) = i
            end do
            call qperm_process(N,M,csdss_Xinit,csdgen_obj,QPerm,X,ecur)
            exit
        end if
        call qperm_generate(N,QPerm)
        call qperm_process(N,M,csdss_Xinit,csdgen_obj,QPerm,X,ecur)
        j = j + 1
    end do
end if
c_mid = MPI_Wtime()

! Initialize solution to initial setup
call csdss_Xcur%copy(csdss_Xinit)
enew = ecur
call csdss_Xnew%copy(csdss_Xinit)
esol = ecur
call csdss_Xsol%copy(csdss_Xinit)
tol = CalcTol(args_obj%TOL_COEFF,ecur)
t0 = tol
e0 = ecur
idx_sol = 0
history(1) = ecur

do i = 1, args_obj%ITER_LIM
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
        tol = CalcTol(args_obj%TOL_COEFF,ecur)
    else
        tol = t0
    end if
end do
if(ecur < esol) then
    call csdss_Xsol%copy(csdss_Xcur)
    esol = ecur
    Perm_sol = Perm
    idx_sol = args_obj%ITER_LIM
end if
! End timing here!
c_end = MPI_Wtime()
c_time1 = c_mid - c_start
c_time2 = c_end - c_mid

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
    write(fhist,'(a,a)')args_obj%fbase(1:args_obj%flength),"_history.dat"
    write(fperm,'(a,a)')args_obj%fbase(1:args_obj%flength),"_perm.dat"
    write(ftex,'(a,a)')args_obj%fbase(1:args_obj%flength),"_circuit.tex"
    ! Assign ftex to the csd_write_handle object
    call csdwh_obj%set_file(ftex)
    ! Write the history of the algorithm to a file
    open(unit=1,file=fhist,action='write')
    do i = 0, args_obj%ITER_LIM
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
            write(2,'(f15.9,a)',advance='no')csdss_Xsol%arr(3)%X(i,j)," "
        end do
        write(2,*)
    end do
    close(2)
    ! Write the .tex files
    call csdss_Xsol%write_circuit(csdwh_obj)
    ! DEBUG CODE - OUTPUT GATES TO FILE FOR MATHEMATICA - ORDER REVERSAL FOR CORRECT CIRCUIT ORDER
    open(unit=13,file="gateseq.txt",action='write')
    do i = 1, csdss_Xsol%arr(4)%csdr_ct
        do j = 1, N+1
            write(13,'(a)')csdss_Xsol%arr(4)%Circuit(j,i)
        end do
    end do
    do i = 1, csdss_Xsol%arr(3)%csdr_ct
        do j = 1, N+1
            write(13,'(a)')csdss_Xsol%arr(3)%Circuit(j,i)
        end do
    end do
    do i = 1, csdss_Xsol%arr(2)%csdr_ct
        do j = 1, N+1
            write(13,'(a)')csdss_Xsol%arr(2)%Circuit(j,i)
        end do
    end do
    close(13)
end if

call flush(6)
call MPI_Barrier(MPI_COMM_WORLD,ierr)
dummy = 0
! Assuming root = 0
if(my_rank /= root) then
    call MPI_Recv(dummy,1,MPI_INTEGER,my_rank-1,101,MPI_COMM_WORLD,mpi_stat,ierr)
end if
write(*,'(a,i4,a,i8,a,i8)')"Process ",my_rank,": Initial number of gates before/after reduction = ",csdss_Xinit%csd_ss_ct,"/",csdss_Xinit%csdr_ss_ct
write(*,'(a,i4,a,i8,a,i8,a,i8,a,i8,a,i8,a,i8,a,i8)')"Process ",my_rank,": Breakdown of solution = ",csdss_Xsol%arr(1)%csdr_ct,"/",csdss_Xsol%arr(2)%csdr_ct,"/",csdss_Xsol%arr(3)%csdr_ct,"/",csdss_Xsol%arr(4)%csdr_ct,"/",csdss_Xsol%arr(5)%csdr_ct,"/",csdss_Xsol%csdr_ss_ct," at step number ",idx_sol
write(*,'(a,i4,a,f15.9,a,f15.9,a)')"Process ",my_rank,": Time taken = ",c_time1,"/",c_time2," seconds."
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
deallocate(COEFF)
deallocate(Perm)
deallocate(Perm_new)
deallocate(Perm_sol)
deallocate(QPerm)
deallocate(type_spec)

! Run object destructors
call csdgen_obj%destructor()
call csdss_Xinit%destructor()
call csdss_Xcur%destructor()
call csdss_Xnew%destructor()
call csdss_Xsol%destructor()
call csdwh_obj%destructor()

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
