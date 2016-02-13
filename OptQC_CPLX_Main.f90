subroutine OptQC_CPLX(root,my_rank,p,args_obj)

use common_module
use csd_tools
use mpi

implicit none
character(len=128) :: finput, fhist, fperm, ftex
double precision :: in_1, in_2
integer :: N, M, d, i, j, Msq
integer :: ecur, enew, esol, idx_sol
integer, allocatable :: history(:)
integer :: hist_idx, m_pos, delta, tol, t0, e0
! Input arrays
double complex, allocatable :: X(:,:)
! Functions
integer :: FindMinPos, CalcTol, ChooseN
! Workspace arrays
integer, allocatable :: index_level(:), index_pair(:,:,:)
double precision, allocatable :: COEFF(:,:)
integer, allocatable :: Perm(:), Perm_new(:), Perm_sol(:)
integer, allocatable :: QPerm(:), QPerm_new(:)

! MPI variables
integer :: my_rank, p, ierr, root
integer, allocatable :: MPI_int_buffer(:)
double complex, allocatable :: X_MPI_temp(:)
integer, allocatable :: r_col_array(:)
integer :: mpi_stat(MPI_STATUS_SIZE)
! Timer variables
double precision :: c_start, c_mid, c_end, c_time1, c_time2

! Objects
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
! Allocate integer buffer
allocate(MPI_int_buffer(1))

if(my_rank == root) then
    write(*,'(a)')"Performing decomposition of a complex unitary matrix."
    write(finput,'(a,a)')args_obj%fbase(1:args_obj%flength),".txt"
    write(*,'(a,a)')"Retrieving unitary matrix from file ",finput
    open(unit=1,file=finput,action='read')
    read(1,*)d
    N = ChooseN(d)
    M = 2**N
    Msq = M*M
    allocate(X(M,M))
    allocate(X_MPI_temp(Msq))
    X = cmplx(0.0d0)
    do i = 1, d
    	do j = 1, d
    		read(1,*)in_1,in_2
    		X(i,j) = cmplx(in_1,in_2)
    	end do
    end do
    if(M /= d) then
    	do i = d+1, M
    		X(i,i) = cmplx(1.0d0,0.0d0)
    	end do
    end if
    close(1)
    MPI_int_buffer(1) = N
    call Flatten_CPLX(M,X,X_MPI_temp)
end if

! Broadcast the matrix X from the root process to all the other processes
! Uses a flattened matrix for broadcasting, which is reshaped into a matrix.
call MPI_Bcast(MPI_int_buffer,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    N = MPI_int_buffer(1)
    M = 2**N
    Msq = M*M
    allocate(X(M,M))
    allocate(X_MPI_temp(Msq))
end if
call MPI_Bcast(X_MPI_temp,Msq,MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    call Box_CPLX(M,X_MPI_temp,X)
end if

! Start timing here so as to include time spent choosing and decomposing the initial permutation
c_start = MPI_Wtime()
! Allocate record arrays
allocate(history(1+args_obj%PERM_ITER_LIM+args_obj%ITER_LIM))
allocate(r_col_array(p))
! Allocate input and output arrays
allocate(index_level(M-1))
allocate(index_pair(M/2,2,N))
allocate(COEFF(M,M))
index_level = 0
index_pair = 0
COEFF = 0.0d0
! index_level and index_pair are invariant per CSD invocation
call CYG_INDEXTABLE(N,M,index_level,index_pair)
! COEFF is also invariant per CSD invocation
call CYGC_COEFF(N,M,COEFF)
! Allocate workspace arrays
allocate(Perm(M))
allocate(Perm_new(M))
allocate(Perm_sol(M))
allocate(QPerm(N))
allocate(QPerm_new(N))
! Start with the identity permutation
do i = 1, M
    Perm(i) = i
    Perm_new(i) = i
    Perm_sol(i) = i
end do
! Initialize all process to the identity qubit permutation first - however, do not perform qubit permutation selection on root process
do i = 1, N
    QPerm(i) = i
end do
! Run object constructors
allocate(type_spec(nset))
do i = 1, nset
    ! Set all matrices to be of real type
    type_spec(i) = 0
end do
! Set the middle matrix to be complex
type_spec(3) = 1
call csdgen_obj%constructor(N,M,1,index_level,index_pair,COEFF)
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
call qperm_process_CPLX(N,M,csdss_Xinit,csdgen_obj,QPerm,X,ecur)
! Initialize history stuffs
history(1) = ecur
hist_idx = 2
! Repeatedly find random permutations until the number of gates has been lowered relative to the original number of gates (i.e. the identity qubit permutation)
if(my_rank /= root) then
    do i = 1, args_obj%PERM_ITER_LIM
        call qperm_generate(N,QPerm_new)
        call qperm_process_CPLX(N,M,csdss_Xinit,csdgen_obj,QPerm_new,X,enew)
        if(enew <= ecur) then
            QPerm = QPerm_new
            ecur = enew
        end if
        history(hist_idx) = ecur
        hist_idx = hist_idx + 1
    end do
    ! Rerun qperm_process_CPLX to ensure consistency with the optimal QPerm
    call qperm_process_CPLX(N,M,csdss_Xinit,csdgen_obj,QPerm,X,ecur)
else
    hist_idx = hist_idx + args_obj%PERM_ITER_LIM
    history(2:hist_idx) = history(1)
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
    history(hist_idx) = ecur
    hist_idx = hist_idx + 1
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
MPI_int_buffer(1) = esol
call MPI_Allgather(MPI_int_buffer,1,MPI_INTEGER,r_col_array,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
m_pos = FindMinPos(p,r_col_array)
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
    write(1,'(i8,a,i8)')args_obj%PERM_ITER_LIM," ",args_obj%ITER_LIM
    do i = 0, args_obj%PERM_ITER_LIM + args_obj%ITER_LIM
        write(1,'(i8,a,i8)')i," ",history(i+1)
    end do
    close(1)
    ! Write the optimal qubit permutation list q and the permutation list p
    open(unit=2,file=fperm,action='write')
    do i = 1, N
        write(2,'(i15,a)',advance='no')QPerm(i)," "
    end do
    write(2,*)
    do i = 1, M
        write(2,'(i15,a)',advance='no')Perm_sol(i)," "
    end do
    write(2,*)
    close(2)
    ! Write the .tex files
    call csdss_Xsol%write_circuit(csdwh_obj)
end if

call flush(6)
call MPI_Barrier(MPI_COMM_WORLD,ierr)
MPI_int_buffer(1) = 0
! Assuming root = 0
if(my_rank /= root) then
    call MPI_Recv(MPI_int_buffer,1,MPI_INTEGER,my_rank-1,101,MPI_COMM_WORLD,mpi_stat,ierr)
end if
write(*,'(a,i4,a,i8,a,i8)')"Process ",my_rank,": Initial number of gates before/after reduction = ",csdss_Xinit%csd_ss_ct,"/",csdss_Xinit%csdr_ss_ct
write(*,'(a,i4,a,i8,a,i8,a,i8,a,i8,a,i8,a,i8,a,i8)')"Process ",my_rank,": Breakdown of solution = ",csdss_Xsol%arr(1)%csdr_ct,"/",csdss_Xsol%arr(2)%csdr_ct,"/",csdss_Xsol%arr(3)%csdr_ct,"/",csdss_Xsol%arr(4)%csdr_ct,"/",csdss_Xsol%arr(5)%csdr_ct,"/",csdss_Xsol%csdr_ss_ct," at step number ",idx_sol
write(*,'(a,i4,a,f15.9,a,f15.9,a)')"Process ",my_rank,": Time taken = ",c_time1,"/",c_time2," seconds."
call flush(6)
if(my_rank /= p-1) then
    call MPI_Send(MPI_int_buffer,1,MPI_INTEGER,my_rank+1,101,MPI_COMM_WORLD,ierr)
end if

! Free up memory
deallocate(MPI_int_buffer)
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
deallocate(QPerm_new)
deallocate(type_spec)

! Run object destructors
call csdgen_obj%destructor()
call csdss_Xinit%destructor()
call csdss_Xcur%destructor()
call csdss_Xnew%destructor()
call csdss_Xsol%destructor()
call csdwh_obj%destructor()

end subroutine OptQC_CPLX

subroutine Box_CPLX(M,source,dest)

implicit none
integer :: M
double complex :: source(M*M), dest(M,M)

integer :: i, j, ct

ct = 1
do i = 1, M
    do j = 1, M
        dest(i,j) = source(ct)
        ct = ct + 1
    end do
end do
return

end subroutine Box_CPLX

subroutine Flatten_CPLX(M,source,dest)

implicit none
integer :: M
double complex :: source(M,M), dest(M*M)

integer :: i, j, ct

ct = 1
do i = 1, M
    do j = 1, M
        dest(ct) = source(i,j)
        ct = ct + 1
    end do
end do
return

end subroutine Flatten_CPLX
