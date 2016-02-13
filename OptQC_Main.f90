subroutine OptQC_REAL(root,my_rank,p,args_obj)

use common_module
use csd_mpi
use csd_perm
use csd_tools
use m_mrgrnk
use mpi
use rng

implicit none
character(len=128) :: finput, fhist, fperm, ftex
integer :: N, M, d, i, j, k, l, Msq, num
integer :: ecur, enew, esol, idx_sol
integer, allocatable :: history(:)
integer :: hist_idx, m_pos, delta, tol, t0, e0
! Input arrays
double precision, allocatable :: X(:,:)
! Functions
integer :: FindMinPos, CalcTol, ChooseN
! Workspace arrays
integer, allocatable :: index_level(:), index_pair(:,:,:)
double precision, allocatable :: COEFF(:,:)
integer, allocatable :: Perm(:), Perm_new(:), Perm_sol(:)
integer, allocatable :: QPerm(:), QPerm_new(:)

! MPI variables
integer :: my_rank, p, ierr, root
integer :: MPI_int_buffer(1)
double precision, allocatable :: X_MPI_temp(:)
integer, allocatable :: r_col_array(:)
integer :: mpi_stat(MPI_STATUS_SIZE)
integer :: mpi_group_world, mpi_group_local, mpi_comm_local
! Timer variables
double precision :: c_start, c_mid, c_end, c_time1, c_time2, c_t1, c_t2, c_tot
! Synchronization variables
integer :: SYNCH_ct, POPTNUM
integer :: p_group
integer, allocatable :: p_rank_pos(:), groupsizes(:)
type(l_arr_int_1) :: subgroups
character, pointer :: synch_buffer(:)

! Objects
integer :: nset
integer, allocatable :: type_spec(:)
type(prog_args) :: args_obj
type(csd_generator) :: csdgen_obj
type(csd_solution_set) :: csdss_Xinit, csdss_Xcur, csdss_Xnew, csdss_Xsol
type(csd_write_handle) :: csdwh_obj
type(csd_mpi_handle) :: csdmh_obj

! Initialize the RNG object from the module
call rng_inst%seed(my_rank)
! Set the number of matrices to be 5
nset = 5
! Get group handle to MPI_COMM_WORLD
call MPI_Comm_group(MPI_COMM_WORLD,mpi_group_world,ierr)

if(my_rank == root) then
    write(*,'(a)')"Performing decomposition of a real orthogonal matrix."
    write(finput,'(a,a)')args_obj%fbase(1:args_obj%flength),".txt"
    write(*,'(a,a)')"Retrieving unitary matrix from file ",finput
    open(unit=1,file=finput,action='read')
    read(1,*)d
    N = ChooseN(d)
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
    		X(i,i) = 1.0d0
    	end do
    end if
    close(1)
    MPI_int_buffer(1) = N
    call Flatten(M,X,X_MPI_temp)
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
call MPI_Bcast(X_MPI_temp,Msq,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
if(my_rank /= root) then
    call Box(M,X_MPI_temp,X)
end if

! Set up modules
call initialize_csd_perm(N,M)

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
call csdgen_obj%constructor(N,M,0,index_level,index_pair,COEFF)
call csdss_Xinit%constructor(N,M,nset,type_spec)
csdss_Xinit%arr(1)%toggle_csd = .false.
csdss_Xinit%arr(2)%toggle_csd = .false.
csdss_Xinit%arr(4)%toggle_csd = .false.
csdss_Xinit%arr(5)%toggle_csd = .false.
call csdss_Xcur%constructor(N,M,nset,type_spec)
call csdss_Xnew%constructor(N,M,nset,type_spec)
call csdss_Xsol%constructor(N,M,nset,type_spec)
call csdwh_obj%constructor(N)
call csdmh_obj%constructor(N,M)
! Convention: U = Q^T P^T U' P Q - CHECKED AND VERIFIED
call permlisttomatrixtr(M,Perm,csdss_Xinit%arr(2)%X)                            ! P^T
call permlisttomatrix(M,Perm,csdss_Xinit%arr(4)%X)                              ! P
call qperm_process(N,M,csdss_Xinit,csdgen_obj,QPerm,X,ecur)
! Initialize history stuffs
history(1) = ecur
hist_idx = 2
! Repeatedly find random permutations until the number of gates has been lowered relative to the original number of gates (i.e. the identity qubit permutation)
if(my_rank /= root) then
    do i = 1, args_obj%PERM_ITER_LIM
        call qperm_generate(N,QPerm_new)
        call qperm_process(N,M,csdss_Xinit,csdgen_obj,QPerm_new,X,enew)
        if(enew <= ecur) then
            QPerm = QPerm_new
            ecur = enew
        end if
        history(hist_idx) = ecur
        hist_idx = hist_idx + 1
    end do
    ! Rerun qperm_process to ensure consistency with the optimal QPerm
    call qperm_process(N,M,csdss_Xinit,csdgen_obj,QPerm,X,ecur)
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
SYNCH_ct = 0

! Note: Might change this later - currently hardcoded to be 10% of the number of processes
POPTNUM = ceiling(0.1d0 * p)
! Allocate synchronization variables
allocate(p_rank_pos(p))
allocate(groupsizes(POPTNUM))
call subgroups%constructor(POPTNUM)
! Determine group sizes
i = floor(p/dble(POPTNUM))
j = p - (i * POPTNUM)
groupsizes = i
if(j > 0) groupsizes(1:j) = i+1
do i = 1, POPTNUM
	call subgroups%l(i)%constructor(groupsizes(i))
end do

c_tot = 0.0d0
do i = 1, args_obj%ITER_LIM
    Perm_new = Perm
    call NeighbourhoodOpt(N,M,csdss_Xinit,csdss_Xcur,csdss_Xnew,Perm_new)
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
    ! Perform synchronization and comparisons after every SYNCH_NUM iterations
    SYNCH_ct = SYNCH_ct + 1
    if(SYNCH_ct == args_obj%SYNCH_NUM .and. POPTNUM < p) then
    	c_t1 = MPI_Wtime()
    	! Wait for all processes to reach this point - possible slowdown issue between fast/slow threads
    	call MPI_Barrier(MPI_COMM_WORLD,ierr)
    	! Collate all current objective function values
    	MPI_int_buffer(1) = ecur
    	call MPI_Allgather(MPI_int_buffer,1,MPI_INTEGER,r_col_array,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    	! Obtain ranking for each process
    	p_rank_pos = -1
    	call mrgrnk(r_col_array,p_rank_pos)
    	! Construct groupings - note index i is already used >_>
    	l = POPTNUM+1
    	do j = 1, POPTNUM
    		num = p_rank_pos(j)
    		if(my_rank+1 == num) p_group = j
    		subgroups%l(j)%arr(1) = num-1
    		do k = 2, groupsizes(j)
    			num = p_rank_pos(l)
    			if(my_rank+1 == num) p_group = j
    			subgroups%l(j)%arr(k) = num-1
    			l = l + 1
    		end do
    	end do
        ! Get handle to the new group
        call MPI_Group_incl(mpi_group_world,groupsizes(p_group),subgroups%l(p_group)%arr,mpi_group_local,ierr)
        ! Create communicator for the new group
        call MPI_Comm_create(MPI_COMM_WORLD,mpi_group_local,mpi_comm_local,ierr)
        ! Broadcast optimal solutions from root processes to the rest of the local groups - note that root process (the most optimal) is always the rank 0 in each group
        call csdmh_obj%synchronize(QPerm,Perm,csdss_Xcur,mpi_group_local,mpi_comm_local,0)
!       ! Check difference
!       Xcpy = matmul(csdss_Xcur%arr(2)%X,csdss_Xcur%arr(3)%X)
! 		Xcpy = matmul(csdss_Xcur%arr(1)%X,Xcpy)
! 		Xcpy = matmul(Xcpy,csdss_Xcur%arr(4)%X)
! 		Xcpy = matmul(Xcpy,csdss_Xcur%arr(5)%X)
! 		write(*,*)"Process ",my_rank," has difference matrix: ",Xcpy-X
        ! Deallocate the new group and corresponding communicator
        call MPI_Group_free(mpi_group_local,ierr)
        call MPI_Comm_free(mpi_comm_local,ierr)
    	! Reset synchronization count
    	SYNCH_ct = 0
    	c_t2 = MPI_Wtime()
    	c_tot = c_tot + c_t2 - c_t1
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
write(*,'(a,i4,a,f15.9,a,f15.9,a,f15.9,a)')"Process ",my_rank,": Time taken = ",c_time1,"/",c_time2,"/",c_tot," seconds."
call flush(6)
if(my_rank /= p-1) then
    call MPI_Send(MPI_int_buffer,1,MPI_INTEGER,my_rank+1,101,MPI_COMM_WORLD,ierr)
end if

! Finalize modules
call finalize_csd_perm()

! Deallocate group handles
call MPI_Group_free(mpi_group_world,ierr)

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
deallocate(QPerm_new)
deallocate(type_spec)
deallocate(p_rank_pos)
deallocate(groupsizes)

! Run object destructors
call subgroups%destructor()
call csdgen_obj%destructor()
call csdss_Xinit%destructor()
call csdss_Xcur%destructor()
call csdss_Xnew%destructor()
call csdss_Xsol%destructor()
call csdwh_obj%destructor()
call csdmh_obj%destructor()

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
