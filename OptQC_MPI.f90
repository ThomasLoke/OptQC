module csd_mpi

use csd_tools
use mpi

implicit none

type csd_mpi_handle
	integer :: N, M
	integer, pointer :: QPerm(:), Perm(:)
	type(csd_solution_set), pointer :: css
	character, allocatable :: buffer(:), circ(:)
	integer, allocatable :: css_info(:)
	integer :: bpos, blen
	integer :: mpi_group, mpi_comm, root, my_rank ! Note: my_rank is determined by the local group
contains
	procedure :: constructor => csd_mpi_handle_constructor
	procedure :: destructor => csd_mpi_handle_destructor
	procedure :: assign_target => csd_mpi_handle_assign_target
	procedure :: nullify_target => csd_mpi_handle_nullify_target
	procedure :: circ_pack_size => csd_mpi_handle_circ_pack_size
	procedure :: pack_size => csd_mpi_handle_pack_size
	procedure :: circ_pack => csd_mpi_handle_circ_pack
	procedure :: circ_unpack => csd_mpi_handle_circ_unpack
	procedure :: pack => csd_mpi_handle_pack
	procedure :: unpack => csd_mpi_handle_unpack
	procedure :: synchronize => csd_mpi_handle_synchronize
end type csd_mpi_handle

contains

subroutine csd_mpi_handle_constructor(this,N,M)

implicit none
class(csd_mpi_handle) :: this
integer :: N, M, my_rank

this%N = N
this%M = M
allocate(this%circ(CLEN*(M*M+1)*(N+1))) ! Maximum circuit size
this%bpos = 0
this%blen = 0
this%mpi_group = 0
this%mpi_comm = 0
this%root = 0
this%my_rank = -1

end subroutine csd_mpi_handle_constructor

subroutine csd_mpi_handle_destructor(this)

implicit none
class(csd_mpi_handle) :: this

nullify(this%QPerm)
nullify(this%Perm)
nullify(this%css)
deallocate(this%circ)
this%bpos = 0
this%blen = 0
this%mpi_group = 0
this%mpi_comm = 0
this%root = 0
this%my_rank = -1

end subroutine csd_mpi_handle_destructor

subroutine csd_mpi_handle_assign_target(this,QPerm,Perm,css_inst,group,comm,root)

implicit none
class(csd_mpi_handle) :: this
integer, target :: QPerm(:), Perm(:)
type(csd_solution_set), target :: css_inst
integer :: group, comm, root

integer :: ierr

this%QPerm => QPerm
this%Perm => Perm
this%css => css_inst
allocate(this%css_info(3*css_inst%nset))
this%bpos = 0
this%blen = 0
this%mpi_group = group
this%mpi_comm = comm
this%root = root
call MPI_Group_rank(group,this%my_rank,ierr)
if(this%my_rank == MPI_UNDEFINED) then
	write(*,'(a)')"Error! Attempted to synchronize with a process that is not a member of a specified group."
    call exit(1)
end if

end subroutine csd_mpi_handle_assign_target

subroutine csd_mpi_handle_nullify_target(this)

implicit none
class(csd_mpi_handle) :: this

nullify(this%QPerm)
nullify(this%Perm)
nullify(this%css)
deallocate(this%css_info)
this%bpos = 0
this%blen = 0
this%mpi_group = 0
this%mpi_comm = 0
this%root = 0
this%my_rank = -1

end subroutine csd_mpi_handle_nullify_target

subroutine csd_mpi_handle_circ_pack_size(this,idx)

implicit none
class(csd_mpi_handle) :: this
integer :: idx

integer :: temp, temp2, ierr

! Check if circuit is of non-zero size - if so then skip, otherwise, don't.
if(this%css%arr(idx)%csdr_ct > 0) then
	! Number of characters = (Length of each circuit element) * (Number of columns) * (Number of rows)
	temp = CLEN * this%css%arr(idx)%csdr_ct * (this%N+1)
	call MPI_Pack_size(temp,MPI_CHARACTER,this%mpi_comm,temp2,ierr)
	this%blen = this%blen + temp2
end if

end subroutine csd_mpi_handle_circ_pack_size

subroutine csd_mpi_handle_pack_size(this)

implicit none
class(csd_mpi_handle) :: this

integer :: nset, i, temp, ierr

!-----------------------------------------------------------------------------------------------------!
! Packing order:
! -To define matrices-
! 1) QPerm - list representation of Q
! 2) Perm - list representation of P
! -To fill the csd_solution_set object-
! 3) (csd_ct,csdr_ct,neg) * 5 - count and stat parameters of each csd_solution
! 4) Circuit * 5
!-----------------------------------------------------------------------------------------------------!
nset = this%css%nset
this%blen = 0
! Calculate required buffer size
! 1)
call MPI_Pack_size(this%N,MPI_INTEGER,this%mpi_comm,temp,ierr)
this%blen = this%blen + temp
! 2)
call MPI_Pack_size(this%M,MPI_INTEGER,this%mpi_comm,temp,ierr)
this%blen = this%blen + temp
! 3)
call MPI_Pack_size(3*nset,MPI_INTEGER,this%mpi_comm,temp,ierr)
this%blen = this%blen + temp
! 4)
do i = 1, nset
	call this%circ_pack_size(i)
end do

end subroutine csd_mpi_handle_pack_size

subroutine csd_mpi_handle_circ_pack(this,idx)

implicit none
class(csd_mpi_handle) :: this
integer :: idx

character(len=CLEN) :: temp
integer :: i, j, k, ct, ierr

! Check if circuit is of non-zero size - if so then skip, otherwise, don't.
if(this%css%arr(idx)%csdr_ct > 0) then
	! Fill the circuit buffer first
	ct = 0
	do i = 1, this%N+1
		do j = 1, this%css%arr(idx)%csdr_ct
			temp = this%css%arr(idx)%Circuit(i,j)
			do k = 1, CLEN
				this%circ(ct+k) = temp(k:k)
			end do
			ct = ct + CLEN
		end do
	end do
	! Push circuit buffer onto packing buffer
	call MPI_Pack(this%circ(1:ct),ct,MPI_CHARACTER,this%buffer,this%blen,this%bpos,this%mpi_comm,ierr)
end if

end subroutine csd_mpi_handle_circ_pack

subroutine csd_mpi_handle_circ_unpack(this,idx)

implicit none
class(csd_mpi_handle) :: this
integer :: idx

character(len=CLEN) :: temp
integer :: i, j, k, ct, ierr, tempi

this%css%arr(idx)%Circuit = ''
! Check if circuit is of non-zero size - if so then skip, otherwise, don't.
if(this%css%arr(idx)%csdr_ct > 0) then
	! Extract circuit buffer from packing buffer - assume csdr_ct has been set beforehand
	tempi = CLEN*this%css%arr(idx)%csdr_ct*(this%N+1)
	call MPI_Unpack(this%buffer,this%blen,this%bpos,this%circ(1:tempi),tempi,MPI_CHARACTER,this%mpi_comm,ierr)
	! Fill the circuit based on circuit buffer
	ct = 0
	do i = 1, this%N+1
		do j = 1, this%css%arr(idx)%csdr_ct
			temp = ''
			do k = 1, CLEN
				temp(k:k) = this%circ(ct+k)
			end do
			this%css%arr(idx)%Circuit(i,j) = temp
			ct = ct + CLEN
		end do
	end do
end if

end subroutine csd_mpi_handle_circ_unpack

subroutine csd_mpi_handle_pack(this)

implicit none
class(csd_mpi_handle) :: this

integer, allocatable :: css_info(:)
integer :: nset, i, temp, ierr

nset = this%css%nset
! Assume that this%blen has been set to the correct value
! Refer to csd_mpi_handle_pack_size for the packing order
this%bpos = 0
! Push information onto packing buffer
! 1)
call MPI_Pack(this%QPerm,this%N,MPI_INTEGER,this%buffer,this%blen,this%bpos,this%mpi_comm,ierr)
! 2)
call MPI_Pack(this%Perm,this%M,MPI_INTEGER,this%buffer,this%blen,this%bpos,this%mpi_comm,ierr)
! 3)
temp = 1
do i = 1, nset
	this%css_info(temp) = this%css%arr(i)%csd_ct
	this%css_info(temp+1) = this%css%arr(i)%csdr_ct
	if(this%css%arr(i)%neg == .true.) then
		this%css_info(temp+2) = 1
	else 
		this%css_info(temp+2) = 0
	end if
	temp = temp + 3
end do
call MPI_Pack(this%css_info,3*nset,MPI_INTEGER,this%buffer,this%blen,this%bpos,this%mpi_comm,ierr)
! 4)
do i = 1, nset
	call this%circ_pack(i)
end do

end subroutine csd_mpi_handle_pack

subroutine csd_mpi_handle_unpack(this)

implicit none
class(csd_mpi_handle) :: this

integer :: nset, i, temp, ierr

nset = this%css%nset
! Assume that this%blen has been set to the correct value, and that the buffer has been received
! Refer to csd_mpi_handle_pack_size for the packing order
this%bpos = 0
! Unpack information from buffer
! 1)
call MPI_Unpack(this%buffer,this%blen,this%bpos,this%QPerm,this%N,MPI_INTEGER,this%mpi_comm,ierr)
! 2)
call MPI_Unpack(this%buffer,this%blen,this%bpos,this%Perm,this%M,MPI_INTEGER,this%mpi_comm,ierr)
! 3) 
call MPI_Unpack(this%buffer,this%blen,this%bpos,this%css_info,3*nset,MPI_INTEGER,this%mpi_comm,ierr)
temp = 1
do i = 1, nset
	this%css%arr(i)%csd_ct = this%css_info(temp)
	this%css%arr(i)%csdr_ct = this%css_info(temp+1)
	if(this%css_info(temp+2) == 1) then
		this%css%arr(i)%neg = .true.
	else
		this%css%arr(i)%neg = .false.
	end if
	temp = temp + 3
end do
! 4)
do i = 1, nset
	call this%circ_unpack(i)
end do

end subroutine csd_mpi_handle_unpack

subroutine csd_mpi_handle_synchronize(this,QPerm,Perm,css_inst,group,comm,root)

implicit none
class(csd_mpi_handle) :: this
integer, target :: QPerm(:), Perm(:)
type(csd_solution_set), target :: css_inst
integer :: group, comm, root

integer :: N, M, MPI_int_buffer(1), ierr
integer, allocatable :: perm_q(:), perm_qinv(:)

N = css_inst%N
M = css_inst%M
! Allocate temporary variables
allocate(perm_q(M))
allocate(perm_qinv(M))
! Assign target
call this%assign_target(QPerm,Perm,css_inst,group,comm,root)
this%bpos = 0
this%blen = 0
if(this%my_rank == root) then
	! Calculate buffer size
	call this%pack_size()
	MPI_int_buffer(1) = this%blen
	! Perform buffer allocation
	allocate(this%buffer(this%blen))
end if
! Broadcast buffer size
call MPI_Bcast(MPI_int_buffer,1,MPI_INTEGER,root,comm,ierr)
if(this%my_rank /= root) then
	this%blen = MPI_int_buffer(1)
	allocate(this%buffer(this%blen))
end if
! Pack information into buffer on root process
if(this%my_rank == root) then
	call this%pack()
end if
! Broadcast buffer
call MPI_Bcast(this%buffer,this%blen,MPI_PACKED,root,comm,ierr)
! Unpack information from buffer
if(this%my_rank /= root) then
	! Recover U matrix - extremely inefficient, but too lazy to figure out the mess of permutations (and inverses), which takes another 20 lines of code or smth
	! U = Q^T P^T U' P Q
	if(css_inst%arr(3)%obj_type == 0) then
		css_inst%arr(3)%X = matmul(css_inst%arr(2)%X,css_inst%arr(3)%X)
		css_inst%arr(3)%X = matmul(css_inst%arr(1)%X,css_inst%arr(3)%X)
		css_inst%arr(3)%X = matmul(css_inst%arr(3)%X,css_inst%arr(4)%X)
		css_inst%arr(3)%X = matmul(css_inst%arr(3)%X,css_inst%arr(5)%X)
	else
		css_inst%arr(3)%Xc = matmul(css_inst%arr(2)%X,css_inst%arr(3)%Xc)
		css_inst%arr(3)%Xc = matmul(css_inst%arr(1)%X,css_inst%arr(3)%Xc)
		css_inst%arr(3)%Xc = matmul(css_inst%arr(3)%Xc,css_inst%arr(4)%X)
		css_inst%arr(3)%Xc = matmul(css_inst%arr(3)%Xc,css_inst%arr(5)%X)
	end if
	call this%unpack()
	! Translate QPerm and Perm (and U') into corresponding matrices
	call qpermtoperm(N,M,QPerm,perm_qinv)
	call invertperm(M,perm_qinv,perm_q) 									! Reminder: MAGIC
	call permlisttomatrixtr(M,perm_q,css_inst%arr(1)%X)						! Q^T
	call permlisttomatrixtr(M,Perm,css_inst%arr(2)%X)						! P^T
	call permlisttomatrix(M,Perm,css_inst%arr(4)%X)							! P
	call permlisttomatrix(M,perm_q,css_inst%arr(5)%X)						! Q
	! U' = P Q U Q^T P^T - also inefficient, but bleh, once in a blue moon pls
	if(css_inst%arr(3)%obj_type == 0) then
		css_inst%arr(3)%X = matmul(css_inst%arr(5)%X,css_inst%arr(3)%X)
		css_inst%arr(3)%X = matmul(css_inst%arr(4)%X,css_inst%arr(3)%X)
		css_inst%arr(3)%X = matmul(css_inst%arr(3)%X,css_inst%arr(1)%X)
		css_inst%arr(3)%X = matmul(css_inst%arr(3)%X,css_inst%arr(2)%X)
	else
		css_inst%arr(3)%Xc = matmul(css_inst%arr(5)%X,css_inst%arr(3)%Xc)
		css_inst%arr(3)%Xc = matmul(css_inst%arr(4)%X,css_inst%arr(3)%Xc)
		css_inst%arr(3)%Xc = matmul(css_inst%arr(3)%Xc,css_inst%arr(1)%X)
		css_inst%arr(3)%Xc = matmul(css_inst%arr(3)%Xc,css_inst%arr(2)%X)
	end if
end if
! Deallocate temporary variables
deallocate(perm_q)
deallocate(perm_qinv)
! Deallocate buffer
deallocate(this%buffer)
! Nullify target
call this%nullify_target()

end subroutine csd_mpi_handle_synchronize

end module csd_mpi