program OptQC

use common_module
use mpi

implicit none
! Program arguments
type(prog_args) :: args_obj
character(len=128) :: fin
integer :: fdum, fstat
! MPI variables
integer(4) :: root, my_rank, p, ierr

! Initialize MPI environment
call MPI_Init(ierr)
! Find process rank
call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
! Find communicator size
call MPI_Comm_size(MPI_COMM_WORLD,p,ierr)
! Choose the root process as 0
root = 0

! Read in matrix from file specified by command line argument, which does not include .txt extension
call get_command_argument(1,args_obj%fbase,args_obj%flength,fstat)
if(len_trim(args_obj%fbase) == 0) then
    write(*,'(a)')"No file specified. Exiting program."
    call exit(1)
end if
if(fstat > 0) then
    write(*,'(a)')"Unable to retrieve filename from command line arguments. Exiting program."
    call exit(1)
end if
if(fstat == -1) then
    write(*,'(a)')"Filename specified on the argument list is too long. Exiting program."
    call exit(1)
end if
! Take parameter values from command-line arguments
call get_command_argument(2,fin,fdum,fstat)
read(fin,*)args_obj%PROG_TYPE
call get_command_argument(3,fin,fdum,fstat)
read(fin,*)args_obj%ITER_LIM
call get_command_argument(4,fin,fdum,fstat)
read(fin,*)args_obj%PERM_ITER_LIM
call get_command_argument(5,fin,fdum,fstat)
read(fin,*)args_obj%TOL_COEFF

if(my_rank == root) then
    ! Output of obtained parameters
    write(*,*)
    write(*,'(a,a,a,i2,a,i8,a,i8,a,f9.6)')"Command-line arguments: ",args_obj%fbase(1:args_obj%flength)," ",args_obj%PROG_TYPE," ",args_obj%ITER_LIM," ",args_obj%PERM_ITER_LIM," ",args_obj%TOL_COEFF
    write(*,'(a,i4)')"Size of communicator: ",p
    write(*,*)
end if

if(args_obj%PROG_TYPE == 0) then
    call OptQC_REAL(root,my_rank,p,args_obj)
else
    call OptQC_CPLX(root,my_rank,p,args_obj)
end if

! Finalize MPI environment
call MPI_Finalize(ierr)

end program OptQC
