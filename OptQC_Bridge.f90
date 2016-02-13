program OptQC

use mpi

implicit none
character(len=128) :: fbase, fin
integer :: flength, fdum, fstat
integer :: PROG_TYPE, ITER_LIM
! MPI variables
integer(4) :: my_rank, p, ierr

! Initialize MPI environment
call MPI_Init(ierr)
! Find process rank
call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
! Find communicator size
call MPI_Comm_size(MPI_COMM_WORLD,p,ierr)

! Read in matrix from file specified by command line argument, which does not include .txt extension
call get_command_argument(1,fbase,flength,fstat)
! Take parameter values from command-line arguments
call get_command_argument(2,fin,fdum,fstat)
read(fin,*)PROG_TYPE
call get_command_argument(3,fin,fdum,fstat)
read(fin,*)ITER_LIM

if(len_trim(fbase) == 0) then
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

if(PROG_TYPE == 0) then
    call OptQC_REAL(my_rank,p,fbase,flength,ITER_LIM)
else
    !call OptQC_CPLX(my_rank,p,fbase,flength,ITER_LIM)
end if

! Finalize MPI environment
call MPI_Finalize(ierr)

end program OptQC
