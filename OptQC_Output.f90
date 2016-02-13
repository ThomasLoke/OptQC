subroutine ReduceSolution(N, M, GATEY, GATEPI, index_level, Circuit, col)

implicit none
integer :: N
integer :: M
double precision :: GATEY(M/2,M-1)
double precision :: GATEPI(M/2,N)
integer :: index_level(M-1)
character(len=15) :: Circuit(N+1,M*M/2+M)
integer :: col

integer :: i, j, k, p, ref, row, lb, rb
character(len=20) :: control, workstr
integer :: IsEmpty

integer :: l
character(len=20), allocatable :: C_Num_Bin(:,:)
double precision, allocatable :: Type_Param(:)
integer, allocatable :: N_Per_Type(:)
integer :: N_Type
integer :: N_Total

l = M/2 ! Note: THIS DOES NOT CHANGE
Circuit = ''
allocate(C_Num_Bin(l,l))
allocate(Type_Param(l))
allocate(N_Per_Type(l))

! Form the circuit for the reduced solution.
! GATEPI
col = 0
do i = 1,N
    lb = 22-i
    call GroupGates(l,2**(i-1),GATEPI(:,i),-1.0d0,C_Num_Bin,Type_Param,N_Per_Type,N_Type,N_Total)
    call ReduceGroups(l,N,C_Num_Bin,Type_Param,N_Per_Type,N_Type,N_Total)
    ! There should only be one type, that is, Type_Param(1) = -1.0d0
    if(N_Type > 1) then
        write(*,*)"Error! GATEPI value different to -1.0d0 encountered. GATEPI section:"
        write(*,*)GATEPI(:,i)
        write(*,*)"Terminating program."
        call exit(1)
    end if
    do k = 1, N_Per_Type(1)
        col = col+1
        if(i > 1) then
            workstr = C_Num_Bin(1,k)
            do row = 1, i-1
                p = 21-i+row
                if(workstr(p:p)=='0') then
                    Circuit(row,col) = '& \ctrlo{1}'
                else if(workstr(p:p)=='1') then
                    Circuit(row,col) = '& \ctrl{1}'
                else
                    ! Case where workstr(p:p) = '*'
                    if(row > 1 .and. IsEmpty(workstr,lb,p-1)==0) then
                        Circuit(row,col) = '& \qw \qwx[1]'
                    else
                        Circuit(row,col) = '& \qw'
                    end if
                end if
            end do
        end if
        Circuit(i,col) = '& \gate{\pi}'
        if(i < N) then
            do row = i+1,N
                Circuit(row,col) = '& \qw'
            end do
        end if
        Circuit(N+1,col) = '&'
    end do
end do
!GATEY
lb = 22-N
rb = 20
do i = M-1, 1, -1
    ref = index_level(i)
    call GroupGates(l,l,GATEY(:,i),0.0d0,C_Num_Bin,Type_Param,N_Per_Type,N_Type,N_Total)
    call ReduceGroups(l,N,C_Num_Bin,Type_Param,N_Per_Type,N_Type,N_Total)
    do j = 1, N_Type
        do k = 1, N_Per_Type(j)
            col = col+1
            workstr = C_Num_Bin(j,k)
            Circuit(ref,col) = '& \gate{R_y}'
            do row = 1,N-1
                p = 21-N+row
				if(row < ref) then
					if (workstr(p:p)=='0') then
						Circuit(row,col) = '& \ctrlo{1}'
					else if(workstr(p:p)=='1') then
						Circuit(row,col) = '& \ctrl{1}'
                    else
                        ! Case where workstr(p:p) = '*'
                        if(row > 1 .and. IsEmpty(workstr,lb,p-1)==0) then
                            Circuit(row,col) = '& \qw \qwx[1]'
                        else
                            Circuit(row,col) = '& \qw'
                        end if
					end if
				else
					if(workstr(p:p)=='0') then
						Circuit(row+1,col) = '& \ctrlo{-1}'
					else if(workstr(p:p)=='1') then
						Circuit(row+1,col) = '& \ctrl{-1}'
                    else
                        ! Case where workstr(p:p) = '*'
                        if(row < N-1 .and. IsEmpty(workstr,p+1,rb)==0) then
                            Circuit(row+1,col) = '& \qw \qwx[-1]'
                        else
                            Circuit(row+1,col) = '& \qw'
                        end if
					end if
				end if
			end do
			write(Circuit(N+1,col),302)Type_Param(j)
        end do
    end do
end do
302 FORMAT('& ',F8.4)

deallocate(C_Num_Bin)
deallocate(Type_Param)
deallocate(N_Per_Type)
return

end subroutine ReduceSolution
