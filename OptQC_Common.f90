! Functions/Subroutines for Main

function CalcTol(ecur)

implicit none
integer :: ecur

integer :: CalcTol

CalcTol = ceiling(0.01d0 * ecur)
return

end function CalcTol

function FindMinPos(p,arr)

implicit none
integer :: p
integer :: arr(p)

integer :: i, curmin
integer :: FindMinPos

FindMinPos = 1
curmin = arr(1)
do i = 2, p
    if(arr(i) < curmin) then
        curmin = arr(i)
        FindMinPos = i
    end if
end do
return

end function FindMinPos

subroutine init_random_seed()

implicit none
integer :: i, n, clock
integer, allocatable :: seed(:)

call random_seed(size = n)
allocate(seed(n))
call system_clock(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(PUT = seed)
deallocate(seed)

end subroutine init_random_seed

! Functions/Subroutines for Perm

! Find the first position (from the left) where an element of used is zero.
function FindUnused(M,used)

implicit none
integer :: M
integer :: used(M)

integer :: i, FindUnused

FindUnused = 1
do i = 1, M
    if(used(i) == 0) exit
    FindUnused = FindUnused + 1
end do

return

end function FindUnused

! Returns the Gray code distance between x and y up to N bits
function GrayCodeDistance(N,x,y)

implicit none
integer :: N, x, y

integer :: i, xor, GrayCodeDistance

GrayCodeDistance = 0
xor = IEOR(x,y)
do i = 1, N
    GrayCodeDistance = GrayCodeDistance + mod(xor,2)
    xor = rshift(xor,1)
end do

return

end function GrayCodeDistance

! Generates a random number RINT such that 1 <= RINT <= maxv
function RINT(maxv)

implicit none
integer :: maxv

integer :: RINT
real :: temp

call random_number(temp)
RINT = floor( maxv * temp ) + 1

return

end function RINT

! Functions/Subroutines for Output

subroutine GroupGates(l,extent,Gate_Param,Sig_Offset,C_Num_Bin,Type_Param,N_Per_Type,N_Type,N_Total)

implicit none
integer :: l
integer :: extent ! Length of Gate_Param section - variable for GATEPI
double precision :: Gate_Param(l)
double precision :: Sig_Offset ! -1 for GatePI, 0 for GATEY
character(len=20) :: C_Num_Bin(l,l) ! Assume N-1 <= 20 - generally a safe assumption.....
double precision :: Type_Param(l)
integer :: N_Per_Type(l)
integer :: N_Type
integer :: N_Total

double precision :: ref, tol
integer :: i, j, match_idx

tol = 0.00000001d0
N_Per_Type = 0
N_Type = 0
N_Total = 0
do i = 1, extent
    ref = Gate_Param(i)
    if(abs(ref+Sig_Offset) > tol) then
        match_idx = 0
        do j = 1, N_Type
            if(abs(ref-Type_Param(j)) <= tol) then
                match_idx = j
                exit
            end if
        end do
        if(match_idx == 0) then
            N_Type = N_Type + 1
            Type_Param(N_Type) = ref
            match_idx = N_Type
        end if
        N_Per_Type(match_idx) = N_Per_Type(match_idx) + 1
        write(C_Num_bin(match_idx,N_Per_Type(match_idx)),"(B20.20)")i-1
        N_Total = N_Total + 1
    end if
end do
return

end subroutine GroupGates

subroutine ReduceGroups(l,N,C_Num_Bin,Type_Param,N_Per_Type,N_Type,N_Total)

implicit none
integer :: l
integer :: N
character(len=20) :: C_Num_Bin(l,l)
double precision :: Type_Param(l)
integer :: N_Per_Type(l)
integer :: N_Type
integer :: N_Total

character(len=20) :: workstr
integer :: i, j, k, row, p, limit
integer :: QSimilar, BitPos

do i = 1, N_Type
    if(N_Per_Type(i) == 1) cycle
    do row = 1, N-1
        j = 1
        do while(j < N_Per_Type(i))
            limit = N_Per_Type(i)
            do k = j+1, limit
                if(QSimilar(N,C_Num_Bin(i,j),C_Num_Bin(i,k),row) == 1) then
                    !write(*,*)"Combining ",C_Num_Bin(i,j)," and ",C_Num_Bin(i,k),"."
                    p = BitPos(N,row)
                    C_Num_Bin(i,j)(p:p) = '*'
                    !write(*,*)"Result: ",C_Num_Bin(i,j)
                    if(k /= limit) C_Num_Bin(i,k) = C_Num_Bin(i,limit)
                    N_Per_Type(i) = limit - 1
                    N_Total = N_Total - 1
                    exit
                end if
            end do
            j = j+1
        end do
    end do
end do
return

end subroutine ReduceGroups

function IsEmpty(workstr,lp,rp)

implicit none
character(len=20) :: workstr
integer :: lp, rp

integer :: i
integer :: IsEmpty

IsEmpty = 1
do i = lp, rp
    if(workstr(i:i)=='0' .or. workstr(i:i)=='1') then
        IsEmpty = 0
        return
    end if
end do
return

end function IsEmpty

function QSimilar(N,str1,str2,e_pos)

implicit none
integer :: N ! Keep in mind that str1 and str2 should be of extent N-1
character(len=20) :: str1
character(len=20) :: str2
integer :: e_pos

integer :: i, p
integer :: QSimilar ! 0 if not similar, 1 if similar
integer :: BitPos

QSimilar = 1
! Must differ at the specified e_pos bit
p = BitPos(N,e_pos)
if(str1(p:p) == str2(p:p)) then
    QSimilar = 0
    return
end if
! Must not differ at elsewhere
if(e_pos > 1) then
    do i = 1, e_pos-1
        p = BitPos(N,i)
        if(str1(p:p) /= str2(p:p)) then
            QSimilar = 0
            return
        end if
    end do
end if
if(e_pos < N) then
    do i = e_pos+1, N-1
        p = BitPos(N,i)
        if(str1(p:p) /= str2(p:p)) then
            QSimilar = 0
            return
        end if
    end do
end if
return

end function QSimilar

function BitPos(N,idx)

implicit none
integer :: N
integer :: idx

integer :: BitPos

BitPos = 21-N+idx

end function BitPos
