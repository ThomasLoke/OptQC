module rng
  
implicit none

private :: ns, default_seed
public :: rng_t, rng_inst

! Dimension of the state
integer, parameter :: ns = 4
! Default seed vector
integer, parameter, dimension(ns) :: default_seed = (/ 521288629, 362436069, 16163801, 1131199299 /)

! A data type for storing the state of the RNG
type rng_t
    integer, dimension(ns) :: state = default_seed
contains
    procedure :: seed => rng_t_seed
    procedure :: r01 => rng_t_r01
    procedure :: rint => rng_t_rint
end type rng_t

! Globally accessible copy of the RNG object
type(rng_t) :: rng_inst

contains

! Seeds the RNG using a single integer and a default seed vector.
subroutine rng_t_seed(this,seed)

implicit none
class(rng_t), intent(inout) :: this
integer, intent(in) :: seed

integer :: clock, i

call system_clock(COUNT=clock)
this%state = clock / (seed + (/ (i, i = 1, ns) /))

end subroutine rng_t_seed

! Draws a uniform real number on [0,1].
function rng_t_r01(this) result(u)

implicit none
class(rng_t), intent(inout) :: this
double precision :: u
integer :: imz

imz = this%state(1) - this%state(3)
if (imz < 0) imz = imz + 2147483579
this%state(1) = this%state(2)
this%state(2) = this%state(3)
this%state(3) = imz
this%state(4) = 69069 * this%state(4) + 1013904243
imz = imz + this%state(4)
u = 0.5d0 + 0.23283064d-9 * imz

end function rng_t_r01

! Draws a uniform real number on [0,1] and then converts it to an integer on [1,maxv]
function rng_t_rint(this,maxv) result(u)

implicit none
class(rng_t), intent(inout) :: this
integer :: maxv
integer :: u
double precision :: rval

rval = this%r01()
u = floor( maxv * rval ) + 1

end function rng_t_rint

end module rng