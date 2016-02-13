module memwork_real

use arrays_real
use common_module

implicit none
! Input arrays
integer, allocatable :: index_level(:), index_pair(:,:,:)
! Output arrays
double precision, allocatable :: GATEY(:,:), GATEPI(:,:)
! Workspace arrays
! CYGR_CSD:
type(l_arr_dp_4) :: Z0, Z1
! CYGR_CUTGATE:
integer, allocatable :: Z_array(:,:), PI_array(:,:), GATEY_sign(:)
! CYGR_BLKCSD:
type(l_arr_dp_2) :: X_blk, X11, X12, X21, X22, U1, U2, V1T, V2T
type(l_arr_dp_1) :: GATEY_blk
type(l_arr_int_1) :: IWORK

end module memwork_real
