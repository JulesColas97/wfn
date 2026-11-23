module 3d_linear_interp
    implicit none

    ! Define the generic interface
    interface bilinear_interp
        module procedure interp_real_1d, interp_real_2d
        module procedure interp_cmplx_1d, interp_cmplx_2d
    end interface

    private :: interp_core_real

contains



end module
