module constants
    use kind_mod, only: wp
    implicit none
    real(wp), PARAMETER :: pi = 4.0_wp * atan(1.0_wp)
    real(wp), parameter :: kappa = 0.41_wp
    complex(wp), parameter :: eye = (0.0_wp, 1.0_wp)
    integer, parameter :: UNSET = -huge(1)   ! sentinel
end module 

