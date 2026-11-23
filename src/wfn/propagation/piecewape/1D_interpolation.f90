module fast_interp_1D
    implicit none

    ! Define the generic interface
    interface linear_interp
        module procedure interp_core_real
        ! module procedure interp_core_complex
        module procedure interp_core_real_nonuniform
        ! module procedure interp_core_complex_nonuniform
    end interface
contains

    ! --- Version for REAL (Double Precision) ---
    subroutine interp_core_real(values, x_min, dx, x, res)
        real, intent(in) :: x_min, dx
        real, intent(in) :: values(:)
        real, intent(in) :: x(:)
        real, intent(out) :: res(:)

        integer :: i, ix, nx, n_pts
        real :: x_rel, fx, inv_dx

        nx = size(values)
        n_pts = size(x)
        inv_dx = 1.0 / dx

        do i = 1, n_pts
            x_rel = (x(i) - x_min) * inv_dx
            ix = int(x_rel) + 1
            if (ix < 1 .or. ix >= nx-1 ) then
                ! res(i) = 0.0
                cycle
            end if

            fx = x_rel - floor(x_rel)
            res(i) = (1.0 - fx) * values(ix) + fx  * values(ix+1)   
        end do
    end subroutine 
    
    subroutine interp_core_real_nonuniform(values, x1, x, res)
        real, intent(in) :: values(:), x1(:)
        real, intent(in) :: x(:)
        real, intent(out) :: res(:)

        integer :: i, ix, nx, n_pts, low, high, mid
        real :: x_rel, fx, inv_dx

        nx = size(values)
        n_pts = size(x)

        do i = 1, n_pts
            if (x(i) < x1(1) .or. x(i) > x1(nx)) then
                cycle
            end if

            low = 1
            high = nx
            do while (high - low > 1)
                mid = (low + high) / 2
                if (x(i) >= x1(mid)) then
                    low = mid
                else
                    high = mid
                end if
            end do
            ix = low
            fx = (x(i) - x1(ix)) / (x1(ix+1) - x1(ix))
            res(i) = (1.0 - fx) * values(ix) + fx  * values(ix+1)   
        end do
    end subroutine 
    
    ! subroutine interp_core_complex_nonuniform(values, x1, x, res)
    !     complex, intent(in) :: values(:)
    !     real, intent(in) :: x(:), x1(:)
    !     complex, intent(out) :: res(:)
    !
    !     integer :: i, ix, nx, n_pts, low, high, mid
    !     real :: x_rel, fx, inv_dx
    !
    !     nx = size(values)
    !     n_pts = size(x)
    !     do i = 1, n_pts
    !         if (x(i) < x1(1) .or. x(i) > x1(nx)) then
    !             cycle
    !         end if
    !
    !         low = 1
    !         high = nx
    !         do while (high - low > 1)
    !             mid = (low + high) / 2
    !             if (x(i) >= x1(mid)) then
    !                 low = mid
    !             else
    !                 high = mid
    !             end if
    !         end do
    !         ix = low
    !         fx = (x(i) - x1(ix)) / (x1(ix+1) - x1(ix))
    !
    !         res(i) = (1.0 - fx) * values(ix) + fx  * values(ix+1)   
    !     end do
    ! end subroutine 
end module fast_interp_1D
