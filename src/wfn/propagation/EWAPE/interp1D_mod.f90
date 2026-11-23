module interp_1D_mod
    use kind_mod, only: wp
    implicit none
    private
    public :: linear_interp

    ! Generic interface: resolved by the type of `values`
    interface linear_interp
        module procedure interp_core_real
        module procedure interp_core_complex
        module procedure interp_core_real_nonuniform
        module procedure interp_core_complex_nonuniform
    end interface

contains

    ! --- uniform grid, REAL ---------------------------------------------------
    subroutine interp_core_real(values, x_min, dx, x, res)
        real(wp), intent(in)  :: values(:), x_min, dx, x(:)
        real(wp), intent(out) :: res(:)

        integer  :: i, ix, nx
        real(wp) :: x_rel, fx, inv_dx

        nx = size(values)
        inv_dx = 1.0_wp / dx
        res = 0.0_wp                          ! points outside the grid stay 0

        do i = 1, size(x)
            x_rel = (x(i) - x_min) * inv_dx
            ix = int(floor(x_rel)) + 1
            if (ix < 1 .or. ix > nx-1) cycle
            fx = x_rel - floor(x_rel)
            res(i) = (1.0_wp - fx) * values(ix) + fx * values(ix+1)
        end do
    end subroutine

    ! --- uniform grid, COMPLEX ------------------------------------------------
    subroutine interp_core_complex(values, x_min, dx, x, res)
        complex(wp), intent(in)  :: values(:)
        real(wp),    intent(in)  :: x_min, dx, x(:)
        complex(wp), intent(out) :: res(:)

        integer  :: i, ix, nx
        real(wp) :: x_rel, fx, inv_dx

        nx = size(values)
        inv_dx = 1.0_wp / dx
        res = (0.0_wp, 0.0_wp)

        do i = 1, size(x)
            x_rel = (x(i) - x_min) * inv_dx
            ix = int(floor(x_rel)) + 1
            if (ix < 1 .or. ix > nx-1) cycle
            fx = x_rel - floor(x_rel)
            res(i) = (1.0_wp - fx) * values(ix) + fx * values(ix+1)
        end do
    end subroutine

    ! --- non-uniform grid, REAL -----------------------------------------------
    subroutine interp_core_real_nonuniform(values, x1, x, res)
        real(wp), intent(in)  :: values(:), x1(:), x(:)
        real(wp), intent(out) :: res(:)

        integer  :: i, ix, nx, low, high, mid
        real(wp) :: fx

        nx = size(values)
        res = 0.0_wp

        do i = 1, size(x)
            if (x(i) < x1(1) .or. x(i) > x1(nx)) cycle
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
            res(i) = (1.0_wp - fx) * values(ix) + fx * values(ix+1)
        end do
    end subroutine

    ! --- non-uniform grid, COMPLEX --------------------------------------------
    subroutine interp_core_complex_nonuniform(values, x1, x, res)
        complex(wp), intent(in)  :: values(:)
        real(wp),    intent(in)  :: x1(:), x(:)
        complex(wp), intent(out) :: res(:)

        integer  :: i, ix, nx, low, high, mid
        real(wp) :: fx

        nx = size(values)
        res = (0.0_wp, 0.0_wp)

        do i = 1, size(x)
            if (x(i) < x1(1) .or. x(i) > x1(nx)) cycle
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
            res(i) = (1.0_wp - fx) * values(ix) + fx * values(ix+1)
        end do
    end subroutine

end module interp_1D_mod
