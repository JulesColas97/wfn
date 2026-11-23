module interp_plane_mod
    use kind_mod, only: wp
    implicit none
    private
    public :: interp_to_plane

    ! One call, resolved by the rank of `values`:
    !   1D : interp_to_plane(f1, xg, dim, r, z, res [, clamp])
    !   2D : interp_to_plane(f2, xg, zg, r, z, res [, clamp])
    !   3D : interp_to_plane(f3, xg, yg, zg, r, z, tau, Xs, Ys, res [, clamp])
    ! res(size(r), size(z)). All source grids must be uniform and ascending.
    ! clamp = .true. (default): points outside the grid take the nearest edge value
    ! clamp = .false.         : points outside the grid are set to 0
    interface interp_to_plane
        module procedure interp_1d_to_plane
        module procedure interp_2d_to_plane
        module procedure interp_3d_to_plane
    end interface

contains

    ! --- bracketing index and weight on a UNIFORM ascending grid -------------
    pure subroutine locate(x, xq, clamp, i, f, ok)
        real(wp), intent(in)  :: x(:), xq
        logical,  intent(in)  :: clamp
        integer,  intent(out) :: i
        real(wp), intent(out) :: f
        logical,  intent(out) :: ok
        integer  :: n
        real(wp) :: s

        n = size(x)
        s = (xq - x(1)) / (x(2) - x(1))       ! fractional index, 0-based
        if (clamp) then
            s  = min(max(s, 0.0_wp), real(n-1, wp))
            ok = .true.
        else
            ok = (s >= 0.0_wp) .and. (s <= real(n-1, wp))
            if (.not. ok) then
                i = 1
                f = 0.0_wp
                return
            end if
        end if

        i = min(int(floor(s)) + 1, n - 1)     ! s = n-1 -> last cell, f = 1
        f = s - real(i - 1, wp)
    end subroutine

    ! --- 1D field, interpolated along `dim` (1 = along r, 2 = along z), ------
    ! --- broadcast along the other direction -----------------------------------
    subroutine interp_1d_to_plane(values, xg, dim, r, z, res, clamp)
        real(wp), intent(in)  :: values(:), xg(:), r(:), z(:)
        integer,  intent(in)  :: dim
        real(wp), intent(out) :: res(:,:)
        logical,  intent(in), optional :: clamp

        real(wp), allocatable :: q(:), line(:)
        real(wp) :: f
        integer  :: k, i
        logical  :: ok, cl

        cl = .true.
        if (present(clamp)) cl = clamp

        if (dim == 1) then
            q = r
        else
            q = z
        end if

        allocate(line(size(q)))
        do k = 1, size(q)
            call locate(xg, q(k), cl, i, f, ok)
            if (ok) then
                line(k) = (1.0_wp - f) * values(i) + f * values(i+1)
            else
                line(k) = 0.0_wp
            end if
        end do

        if (dim == 1) then
            res = spread(line, dim=2, ncopies=size(z))
        else
            res = spread(line, dim=1, ncopies=size(r))
        end if
    end subroutine

    ! --- 2D field values(nx,nz) on (xg, zg) -> plane (r, z): bilinear --------
    subroutine interp_2d_to_plane(values, xg, zg, r, z, Xs, res, clamp)
        real(wp), intent(in)  :: values(:,:), xg(:), zg(:), r(:), z(:)
        real(wp), intent(in) :: Xs
        real(wp), intent(out) :: res(:,:)
        logical,  intent(in), optional :: clamp

        integer  :: i, j, ix, jz
        integer,  allocatable :: ir(:)
        real(wp), allocatable :: fr(:), xq(:)
        logical,  allocatable :: okr(:)
        real(wp) :: fz
        logical  :: okz, cl

        allocate(xq(size(xg)))
        xq = xg-Xs
        cl = .true.
        if (present(clamp)) cl = clamp

        allocate(ir(size(r)), fr(size(r)), okr(size(r)))
        do i = 1, size(r)
            call locate(xq, r(i), cl, ir(i), fr(i), okr(i))
        end do

        do j = 1, size(z)
            call locate(zg, z(j), cl, jz, fz, okz)
            do i = 1, size(r)
                if (.not. (okz .and. okr(i))) then
                    res(i,j) = 0.0_wp
                    cycle
                end if
                ix = ir(i)
                res(i,j) = (1.0_wp - fr(i)) * ((1.0_wp - fz) * values(ix,  jz) + fz * values(ix,  jz+1)) & 
                         +          fr(i)   * ((1.0_wp - fz) * values(ix+1,jz) + fz * values(ix+1,jz+1))
            end do
        end do
    end subroutine

    ! --- 3D field values(nx,ny,nz) -> vertical plane through (Xs,Ys) --------
    ! plane point (r,z) <-> (x,y,z) = (Xs + r cos(tau), Ys + r sin(tau), z)
    ! tau in radians, measured from the +x axis towards +y.
    subroutine interp_3d_to_plane(values, xg, yg, zg, r, z, tau, Xs, Ys, res, clamp)
        real(wp), intent(in)  :: values(:,:,:), xg(:), yg(:), zg(:), r(:), z(:)
        real(wp), intent(in)  :: tau, Xs, Ys
        real(wp), intent(out) :: res(:,:)
        logical,  intent(in), optional :: clamp

        integer  :: i, j, ix, iy
        integer,  allocatable :: jz(:)
        real(wp), allocatable :: fz(:)
        logical,  allocatable :: okz(:)
        real(wp) :: ct, st, xq, yq, fx, fy
        real(wp) :: c00, c10, c01, c11, c0, c1
        logical  :: okx, oky, cl

        cl = .true.
        if (present(clamp)) cl = clamp
        ct = cosd(tau)
        st = sind(tau)

        ! z lookups do not depend on r: do them once
        allocate(jz(size(z)), fz(size(z)), okz(size(z)))
        do j = 1, size(z)
            call locate(zg, z(j), cl, jz(j), fz(j), okz(j))
        end do

        do i = 1, size(r)
            xq = Xs + r(i) * ct
            yq = Ys + r(i) * st
            call locate(xg, xq, cl, ix, fx, okx)
            call locate(yg, yq, cl, iy, fy, oky)

            if (.not. (okx .and. oky)) then
                res(i,:) = 0.0_wp
                cycle
            end if

            do j = 1, size(z)
                if (.not. okz(j)) then
                    res(i,j) = 0.0_wp
                    cycle
                end if
                ! trilinear: collapse x, then y, then z
                c00 = (1.0_wp - fx) * values(ix,iy,  jz(j)  ) + fx * values(ix+1,iy,  jz(j)  )
                c10 = (1.0_wp - fx) * values(ix,iy+1,jz(j)  ) + fx * values(ix+1,iy+1,jz(j)  )
                c01 = (1.0_wp - fx) * values(ix,iy,  jz(j)+1) + fx * values(ix+1,iy,  jz(j)+1)
                c11 = (1.0_wp - fx) * values(ix,iy+1,jz(j)+1) + fx * values(ix+1,iy+1,jz(j)+1)
                c0  = (1.0_wp - fy) * c00 + fy * c10
                c1  = (1.0_wp - fy) * c01 + fy * c11
                res(i,j) = (1.0_wp - fz(j)) * c0 + fz(j) * c1
            end do
        end do
    end subroutine

end module interp_plane_mod
