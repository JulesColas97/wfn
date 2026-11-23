module fast_interp
    implicit none

    ! Define the generic interface
    interface bilinear_interp
        module procedure interp_real_1d, interp_real_2d
        module procedure interp_cmplx_1d, interp_cmplx_2d
    end interface

    private :: interp_core_real, interp_core_complex

contains
    ! --- 2D WRAPPERS ---
    subroutine interp_real_2d(values, x_min, y_min, dx, dy, x, y, res,extrap)
        real, intent(in)  :: values(:,:), x(:,:), y(:,:), x_min, y_min, dx, dy
        real, intent(out) :: res(:,:)
        logical, intent(in), optional :: extrap
        call interp_2d_keys_real(values, size(values,1), size(values,2), x_min, y_min, dx, dy, &
                              x, y, size(x), res,extrap)
    end subroutine

    subroutine interp_cmplx_2d(values, x_min, y_min, dx, dy, x, y, res, extrap)
        real, intent(in)     :: x_min, y_min, dx, dy, x(:,:), y(:,:)
        complex, intent(in)  :: values(:,:)
        complex, intent(out) :: res(:,:)
        logical, intent(in), optional :: extrap
        call interp_2d_keys_complex(values, size(values,1), size(values,2), x_min, y_min, dx, dy, &
                               x, y, size(x), res, extrap)
    end subroutine

    ! --- 1D WRAPPERS ---
    subroutine interp_real_1d(values, x_min, y_min, dx, dy, x, y, res,extrap)
        real, intent(in)  :: values(:,:), x(:), y(:), x_min, y_min, dx, dy
        real, intent(out) :: res(:)
        logical, intent(in), optional :: extrap
        call interp_2d_keys_real(values, size(values,1), size(values,2), x_min, y_min, dx, dy, &
                              x, y, size(x), res,extrap)
    end subroutine

    subroutine interp_cmplx_1d(values, x_min, y_min, dx, dy, x, y, res, extrap)
        real, intent(in)     :: x_min, y_min, dx, dy, x(:), y(:)
        complex, intent(in)  :: values(:,:)
        complex, intent(out) :: res(:)
        logical, intent(in), optional :: extrap
        call interp_2d_keys_complex(values, size(values,1), size(values,2), x_min, y_min, dx, dy, &
                               x, y, size(x), res,extrap)
    end subroutine

    ! --- Version for REAL (Double Precision) ---
    subroutine interp_core_real(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res)
        integer, intent(in) :: nx, ny, n_pts
        real, intent(in) :: x_min, y_min, dx, dy
        real, intent(in) :: values(nx, ny)
        real, intent(in) :: x(n_pts), y(n_pts)
        real, intent(out) :: res(n_pts)

        integer :: i, ix, iy
        real :: x_rel, y_rel, fx, fy, inv_dx, inv_dy

        inv_dx = 1.0 / dx
        inv_dy = 1.0 / dy

        do i = 1, n_pts
            x_rel = (x(i) - x_min) * inv_dx
            y_rel = (y(i) - y_min) * inv_dy

            ix = int(x_rel) + 1
            iy = int(y_rel) + 1
            if (ix < 1 .or. ix >= nx-1 .or. iy < 1 .or. iy >= ny-1) then
                ! res(i) = 0.0
                cycle
            end if

            fx = x_rel - floor(x_rel)
            fy = y_rel - floor(y_rel)

            res(i) = (1.0 - fx) * (1.0 - fy) * values(ix, iy)     + &
                     fx * (1.0 - fy) * values(ix+1, iy)   + &
                     (1.0 - fx) * fy * values(ix, iy+1)   + &
                     fx * fy * values(ix+1, iy+1)
        end do
    end subroutine 

    ! --- Version for COMPLEX (Double Precision) ---
    subroutine interp_core_complex(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res)
        integer, intent(in) :: nx, ny, n_pts
        real, intent(in) :: x_min, y_min, dx, dy
        complex, intent(in) :: values(nx, ny)
        real, intent(in) :: x(n_pts), y(n_pts)
        complex, intent(out) :: res(n_pts)

        integer :: i, ix, iy
        real :: x_rel, y_rel, fx, fy, inv_dx, inv_dy

        inv_dx = 1.0 / dx
        inv_dy = 1.0 / dy

        do i = 1, n_pts
            x_rel = (x(i) - x_min) * inv_dx
            y_rel = (y(i) - y_min) * inv_dy

            ix = int(x_rel) + 1
            iy = int(y_rel) + 1
            if (ix < 1 .or. ix >= nx-1 .or. iy < 1 .or. iy >= ny-1) then
                ! res(i) = 0.0
                cycle
            end if

            fx = x_rel - floor(x_rel)
            fy = y_rel - floor(y_rel)

            ! Weights are real, values values are complex
            res(i) = (1.0 - fx) * (1.0 - fy) * values(ix, iy)   + &
                     fx * (1.0 - fy) * values(ix+1, iy)   + &
                     (1.0 - fx) * fy * values(ix, iy+1)   + &
                     fx * fy * values(ix+1, iy+1)
        end do
    end subroutine

    subroutine interp_2d_keys_real(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res, extrap)
    integer, intent(in) :: nx, ny, n_pts
    real, intent(in)    :: x_min, y_min, dx, dy
    real, intent(in) :: values(nx, ny)
    real, intent(in)    :: x(n_pts), y(n_pts)
    real, intent(out):: res(n_pts)
    logical, intent(in), optional :: extrap

    integer :: i, j, k, ix, iy, idx, idy
    integer :: sx(4), sy(4) ! Stencil indices
    real    :: x_rel, y_rel, tx, ty, cur_x, cur_y
    real    :: wx(4), wy(4)
    logical :: extrapolation_allowed
    
    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    ! Constants for Keys' kernel (a = -0.5)
    do i = 1, n_pts
        if (extrapolation_allowed) then 
        cur_x = max(x_min, min(x(i), x_min + (nx-1)*dx))
        cur_y = max(y_min, min(y(i), y_min + (ny-1)*dy))
        else
            cur_x = x(i)
            cur_y = y(i)
        endif 

        ! Coordinate mapping
        x_rel = (cur_x - x_min) / dx
        y_rel = (cur_y - y_min) / dy
        ix = floor(x_rel) + 1
        iy = floor(y_rel) + 1
        ! cycle if outside of bounds

        if (ix < 1 .or. ix > nx .or. iy < 1 .or. iy > ny) then
                ! res(i) = 0.0
            cycle
        endif
        
        ! Fractional distance
        tx = x_rel - floor(x_rel)
        ty = y_rel - floor(y_rel)

        ! Generate Clamped Stencil Indices
        ! This ensures that even at the very edge, we always have 4 points.
        ! If index 0, take index 1. If index nx+1, take nx.
        sx = [ix-1, ix, ix+1, ix+2]
        sy = [iy-1, iy, iy+1, iy+2]
        do j = 1, 4
            sx(j) = max(1, min(sx(j), nx))
            sy(j) = max(1, min(sy(j), ny))
        end do

        ! Calculate Keys' Weights (Catmull-Rom)
        wx(1) = -0.5*tx**3 + tx**2 - 0.5*tx
        wx(2) =  1.5*tx**3 - 2.5*tx**2 + 1.0
        wx(3) = -1.5*tx**3 + 2.0*tx**2 + 0.5*tx
        wx(4) =  0.5*tx**3 - 0.5*tx**2

        wy(1) = -0.5*ty**3 + ty**2 - 0.5*ty
        wy(2) =  1.5*ty**3 - 2.5*ty**2 + 1.0
        wy(3) = -1.5*ty**3 + 2.0*ty**2 + 0.5*ty
        wy(4) =  0.5*ty**3 - 0.5*ty**2

        ! 2D Accumulation
        res(i) = (0.0, 0.0)
        do k = 1, 4
            do j = 1, 4
                res(i) = res(i) + wx(j) * wy(k) * values(sx(j), sy(k))
            end do
        end do
    end do
    end subroutine

    subroutine interp_2d_keys_complex(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res,extrap)
    integer, intent(in) :: nx, ny, n_pts
    real, intent(in)    :: x_min, y_min, dx, dy
    complex, intent(in) :: values(nx, ny)
    real, intent(in)    :: x(n_pts), y(n_pts)
    complex, intent(out):: res(n_pts)
    logical, intent(in), optional :: extrap

    integer :: i, j, k, ix, iy, idx, idy
    integer :: sx(4), sy(4) ! Stencil indices
    real    :: x_rel, y_rel, tx, ty, cur_x, cur_y
    real    :: wx(4), wy(4)
    logical :: extrapolation_allowed
    
    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if
    
    ! Constants for Keys' kernel (a = -0.5)
    do i = 1, n_pts

        if (extrapolation_allowed) then 
        cur_x = max(x_min, min(x(i), x_min + (nx-1)*dx))
        cur_y = max(y_min, min(y(i), y_min + (ny-1)*dy))
        else
            cur_x = x(i)
            cur_y = y(i)
        endif 
        ! Coordinate mapping
        x_rel = (cur_x - x_min) / dx
        y_rel = (cur_y - y_min) / dy
        ix = floor(x_rel) + 1
        iy = floor(y_rel) + 1
        ! cycle if outside of bounds
        if (ix < 1 .or. ix > nx .or. iy < 1 .or. iy > ny) then
            ! res(i) = 0.0
            cycle
        end if
        
        ! Fractional distance
        tx = x_rel - floor(x_rel)
        ty = y_rel - floor(y_rel)

        ! Generate Clamped Stencil Indices
        ! This ensures that even at the very edge, we always have 4 points.
        ! If index 0, take index 1. If index nx+1, take nx.
        sx = [ix-1, ix, ix+1, ix+2]
        sy = [iy-1, iy, iy+1, iy+2]
        do j = 1, 4
            sx(j) = max(1, min(sx(j), nx))
            sy(j) = max(1, min(sy(j), ny))
        end do

        ! Calculate Keys' Weights (Catmull-Rom)
        wx(1) = -0.5*tx**3 + tx**2 - 0.5*tx
        wx(2) =  1.5*tx**3 - 2.5*tx**2 + 1.0
        wx(3) = -1.5*tx**3 + 2.0*tx**2 + 0.5*tx
        wx(4) =  0.5*tx**3 - 0.5*tx**2

        wy(1) = -0.5*ty**3 + ty**2 - 0.5*ty
        wy(2) =  1.5*ty**3 - 2.5*ty**2 + 1.0
        wy(3) = -1.5*ty**3 + 2.0*ty**2 + 0.5*ty
        wy(4) =  0.5*ty**3 - 0.5*ty**2

        ! 2D Accumulation
        res(i) = (0.0, 0.0)
        do k = 1, 4
            do j = 1, 4
                res(i) = res(i) + wx(j) * wy(k) * values(sx(j), sy(k))
            end do
        end do
    end do
    end subroutine

    subroutine interp_2d_cubic_extrap(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res)
    integer, intent(in) :: nx, ny, n_pts
    real, intent(in)    :: x_min, y_min, dx, dy
    complex, intent(in) :: values(nx, ny)
    real, intent(in)    :: x(n_pts), y(n_pts)
    complex, intent(out):: res(n_pts)

    integer :: i, j, k, i_start, j_start, ic, jc
    real    :: x_rel, y_rel, fx, fy, cur_x, cur_y
    real    :: wx(4), wy(4)
    real, parameter :: eps = 1.0e-12
    complex :: local_nodes(2,2)
    real    :: r_min, r_max, i_min, i_max

    do i = 1, n_pts
        ! 1. Fix numerical aberrations (-8e-16)
        ! cur_x = max(x_min, min(x(i), x_min + (nx-1)*dx))
        ! cur_y = max(y_min, min(y(i), y_min + (ny-1)*dy))

        ! 2. Coordinate mapping (0 to nx-1)
        x_rel = (x(i) - x_min) / dx
        y_rel = (y(i) - y_min) / dy

        ! 3. Shifted Window Selection
        ! We pick the start of a 4-point block. 
        ! We want the target to be in the middle (between index 2 and 3), 
        ! but if we are at the edge, we clamp the window, not the indices.
        i_start = floor(x_rel)  ! Base index (0-based)
        
        j_start = floor(y_rel)

        if (i_start < 0 .or. i_start >= nx-2 .or. j_start < 0 .or. j_start >= ny-2) then
            ! res(i) = 0.0
            cycle
        end if
        i_start = max(0, min(i_start - 1, nx - 4)) 
        j_start = max(0, min(j_start - 1, ny - 4))
        ! 4. Local coordinates within the 4-point window
        ! fx and fy will range from 0.0 to 3.0
        fx = x_rel - real(i_start)
        fy = y_rel - real(j_start)

        ! 5. Lagrange Cubic Weights
        ! These weights represent the cubic polynomial trend of the 4 points
        ! If fx=0, it picks point 1. If fx=1, point 2. If fx=0.5, it 'extrapolates' between 1 & 2.
        wx(1) = -(fx - 1.0) * (fx - 2.0) * (fx - 3.0) / 6.0
        wx(2) =  fx * (fx - 2.0) * (fx - 3.0) / 2.0
        wx(3) = -fx * (fx - 1.0) * (fx - 3.0) / 2.0
        wx(4) =  fx * (fx - 1.0) * (fx - 2.0) / 6.0

        wy(1) = -(fy - 1.0) * (fy - 2.0) * (fy - 3.0) / 6.0
        wy(2) =  fy * (fy - 2.0) * (fy - 3.0) / 2.0
        wy(3) = -fy * (fy - 1.0) * (fy - 3.0) / 2.0
        wy(4) =  fy * (fy - 1.0) * (fy - 2.0) / 6.0

        ! 6. 2D Accumulation (1-based indexing for values)
        res(i) = (0.0, 0.0)
        do k = 1, 4
            do j = 1, 4
                res(i) = res(i) + wx(j) * wy(k) * values(i_start + j, j_start + k)
            end do
        end do
        ! 5. MONOTONICITY LIMITER
        ! Identify the 4 immediate neighbors (the 2x2 cell containing the point)
        ic = max(1, min(int(floor(x_rel)) + 1, nx - 1))
        jc = max(1, min(int(floor(y_rel)) + 1, ny - 1))

        local_nodes(1,1) = values(ic, jc)
        local_nodes(2,1) = values(ic+1, jc)
        local_nodes(1,2) = values(ic, jc+1)
        local_nodes(2,2) = values(ic+1, jc+1)

        ! Limit Real Part
        r_min = minval(real(local_nodes))
        r_max = maxval(real(local_nodes))
        res(i) = cmplx(max(r_min, min(r_max, real(res(i)))), aimag(res(i)))

        ! Limit Imaginary Part
        i_min = minval(aimag(local_nodes))
        i_max = maxval(aimag(local_nodes))
        res(i) = cmplx(real(res(i)), max(i_min, min(i_max, aimag(res(i)))))
    end do
    end subroutine

    subroutine interp_2d_bspline_decentered(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res)
    integer, intent(in) :: nx, ny, n_pts
    real, intent(in)    :: x_min, y_min, dx, dy
    complex, intent(in) :: values(nx, ny)
    real, intent(in)    :: x(n_pts), y(n_pts)
    complex, intent(out):: res(n_pts)

    integer :: i, j, k, i_start, j_start, i_cell, j_cell
    real    :: x_rel, y_rel, fx, fy, cur_x, cur_y
    real    :: wx(4), wy(4)
    real    :: r_min, r_max, i_min, i_max
    complex :: local_2x2(2,2)
    real, parameter :: one_sixth = 1.0 / 6.0

    do i = 1, n_pts
        x_rel = (x(i) - x_min) / dx
        y_rel = (y(i) - y_min) / dy

        ! 2. Identify Current Cell (for the Limiter)
        i_cell = int(floor(x_rel)) + 1
        j_cell = int(floor(y_rel)) + 1
        if (i_cell < 1 .or. i_cell >= nx-1 .or. j_cell < 1 .or. j_cell >= ny-1) then
            ! res(i) = 0.0
            cycle
        end if

        ! 3. Shifted Window Selection (1 to nx-3)
        ! This allows the stencil to stay fixed at the edges for extrapolation
        i_start = max(1, min(i_cell - 1, nx - 3))
        j_start = max(1, min(j_cell - 1, ny - 3))

        ! 4. Local coordinates relative to the START of the 4-point window
        ! t ranges from [0.0, 3.0]
        fx = x_rel - real(i_start - 1)
        fy = y_rel - real(j_start - 1)

        ! 5. Cubic B-Spline Basis Evaluation
        ! We evaluate the 4 basis functions based on the distance t
        wx = get_bspline_weights(fx)
        wy = get_bspline_weights(fy)

        ! 6. Accumulate Result
        res(i) = (0.0, 0.0)
        do k = 1, 4
            do j = 1, 4
                res(i) = res(i) + wx(j) * wy(k) * values(i_start + j - 1, j_start + k - 1)
            end do
        end do

        ! 7. MONOTONICITY LIMITER (Applied to the 2x2 physical cell)
        ! local_2x2(1,1) = values(i_cell,   j_cell)
        ! local_2x2(2,1) = values(i_cell+1, j_cell)
        ! local_2x2(1,2) = values(i_cell,   j_cell+1)
        ! local_2x2(2,2) = values(i_cell+1, j_cell+1)
        !
        ! r_min = minval(real(local_2x2)); r_max = maxval(real(local_2x2))
        ! i_min = minval(aimag(local_2x2)); i_max = maxval(aimag(local_2x2))
        !
        ! res(i) = cmplx(max(r_min, min(r_max, real(res(i)))), &
        !                max(i_min, min(i_max, aimag(res(i)))))
    end do

    contains

    pure function get_bspline_weights(t) result(w)
        real, intent(in) :: t
        real :: w(4)
        ! Segment-based weight calculation for t in [0, 3]
        ! B0, B1, B2, B3 are the four basis functions overlapping this window
        w(1) = b_kern(t)
        w(2) = b_kern(t - 1.0)
        w(3) = b_kern(t - 2.0)
        w(4) = b_kern(t - 3.0)
    end function

    pure real function b_kern(xx)
        real, intent(in) :: xx
        real :: ax
        ax = abs(xx - 1.0) ! Kernel is centered at relative index 1 (the 2nd node)
        if (ax < 1.0) then
            b_kern = one_sixth * (3.0*ax**3 - 6.0*ax**2 + 4.0)
        else if (ax < 2.0) then
            b_kern = one_sixth * (2.0 - ax)**3
        else
            b_kern = 0.0
        end if
    end function

    end subroutine




    ! bi cubic interpolation
!     subroutine interp_core_complex(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res)
!     integer, intent(in) :: nx, ny, n_pts
!     real, intent(in) :: x_min, y_min, dx, dy
!     complex, intent(in) :: values(nx, ny)
!     real, intent(in) :: x(n_pts), y(n_pts)
!     complex, intent(out) :: res(n_pts)
!
!     integer :: i, j, k, ix, iy
!     real :: x_rel, y_rel, fx, fy
!     real :: wx(4), wy(4) ! 4x4 stencil for smoother results
!
!     do i = 1, n_pts
!         x_rel = (x(i) - x_min) / dx
!         y_rel = (y(i) - y_min) / dy
!
!         ! Get the "left-most" index of the 4-point block
!         ix = floor(x_rel) 
!         iy = floor(y_rel)
!
!         ! Check boundaries for 4x4 (needs nodes: ix-1, ix, ix+1, ix+2)
!         if (ix < 2 .or. ix > nx-2 .or. iy < 2 .or. iy > ny-2) then
!             cycle
!         end if
!
!         ! Distance from the node 'ix'
!         fx = x_rel - real(ix)
!         fy = y_rel - real(iy)
!
!         ! Catmull-Rom Cubic Weights (Smoother than Lagrange)
!         ! These reduce the "ringing" (oscillations) significantly
!         wx(1) = -0.5*fx + fx**2 - 0.5*fx**3
!         wx(2) = 1.0 - 2.5*fx**2 + 1.5*fx**3
!         wx(3) = 0.5*fx + 2.0*fx**2 - 1.5*fx**3
!         wx(4) = -0.5*fx**2 + 0.5*fx**3
!
!         wy(1) = -0.5*fy + fy**2 - 0.5*fy**3
!         wy(2) = 1.0 - 2.5*fy**2 + 1.5*fy**3
!         wy(3) = 0.5*fy + 2.0*fy**2 - 1.5*fy**3
!         wy(4) = -0.5*fy**2 + 0.5*fy**3
!
!         res(i) = 0.0
!         do k = 1, 4
!             do j = 1, 4
!                 res(i) = res(i) + wx(j) * wy(k) * values(ix + j - 2, iy + k - 2)
!             end do
!         end do
!     end do
! end subroutine
!
! subroutine interp_core_complex(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res)
!     integer, intent(in) :: nx, ny, n_pts
!     real, intent(in) :: x_min, y_min, dx, dy
!     complex, intent(in) :: values(nx, ny)
!     real, intent(in) :: x(n_pts), y(n_pts)
!     complex, intent(out) :: res(n_pts)
!
!     integer :: i, j, k, ix, iy
!     real :: x_rel, y_rel, fx, fy
!     real :: wx(4), wy(4)
!     real, parameter :: one_sixth = 1.0 / 6.0
!
!     do i = 1, n_pts
!         x_rel = (x(i) - x_min) / dx
!         y_rel = (y(i) - y_min) / dy
!
!         ! ix is the index of the grid point to the immediate left of the target x
!         ix = floor(x_rel) + 1
!         iy = floor(y_rel) + 1
!
!         ! Boundary check for a 4x4 stencil: 
!         ! We need to access (ix-1) through (ix+2).
!         ! So ix-1 >= 1  =>  ix >= 2
!         ! And ix+2 <= nx =>  ix <= nx-2
!         if (ix < 2 .or. ix > nx - 2 .or. iy < 2 .or. iy > ny - 2) then
!             ! res(i) = (0.0, 0.0)
!             cycle
!         end if
!
!         ! fx is the distance from the node 'ix' (ranges 0 to 1)
!         fx = x_rel - floor(x_rel)
!         fy = y_rel - floor(y_rel)
!
!         ! Cubic B-Spline Weights
!         wx(1) = one_sixth * (1.0 - fx)**3
!         wx(2) = one_sixth * (3.0*fx**3 - 6.0*fx**2 + 4.0)
!         wx(3) = one_sixth * (-3.0*fx**3 + 3.0*fx**2 + 3.0*fx + 1.0)
!         wx(4) = one_sixth * fx**3
!
!         wy(1) = one_sixth * (1.0 - fy)**3
!         wy(2) = one_sixth * (3.0*fy**3 - 6.0*fy**2 + 4.0)
!         wy(3) = one_sixth * (-3.0*fy**3 + 3.0*fy**2 + 3.0*fy + 1.0)
!         wy(4) = one_sixth * fy**3
!
!         ! Accumulate 4x4 sum using the corrected offset
!         res(i) = (0.0, 0.0)
!         do k = 1, 4
!             do j = 1, 4
!                 ! Accessing: (ix-1, iy-1), (ix, iy), (ix+1, iy+1), (ix+2, iy+2)
!                 res(i) = res(i) + wx(j) * wy(k) * values(ix + j - 2, iy + k - 2)
!             end do
!         end do
!     end do
! end subroutine
!
!
!     subroutine interp_core_complex(values, nx, ny, x_min, y_min, dx, dy, x, y, n_pts, res)
!     integer, intent(in) :: nx, ny, n_pts
!     real, intent(in) :: x_min, y_min, dx, dy
!     complex, intent(in) :: values(nx, ny)
!     real, intent(in) :: x(n_pts), y(n_pts)
!     complex, intent(out) :: res(n_pts)
!
!     integer :: i, j, k, ix, iy
!     real :: x_rel, y_rel, fx, fy
!     real :: wx(3), wy(3)  ! Weights for the 3x3 stencil
!
!     do i = 1, n_pts
!         x_rel = (x(i) - x_min) / dx
!         y_rel = (y(i) - y_min) / dy
!
!         ! Get the central index (rounding to nearest)
!         ix = nint(x_rel) + 1
!         iy = nint(y_rel) + 1
!
!         ! Check boundaries for a 3x3 stencil (needs 1 node buffer on all sides)
!         if (ix < 2 .or. ix > nx-1 .or. iy < 2 .or. iy > ny-1) then
!             ! res(i) = 0.0 ! Or handle extrapolation
!             cycle
!         end if
!
!         ! Distance from the central node ix, iy
!         fx = x_rel - (ix - 1)
!         fy = y_rel - (iy - 1)
!
!         ! Quadratic Lagrange weights for X direction
!         wx(1) = 0.5 * fx * (fx - 1.0)             ! Weight for ix-1
!         wx(2) = (1.0 - fx**2)                    ! Weight for ix
!         wx(3) = 0.5 * fx * (fx + 1.0)             ! Weight for ix+1
!
!         ! Quadratic Lagrange weights for Y direction
!         wy(1) = 0.5 * fy * (fy - 1.0)             ! Weight for iy-1
!         wy(2) = (1.0 - fy**2)                    ! Weight for iy
!         wy(3) = 0.5 * fy * (fy + 1.0)             ! Weight for iy+1
!
!         ! Accumulate the 3x3 weighted sum
!         res(i) = 0.0
!         do k = 1, 3
!             do j = 1, 3
!                 res(i) = res(i) + wx(j) * wy(k) * values(ix + j - 2, iy + k - 2)
!             end do
!         end do
!     end do
! end subroutine

end module fast_interp
