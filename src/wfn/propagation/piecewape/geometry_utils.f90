module geometry_utils
    implicit none
    integer, parameter :: dp = 8

contains

    ! Returns .true. if segment (p1,q1) and (p2,q2) intersect
    function do_intersect(p1x, p1y, q1x, q1y, p2x, p2y, q2x, q2y) result(intersect)
        real(dp), intent(in) :: p1x, p1y, q1x, q1y ! Segment 1
        real(dp), intent(in) :: p2x, p2y, q2x, q2y ! Segment 2
        logical :: intersect

        integer :: o1, o2, o3, o4

        ! Find the four orientations needed for general and special cases
        o1 = orientation(p1x, p1y, q1x, q1y, p2x, p2y)
        o2 = orientation(p1x, p1y, q1x, q1y, q2x, q2y)
        o3 = orientation(p2x, p2y, q2x, q2y, p1x, p1y)
        o4 = orientation(p2x, p2y, q2x, q2y, q1x, q1y)

        ! General case: Points are on opposite sides
        if (o1 /= o2 .and. o3 /= o4) then
            intersect = .true.
            return
        end if

        ! Special Cases: Collinear points 
        ! (Check if a point lies on the other segment)
        intersect = .false.
        if (o1 == 0 .and. on_segment(p1x, p1y, p2x, p2y, q1x, q1y)) intersect = .true.
        if (o2 == 0 .and. on_segment(p1x, p1y, q2x, q2y, q1x, q1y)) intersect = .true.
        if (o3 == 0 .and. on_segment(p2x, p2y, p1x, p1y, q2x, q2y)) intersect = .true.
        if (o4 == 0 .and. on_segment(p2x, p2y, q1x, q1y, q2x, q2y)) intersect = .true.

    end function do_intersect

    integer function orientation(px, py, qx, qy, rx, ry)
        real(dp), intent(in) :: px, py, qx, qy, rx, ry
        real(dp) :: val
        val = (qy - py) * (rx - qx) - (qx - px) * (ry - qy)
        if (abs(val) < 1e-12_dp) then
            orientation = 0  ! Collinear
        else if (val > 0) then
            orientation = 1  ! Clockwise
        else
            orientation = 2  ! Counter-clockwise
        end if
    end function orientation

    logical function on_segment(px, py, rx, ry, qx, qy)
        ! Checks if point R lies on segment PQ
        real(dp), intent(in) :: px, py, qx, qy, rx, ry
        on_segment = (rx <= max(px, qx) .and. rx >= min(px, qx) .and. &
                      ry <= max(py, qy) .and. ry >= min(py, qy))
    end function on_segment


    subroutine find_height_to_intersect(ax, ay, cx, cy, dx, dy, theta, length, possible)
        implicit none
        ! Inputs
        real, intent(in)  :: ax, ay          ! Starting point A
        real, intent(in)  :: cx, cy, dx, dy  ! Target segment C->D
        real, intent(in)  :: theta        !  absosulte Angle from point A (in radians)
        ! Outputs
        real, intent(out) :: length           ! How long the segment must be
        logical, intent(out)  :: possible     ! False if line never intersects
    
        real :: alpha,  denom
        real :: ux, uy, vx, vy, T 
    
        ! 1. Direction vector of the target line (v)
        vx = dx - cx
        vy = dy - cy
        
        ! 2. Absolute angle of the target line
        alpha = atan2(vy, vx)
        
        ! 4. Direction vector of our segment (u)
        ux = cos(theta)
        uy = sin(theta)
    
        ! 5. Solve for intersection using Cramer's rule / 2D Cross product
        ! Equation: A + L*u = C + T*v  => L*u - T*v = C - A
        ! Denominator is the cross product of direction vectors
        denom = ux * vy - uy * vx
    
        if (abs(denom) < 1e-12) then
            ! Lines are parallel
            possible = .false.
            length = 0.0
        else
            ! Length L = det([Cx-Ax, -vx], [Cy-Ay, -vy]) / denom
            length = ((cx - ax) * vy - (cy - ay) * vx) / denom
            T = (((cx - ax) * uy - (cy - ay) * ux) / denom)
            if (0 < T .and. T < 1) then 
                possible = .true.
            else
                possible = .false.
            end if 
        end if
    end subroutine

    subroutine find_first_xc(x, hknot, htx, n, z0, x0, found)
        use interpoland_1D 
        real, dimension(:), intent(in) :: x, hknot, htx
        integer, intent(in) :: n
        real, intent(in) :: z0
        real, intent(out) :: x0
        logical, intent(out) :: found

        real :: h0, z_target , x_step, x_max, x_min, x_left, x_right
        real :: slope, g_left, g_right, WINDOW
        integer :: steps, ii
        logical :: found_bracket

        WINDOW = 5
        steps = 1000

        x0 = 0.0
        h0 = interpolate(hknot, htx, x0)
        z_target = z0 + h0 
        
        x_min = max(x(1) - 40.0, -50.0)
        x_max = x(n)

        steps = 1000  ! Number of intervals for the global scan
        x_step = (x_max - x_min) / real(steps)
        found = .false.

        ! 1. Global Scan: Look for a sign change in g_func
        do ii = 0, steps - 1
            x_left = x_min + ii * x_step
            x_right = x_left + x_step
            slope = (interpolate(hknot, htx, x_left+WINDOW) - interpolate(hknot, htx,x_left))/WINDOW
            g_left = (z_target - interpolate(hknot, htx,x_left)) * slope - x_left
            
            slope = (interpolate(hknot, htx,x_right+WINDOW) - interpolate(hknot, htx, x_right))/WINDOW
            g_right = (z_target - interpolate(hknot, htx,x_right)) * slope - x_right

            if (g_left * g_right <= 0.0) then
                found = .true.
                x0 = x_left
                exit ! Found the segment containing the root
            end if
        end do



    end subroutine 

    ! subroutine find_first_xc(x, hknot, htx, n, z0, x0, found)
    !     use  interpoland_1D
    !     integer, intent(in) :: n
    !     real, intent(in) :: x(n),z0
    !     real,  dimension (:),intent(in)         :: hKnot
    !     real,  dimension (:),intent(in)         :: htx
    !     real, intent(out) :: x0
    !     logical, intent(out) :: found
    !
    !     integer :: i
    !     real, dimension(:), allocatable :: x_stencil, h_stencil
    !     real :: slope, intercept, a, b, candidate_x, z0_abs
    !     real, allocatable, dimension (:) :: x0_dummy, h0_dummy
    !
    !     found = .false.
    !     h0 = interpolate_1D(hKnot, htx, x0, h0)
    !
    !     z0_abs = h0 + z0 
    !     print*, z0
    !     print*, h0_dummy(1)
    !     print*, z0_abs
    !
    !     ! Check each segment [x_i, x_{i+1}]
    !     do i = 1, n-1
    !         x_stencil(1) = x(i)
    !         x_stencil(2) = x(i) + 5
    !         call interpolate_1D(hKnot, htx, x_stencil, h_stencil)
    !         slope = (h_stencil(2) - h_stencil(1)) / (5)
    !
    !         ! candidate_x = (slope * (z0 - h(i) + slope * x(i))) / (1.0 + slope**2)
    !         candidate_x = slope * (z0_abs - h_stencil(1)) 
    !
    !         ! Check if the calculated x falls within the current segment
    !         if (candidate_x >= x(i) .and. candidate_x <= x(i+1)) then
    !             x0 = candidate_x
    !             found = .true.
    !             return
    !         end if
    !     end do
    ! end subroutine
end module geometry_utils
