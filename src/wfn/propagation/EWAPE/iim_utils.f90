module iim_module
    use kind_mod, only: wp
    implicit none

    public :: iim_stencil, iim_stencil_local, reproject_previous_nodes

contains 
    subroutine iim_stencil(dist, sides, J, rho_m, rho_p, row, coef)
        use kind_mod, only: wp
        implicit none
        real(wp),    intent(in)  :: dist(3), rho_m, rho_p
        integer,     intent(in)  :: sides(3), row
        complex(wp), intent(in)  :: J(:)
        complex(wp), intent(out) :: coef(size(J), 3)

        complex(wp) :: M(5,5), B(5,3), B0(5,3)
        integer     :: ipiv(5), info, k, i
        real(wp)    :: d

        ! right-hand side N: zeros on top, identity on rows 3..5 (same for all k)
        B0 = (0.0_wp, 0.0_wp)
        do i = 1, 3
        B0(2+i, i) = (1.0_wp, 0.0_wp)
        end do

        do k = 1, size(J)
        M = (0.0_wp, 0.0_wp)

        ! jump conditions
        M(1,2) = cmplx( 1.0_wp/rho_m, 0.0_wp, wp)
        M(1,3) = cmplx(-1.0_wp/rho_p, 0.0_wp, wp)
        M(2,1) = -J(k)
        M(2,4) = (-1.0_wp, 0.0_wp)
        M(2,5) = ( 1.0_wp, 0.0_wp)

        ! Taylor expansions at the three stencil points
        do i = 1, 3
        d = dist(i)
        M(2+i,1) = (1.0_wp, 0.0_wp)
        if (sides(i) < 0) then
            M(2+i,2) = cmplx(d,           0.0_wp, wp)
            M(2+i,4) = cmplx(0.5_wp*d*d,  0.0_wp, wp)
        else
            M(2+i,3) = cmplx(d,           0.0_wp, wp)
            M(2+i,5) = cmplx(0.5_wp*d*d,  0.0_wp, wp)
        end if
        end do

        B = B0
        call zgesv(5, 3, M, 5, ipiv, B, 5, info)
        if (info /= 0) then
            print *, "iim_stencil: zgesv failed, info =", info, " k =", k
            stop
        end if
        coef(k,:) = B(row,:)
        end do
    end subroutine
    
    subroutine iim_stencil_local(dist, sides, J, rho_m, rho_p, xr, xrr, kG, S)
        
        real(wp),    intent(in)  :: dist(3), rho_m, rho_p, xr, xrr, kG
        integer,     intent(in)  :: sides(3)
        complex(wp), intent(in)  :: J
        complex(wp), intent(out) :: S(5,3)

        complex(wp) :: M(5,5), g
        real(wp)    :: s2, d
        integer     :: ipiv(5), info, i

        s2 = 1.0_wp + xr**2
        g  = cmplx(xrr, 2.0_wp*kG*xr, wp)               ! 2 i kG xr + xrr

        M = (0.0_wp, 0.0_wp)
        ! jump conditions
        M(1,1) = cmplx(0.0_wp, kG*xr*(1.0_wp/rho_p - 1.0_wp/rho_m), wp)
        M(1,2) = cmplx( s2/rho_m, 0.0_wp, wp)
        M(1,3) = cmplx(-s2/rho_p, 0.0_wp, wp)
        M(2,1) = -J
        M(2,2) =  g
        M(2,3) = -g
        M(2,4) = cmplx(-s2, 0.0_wp, wp)
        M(2,5) = cmplx( s2, 0.0_wp, wp)

        ! Taylor expansions at the three stencil points
        do i = 1, 3
            d = dist(i)
            M(2+i,1) = (1.0_wp, 0.0_wp)
            if (sides(i) < 0) then
                M(2+i,2) = cmplx(d,          0.0_wp, wp)
                M(2+i,4) = cmplx(0.5_wp*d*d, 0.0_wp, wp)
            else
                M(2+i,3) = cmplx(d,          0.0_wp, wp)
                M(2+i,5) = cmplx(0.5_wp*d*d, 0.0_wp, wp)
            end if
        end do

        ! right-hand side: zeros on top, identity on rows 3..5  ->  S directly
        S = (0.0_wp, 0.0_wp)
        do i = 1, 3
            S(2+i,i) = (1.0_wp, 0.0_wp)
        end do

        call zgesv(5, 3, M, 5, ipiv, S, 5, info)
        if (info /= 0) then
            print *, "iim_stencil_local: zgesv failed, info =", info
            stop
        end if
    end subroutine


    ! subroutine iim_stencil_local(z, js, sides, J, rho_m, rho_p,xr, xrr, kG, lo, coef)
    !     use kind_mod, only: wp
    !     use constants, only: eye
    !     implicit none
    !     real(wp),    intent(in)  :: z(:),xr, xrr, rho_m, rho_p, kG
    !     integer, intent(in) :: js(:)
    !     logical, intent(in) :: lo
    !     integer,     intent(in)  :: sides(3)
    !     complex(wp), intent(in)  :: J(:)
    !     complex(wp), intent(out) :: coef(size(J), 3)
    !
    !     real(wp)    :: dist(3)
    !     complex(wp) :: M(5,5), B(5,3), B0(5,3), s2, g2
    !     integer     :: ipiv(5), info, k, i
    !     real(wp)    :: d
    !     integer :: row 
    !
    !     ! right-hand side N: zeros on top, identity on rows 3..5 (same for all k)
    !     B0 = (0.0_wp, 0.0_wp)
    !     do i = 1, 3
    !     B0(2+i, i) = (1.0_wp, 0.0_wp)
    !     end do
    !
    !     if (lo) row = 4 else row = 5
    !     do k = 1, size(J)
    !         if (lo) then 
    !         dist = z(js(kk)-1:js(kk)+1)
    !         else 
    !         dist = z(js(kk):js(kk)+2)
    !         endif 
    !
    !         M = (0.0_wp, 0.0_wp)
    !
    !     ! jump conditions
    !     M(1,1) = eye * kG * xr(kk) * (1/rho_p - 1/rho_m)
    !     s2 = complex((1.0_wp + xr(kk)**2), 0.0_wp)
    !     M(1,2) = s2/rho_m
    !     M(1,3) = -s2/rho_p
    !     M(2,1) = -J(k)
    !     g2 = complex(xrr(k), kG*xr(kk))
    !     M(2,2) = g2
    !     M(2,3) = -g2
    !     M(2,4) = -s2
    !     M(2,5) = s2
    !
    !     ! Taylor expansions at the three stencil points
    !     do i = 1, 3
    !     d = dist(i)
    !     M(2+i,1) = (1.0_wp, 0.0_wp)
    !     if (sides(i) < 0) then
    !         M(2+i,2) = cmplx(d,           0.0_wp, wp)
    !         M(2+i,4) = cmplx(0.5_wp*d*d,  0.0_wp, wp)
    !     else
    !         M(2+i,3) = cmplx(d,           0.0_wp, wp)
    !         M(2+i,5) = cmplx(0.5_wp*d*d,  0.0_wp, wp)
    !     end if
    !     end do
    !
    !     B = B0
    !     call zgesv(5, 3, M, 5, ipiv, B, 5, info)
    !     if (info /= 0) then
    !         print *, "iim_stencil: zgesv failed, info =", info, " k =", k
    !         stop
    !     end if
    !     coef(k,:) = B(row,:)
    !     end do
    ! end subroutine


    subroutine reproject_previous_nodes(phi, S_prev, js_prev, z, xi_prev, xi_new)
        use kind_mod, only: wp
        implicit none
        complex(wp), intent(inout) :: phi(:)
        complex(wp), intent(in)    :: S_prev(5,3)
        integer,     intent(in)    :: js_prev
        real(wp),    intent(in)    :: z(:), xi_prev, xi_new

        complex(wp) :: u(5), p_minus, p_plus
        real(wp)    :: lo, hi, dz, d
        integer     :: j, j1, j2

        if (xi_new == xi_prev) return

        lo = min(xi_prev, xi_new)
        hi = max(xi_prev, xi_new)
        dz = z(2) - z(1)

        ! candidate nodes from the uniform grid, then the exact test lo < z <= hi
        j1 = max(1,         int(floor((lo - z(1))/dz)) + 1)
        j2 = min(size(z),   int(floor((hi - z(1))/dz)) + 3)

        u = matmul(S_prev, phi(js_prev-1:js_prev+1))

        do j = j1, j2
        if (.not. (z(j) > lo .and. z(j) <= hi)) cycle
        d = z(j) - xi_prev
        p_minus = u(1) + d*u(2) + 0.5_wp*d*d*u(4)      ! water-side expansion
        p_plus  = u(1) + d*u(3) + 0.5_wp*d*d*u(5)      ! air-side expansion
        if (xi_new > xi_prev) then                      ! node enters the water
            phi(j) = phi(j) + p_minus - p_plus
        else                                            ! node enters the air
            phi(j) = phi(j) + p_plus - p_minus
        end if
        end do
    end subroutine
end module 
