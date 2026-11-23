module trisolve_mod
    use kind_mod, only: wp
    implicit none
    private :: tri_step_fused
    public :: tri_multiply, thomas_solve, tri_step, apply_rational, apply_rational_new

    ! Tridiagonal storage convention (all arrays of size n, 1-based):
    !   lo(i) multiplies phi(i-1)   (lo(1) unused)
    !   di(i) multiplies phi(i)
    !   up(i) multiplies phi(i+1)   (up(n) unused)

contains

    ! =========================================================================
    ! MULTIPLICATION : rhs = B * phi
    ! =========================================================================
    subroutine tri_multiply(lb, db, ub, phi, rhs)
        complex(wp), intent(in)  :: lb(:), db(:), ub(:), phi(:)
        complex(wp), intent(out) :: rhs(:)
        integer :: n, i

        n = size(phi)
        rhs(1) = db(1)*phi(1) + ub(1)*phi(2)
        do i = 2, n-1
            rhs(i) = lb(i)*phi(i-1) + db(i)*phi(i) + ub(i)*phi(i+1)
        end do
        rhs(n) = lb(n)*phi(n-1) + db(n)*phi(n)
    end subroutine

    ! =========================================================================
    ! INVERSION : solve A * x = rhs with the Thomas algorithm
    !   cp is a work array of size n (modified coefficients of the sweep).
    !   x may be the same array as rhs' target only if it is a separate
    !   argument: rhs is left untouched here, x is the solution.
    ! =========================================================================
    subroutine thomas_solve(la, da, ua, rhs, x, cp)
        complex(wp), intent(in)  :: la(:), da(:), ua(:), rhs(:)
        complex(wp), intent(out) :: x(:)
        complex(wp), intent(out) :: cp(:)
        complex(wp) :: m
        integer :: n, i

        n = size(rhs)

        ! forward elimination
        cp(1) = ua(1)/da(1)
        x(1)  = rhs(1)/da(1)
        do i = 2, n
            m = da(i) - la(i)*cp(i-1)
            if (i < n) cp(i) = ua(i)/m
            x(i) = (rhs(i) - la(i)*x(i-1))/m
        end do

        ! back substitution
        do i = n-1, 1, -1
            x(i) = x(i) - cp(i)*x(i+1)
        end do
    end subroutine

    ! =========================================================================
    ! One Pade factor, in place:  A phi_new = B phi
    ! =========================================================================
    subroutine tri_step(phi, la, da, ua, lb, db, ub, rhs, cp)
        complex(wp), intent(inout) :: phi(:)
        complex(wp), intent(in)    :: la(:), da(:), ua(:), lb(:), db(:), ub(:)
        complex(wp), intent(inout) :: rhs(:), cp(:)

        ! multiplication
        call tri_multiply(lb, db, ub, phi, rhs)
        ! inversion
        call thomas_solve(la, da, ua, rhs, phi, cp)
    end subroutine


    subroutine tri_step_fused(phi, sub, dia, sup, a0, a1, b0, b1, rhs, cp)
        complex(wp), intent(inout), contiguous :: phi(:)
        complex(wp), intent(in),    contiguous :: sub(:), dia(:), sup(:)
        complex(wp), intent(in)    :: a0, a1, b0, b1
        complex(wp), intent(inout), contiguous :: rhs(:), cp(:)

        complex(wp) :: r, d, la, m, inv
        integer :: n, i

        n = size(phi)

        ! ---- forward sweep fused with the multiplication ---------------------
        ! row 1
        r   = (a0 + a1*dia(1))*phi(1) + a1*sup(1)*phi(2)
        inv = 1.0_wp / (b0 + b1*dia(1))
        cp(1)  = b1*sup(1)*inv
        rhs(1) = r*inv

        ! rows 2..n-1
        do i = 2, n-1
            r   = (a0 + a1*dia(i))*phi(i) + a1*(sub(i)*phi(i-1) + sup(i)*phi(i+1))
            la  = b1*sub(i)
            d   = b0 + b1*dia(i)
            m   = d - la*cp(i-1)
            inv = 1.0_wp / m
            cp(i)  = b1*sup(i)*inv
            rhs(i) = (r - la*rhs(i-1))*inv
        end do

        ! row n
        r   = (a0 + a1*dia(n))*phi(n) + a1*sub(n)*phi(n-1)
        la  = b1*sub(n)
        m   = (b0 + b1*dia(n)) - la*cp(n-1)
        rhs(n) = (r - la*rhs(n-1)) / m

        ! ---- back substitution (phi overwritten here) ------------------------
        phi(n) = rhs(n)
        do i = n-1, 1, -1
            phi(i) = rhs(i) - cp(i)*phi(i+1)
        end do
    end subroutine

    ! =========================================================================
    ! phi <- h0 * prod_j (B0_j + B1_j X)^-1 (A0_j + A1_j X) phi
    !   X given by (sub, dia, sup)
    ! =========================================================================
    subroutine apply_rational(phi, h0, A0, A1, B0, B1, sub, dia, sup, rhs, cp)
        complex(wp), intent(inout) :: phi(:), rhs(:), cp(:)
        complex(wp), intent(in)    :: h0, A0(:), A1(:), B0(:), B1(:)
        complex(wp), intent(in)    :: sub(:), dia(:), sup(:)
        integer :: j

        ! do j = 1, size(A0)
        !     call tri_step(phi,                                    &
        !                   B1(j)*sub, B0(j) + B1(j)*dia, B1(j)*sup, &   ! A (inverted)
        !                   A1(j)*sub, A0(j) + A1(j)*dia, A1(j)*sup, &   ! B (multiplied)
        !                   rhs, cp)
        ! end do

        do j = 1, size(A0)
            call tri_step_fused(phi, sub, dia, sup, A0(j), A1(j), B0(j), B1(j), rhs, cp)
        end do
        phi = h0*phi
    end subroutine

    subroutine apply_rational_new(phi, h0, A0, A1, B0, B1, sub, dia, sup, rhs, cp)
        complex(wp), intent(inout) :: phi(:), rhs(:), cp(:)
        complex(wp), intent(in)    :: h0, A0(:), A1(:), B0(:), B1(:)
        complex(wp), intent(in)    :: sub(:), dia(:), sup(:)
        complex(wp), allocatable :: la(:), da(:), ua(:), lb(:), db(:), ub(:)
        complex(wp) :: m
        integer :: j,k, n
        n = size(sub)
        allocate(la(n))
        allocate(da(n))
        allocate(ua(n))
        allocate(lb(n))
        allocate(db(n))
        allocate(ub(n))

        do j = 1, size(A0)

            la = B1(j) * sub 
            da = B0(j) + B1(j) * dia
            ua = B1(j) * sup 

            lb = A1(j) * sub 
            db = A0(j) + A1(j) * dia
            ub = A1(j) * sup 
            
            ! multiply 
            rhs(1) = db(1)*phi(1) +ub(1)*phi(2)
            do  k=2,n-1
                rhs(k) = lb(k)*phi(k-1) +db(k)*phi(k) +ub(k)*phi(k+1)
            enddo
            rhs(n) = lb(n)*phi(n-1) +db(n)*phi(n)

            ! solve
            cp(1) = ua(1)/da(1)
            phi(1) = rhs(1)/da(1)
            do k=2,n
                m = 1/( da(k) -la(k)*cp(k-1))
                if (k<n) then 
                    cp(k) = ua(k) * m
                endif
                phi(k) = (rhs(k) -la(k)*phi(k-1)) * m
            enddo
            do k=n-1, 1, -1
                phi(k) = phi(k) - cp(k)*phi(k+1)
            enddo 
        end do

        phi = h0*phi
    end subroutine
end module trisolve_mod
