module pade_mod                                                                            
    use kind_mod, only: wp                                                                 
    implicit none                                                                          
    private                                                                                
    public :: pade_factors, qp, pade_propagator, pade_starter                                                         
                                                                                           
    integer, parameter :: qp = selected_real_kind(33)                                      
                                                                                           
contains                                                                                   


    ! Rotated Pade factors of exp(i sigma (sqrt(1+X) - 1))
    subroutine pade_propagator(sigma, m, theta, h0, A0, A1, B0, B1)
        real(qp),    intent(in)  :: sigma, theta
        integer,     intent(in)  :: m
        complex(wp), intent(out) :: h0, A0(m), A1(m), B0(m), B1(m)

        complex(qp) :: c, f(0:2*m)

        c = exp(cmplx(0.0_qp, 0.5_qp*theta, qp))
        f = series_exp_sqrt(sigma, c, 2*m)
        call pade_factors(f, m, theta, h0, A0, A1, B0, B1)
    end subroutine                                                                                        

    function series_exp_sqrt(sigma, c, N) result(g)
        ! sigma = k0*dx, c= exp(itheta/2), N=2*padde_n
        ! compute g_k for g = exp(u(Y)), using g' = u'g
        ! first compute u_k, then g_k
        real(qp),    intent(in) :: sigma
        complex(qp), intent(in) :: c
        integer,     intent(in) :: N
        complex(qp) :: g(0:N)

        complex(qp), parameter :: eye_q = (0.0_qp, 1.0_qp)
        real(qp)    :: s(0:N)
        complex(qp) :: u(0:N)
        integer     :: k, j

        ! s_k = binomial(1/2, k)
        s(0) = 1.0_qp
        do k = 1, N
            s(k) = s(k-1) * (0.5_qp - real(k-1, qp)) / real(k, qp)
        end do

        u(0) = eye_q*sigma*(c*s(0) - 1.0_qp)
        do k = 1, N
            u(k) = eye_q*sigma*c*s(k)
        end do

        ! h' = u' h  ->  k h_k = sum_{j=1..k} j u_j h_{k-j}
        g(0) = exp(u(0))
        do k = 1, N
            g(k) = (0.0_qp, 0.0_qp)
        do j = 1, k
        g(k) = g(k) + real(j, qp)*u(j)*g(k-j)
        end do
            g(k) = g(k) / real(k, qp)
        end do
    end function
   

    ! Rotated Pade factors of the starter
    !   T(X)^2 (1+X)^(-1/4) exp(i sigma (sqrt(1+X) - 1)),  T(X) = 1 - X/Xs
    subroutine pade_starter(sigma, m, theta, Xs, h0, A0, A1, B0, B1)
        real(qp),    intent(in)  :: sigma, theta, Xs
        integer,     intent(in)  :: m
        complex(wp), intent(out) :: h0, A0(m), A1(m), B0(m), B1(m)

        complex(qp) :: c, et, t0, t1, rot
        complex(qp) :: h(0:2*m), q4(0:2*m), T2(0:2), f(0:2*m)
        real(qp)    :: bin(0:2*m)
        integer     :: k, N

        N   = 2*m
        c   = exp(cmplx(0.0_qp, 0.5_qp*theta, qp))
        et  = exp(cmplx(0.0_qp,        theta, qp))
        rot = exp(cmplx(0.0_qp, -0.25_qp*theta, qp))

        h = series_exp_sqrt(sigma, c, N)

        ! (1+Y)^(-1/4): binomial(-1/4, k), with the rotation phase exp(-i theta/4)
        bin(0) = 1.0_qp
        do k = 1, N
            bin(k) = bin(k-1) * (-0.25_qp - real(k-1, qp)) / real(k, qp)
        end do
        q4 = rot * bin

        ! T(X)^2 expanded in Y
        t0 = 1.0_qp - (et - 1.0_qp)/Xs
        t1 = -et/Xs
        T2(0) = t0*t0
        T2(1) = 2.0_qp*t0*t1
        T2(2) = t1*t1

        f = series_mul(series_mul(h, q4, N), T2, N)
        call pade_factors(f, m, theta, h0, A0, A1, B0, B1)
    end subroutine
    

    ! Cauchy product of two series, truncated at order N
    function series_mul(a, b, N) result(c)
        complex(qp), intent(in) :: a(0:), b(0:)
        integer,     intent(in) :: N
        complex(qp) :: c(0:N)
        integer :: k, j

        do k = 0, N
            c(k) = (0.0_qp, 0.0_qp)
            do j = max(0, k - ubound(b,1)), min(k, ubound(a,1))
                c(k) = c(k) + a(j)*b(k-j)
            end do
        end do
    end function

    ! [n/n] Pade of Taylor series f(0:2n) in Y, as factors in X with                       
    ! Y = exp(-i theta)(1+X) - 1 :  f ~ h0 * prod (A0+A1 X)/(B0+B1 X)                      
    ! A0,A1,B0,B1 have size n; unused factors are identity (1,0).                          
    subroutine pade_factors(f, n, theta, h0, A0, A1, B0, B1)                               
        complex(qp), intent(in)  :: f(0:)                                                  
        integer,     intent(in)  :: n                                                      
        real(qp),    intent(in)  :: theta                                                  
        complex(wp), intent(out) :: h0, A0(n), A1(n), B0(n), B1(n)                         
                                                                                           
        complex(qp) :: M(n,n), q(0:n), p(0:n), e, rhs(n)                                   
        complex(qp), allocatable :: a(:), b(:)                                             
        integer :: i, j, k                                                                 
                                                                                           
        ! Pade denominator: sum_j q_j f_{k-j} = -f_k, k = n+1..2n                          
        do i = 1, n                                                                        
            k = n + i
            do j = 1, n                                                                    
                M(i,j) = f(k-j)
            end do                                                                         
            rhs(i) = -f(k)                                                                 
        end do                                                                             
        call lu_solve_q(M, rhs, n)
        q(0) = (1.0_qp, 0.0_qp)                                                            
        q(1:n) = rhs                                                                       

        ! numerator                                                                        
        do k = 0, n                                                                        
            p(k) = (0.0_qp, 0.0_qp)                                                        
            do j = 0, k                                                                    
                p(k) = p(k) + q(j)*f(k-j)                                                  
            end do                                                                         
        end do                                                                             
                                                                                           
        call factors(p, n, a)
        call factors(q, n, b)                                                              

        e = exp(cmplx(0.0_qp, -theta, qp))                                                 
        A0 = (1.0_wp, 0.0_wp);  A1 = (0.0_wp, 0.0_wp)                                      
        B0 = (1.0_wp, 0.0_wp);  B1 = (0.0_wp, 0.0_wp)                                      
        do j = 1, size(a)
            A0(j) = cmplx(1 + a(j)*(e - 1), kind=wp)                                       
            A1(j) = cmplx(a(j)*e,           kind=wp)
        end do                                                                             
        do j = 1, size(b)                                                                  
            B0(j) = cmplx(1 + b(j)*(e - 1), kind=wp)                                       
            B1(j) = cmplx(b(j)*e,           kind=wp)
        end do                                                                             
        h0 = cmplx(p(0), kind=wp)                                                          
    end subroutine
                                                                                           
    ! (1 + a Y) factors of a polynomial given by ascending coefficients c(0:n)             
    subroutine factors(c, n, a)                                                            
        complex(qp), intent(in) :: c(0:n)                                                  
        integer, intent(in) :: n                                                           
        complex(qp), allocatable, intent(out) :: a(:)                                      
        complex(qp), allocatable :: roots(:)                                               
        real(qp) :: scale                                                                  
        integer :: deg
                                                                                           
        scale = maxval(abs(c))
        deg = n                                                                            
        do while (deg > 0)                                                                 
            if (abs(c(deg)) >= 1.0e-30_qp*scale) exit                                      
            deg = deg - 1
        end do                                                                             
        if (deg == 0) then
            allocate(a(0)); return                                                         
        end if                                                                             
        allocate(roots(deg))                                                               
        call poly_roots(c(0:deg), deg, roots)
        a = -1.0_qp/roots                                                                  
    end subroutine
                                                                                           
    ! Roots of sum c(j) x^j by Laguerre + deflation + polish on the full polynomial        
    subroutine poly_roots(c, deg, roots)                                                   
        integer,     intent(in)  :: deg                                                    
        complex(qp), intent(in)  :: c(0:deg)
        complex(qp), intent(out) :: roots(deg)
        complex(qp) :: w(0:deg), z, bb, tmp                                                
        integer :: k, j                                                                    
                                                                                           
        w = c
        do k = deg, 1, -1
            z = (0.0_qp, 0.0_qp)                                                           
            call laguerre(w(0:k), k, z)                                                    
            call laguerre(c, deg, z)               ! polish                                
            roots(deg-k+1) = z
            bb = w(k)                              ! deflate
            do j = k-1, 0, -1                                                              
                tmp  = w(j)                                                                
                w(j) = bb                                                                  
                bb   = z*bb + tmp
            end do
        end do                                                                             
    end subroutine                                                                         
                                                                                           
    subroutine laguerre(a, m, x)
        complex(qp), intent(in)    :: a(0:m)                                               
        integer,     intent(in)    :: m                                                    
        complex(qp), intent(inout) :: x                                                    
        integer, parameter :: MAXIT = 400
        real(qp), parameter :: EPS = epsilon(1.0_qp)                                       
        complex(qp) :: b, d, f, g, g2, h, sq, gp, gm, dx, x1                               
        real(qp) :: err, abx                                                               
        integer :: it, j
                                                                                           
        do it = 1, MAXIT                                                                   
            b = a(m); err = abs(b); d = 0; f = 0; abx = abs(x)                             
            do j = m-1, 0, -1
                f = x*f + d                                                                
                d = x*d + b                                                                
                b = x*b + a(j)
                err = abs(b) + abx*err                                                     
            end do                                                                         
            if (abs(b) <= EPS*err) return
            g  = d/b                                                                       
            g2 = g*g                                                                       
            h  = g2 - 2*f/b
            sq = sqrt((m-1)*(m*h - g2))
            gp = g + sq                                                                    
            gm = g - sq                                                                    
            if (abs(gp) < abs(gm)) gp = gm                                                 
            if (abs(gp) > 0) then
                dx = m/gp                                                                  
            else                                                                           
                dx = (1 + abx)*exp(cmplx(0.0_qp, real(it,qp), qp))                         
            end if
            x1 = x - dx
            if (x1 == x) return
            if (mod(it, 10) /= 0) then                                                     
                x = x1
            else                                                                           
                x = x - dx*real(it/10, qp)/10     ! break limit cycles                     
            end if                                                                         
        end do
        print *, "laguerre: too many iterations"                                           
    end subroutine                                                                         
                                                                                           
    ! Dense Gaussian elimination, partial pivoting; solution overwrites rhs
    subroutine lu_solve_q(M, rhs, n)                                                       
        complex(qp), intent(inout) :: M(n,n), rhs(n)                                       
        integer, intent(in) :: n                                                           
        complex(qp) :: row(n), t
        integer :: i, j, k, piv                                                            
                                                                                           
        do k = 1, n
            piv = k - 1 + maxloc(abs(M(k:n,k)), 1)                                         
            if (piv /= k) then                                                             
                row = M(k,:); M(k,:) = M(piv,:); M(piv,:) = row
                t = rhs(k); rhs(k) = rhs(piv); rhs(piv) = t                                
            end if                                                                         
            do i = k+1, n
                t = M(i,k)/M(k,k)                                                          
                M(i,k:n) = M(i,k:n) - t*M(k,k:n)                                           
                rhs(i) = rhs(i) - t*rhs(k)
            end do                                                                         
        end do
        do i = n, 1, -1                                                                    
            do j = i+1, n
                rhs(i) = rhs(i) - M(i,j)*rhs(j)                                            
            end do
            rhs(i) = rhs(i)/M(i,i)                                                         
        end do
    end subroutine

end module pade_mod
