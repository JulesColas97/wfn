module simple_wape
implicit none
public :: propa_1angle_1freq_1domain
contains

subroutine propa_1angle_1freq_1domain(phi0, phi, p_over_p0, f, k0, beta, c0,  &
        pe_model,  h10, h12, h20, h22, mach, step,  pml, dx, dz, x_grid, z_grid, Zs)
    use types
    use utils
    use interpo
    use read_write
    use hdf5
    use splines
    implicit none
    ! input parameters
    ! type(complex),dimension(:), intent(in)                  :: phi_m1
    type(complex),dimension(:), intent(in)                  :: phi0
    type(complex),dimension(:,:), intent(inout)             :: phi
    type(real),dimension(:,:), intent(inout)                :: p_over_p0
    type(real), intent(in)                                  :: f, k0
    type(complex), intent(in)                               :: beta
    type(real), intent(in)                                  :: c0
    integer, intent(in)                                     :: pe_model
    type(real),dimension(:,:), intent(in)                   :: h10,h12,h20,h22,mach
    type(integer), intent(in)                               :: step
    type(perfectly_matched_layer),intent(inout)             :: pml
    type(real), intent(in)                                  :: dx, dz, Zs
    type(real),dimension(:,:), intent(in)                   :: x_grid, z_grid


    ! local variables
    type(complex),allocatable, dimension(:)                 :: u
    type(solve)                                             :: m
    real                                                    :: omega
    complex                                                 :: kappa
    complex,allocatable,dimension(:)                        :: AlphamBeta, AlphapBeta


    ! generel variables
    real                                                    :: pi = 3.141592653589793
    complex                                                 :: eye = (0.0,1.0)
    integer                                                 :: ii, kk, index
    integer                                                 :: nx, nz


    ! time computation 
    real                                                    :: start, finish, t_tot
    real                                                    :: start_propa, finish_propa, t_propa
    real                                                    :: start_dl, finish_dl, t_dl
    real                                                    :: start_interp, finish_interp, t_interp
    real                                                    :: start_interpr, finish_interpr, t_interpr
    real                                                    :: start_interps, finish_interps, t_interps
    real                                                    :: start_interpp, finish_interpp, t_interpp


    nx = size(phi,1)
    nz = size(phi,2)
    ! init frequency param
    omega = 2*pi*f
    ! k0 = 2*pi*f/c0

    !totally reflective ground hard code
    ! PML definition
    !-------------------------------------------------------------------------
    allocate(pml%sigma(nz))
    allocate(pml%sigma_de(nz))
    pml%sigma = 0.
    pml%sigma_de = 0.

    do ii = nz - pml%size, nz
        pml%sigma(ii) = pml%param * (real(ii- nz + pml%size)/real(pml%size))**pml%n
        pml%sigma_de(ii) = pml%param*pml%n*real((ii- nz + pml%size))**(pml%n-1)/(dz*(pml%size**pml%n))
    enddo


    ! init vector for matrix resolution
    !-------------------------------------------------------------------------
    allocate(m%h10(nz))
    allocate(m%h12(nz))
    allocate(m%h20(nz))
    allocate(m%h22(nz))

    kappa = eye*k0*dx/2

    m%h10 = h10(1,1:nz)
    m%h12 = h12(1,1:nz)
    m%h20 = h20(1,1:nz)
    m%h22 = h22(1,1:nz)
    
    allocate(m%a_vect(nz))
    allocate(m%b_vect(nz))
    allocate(m%c_vect(nz))

    allocate(m%d_vect(nz))
    allocate(m%e_vect(nz))
    allocate(m%f_vect(nz))
    allocate(AlphamBeta(nz))
    allocate(AlphapBeta(nz))

    AlphapBeta = ((1/(((k0+eye*pml%sigma/c0)**2)*dz**2)) - &
                & (eye*pml%sigma_de)/(2*(omega+eye*pml%sigma)*((k0+eye*pml%sigma/c0)**2)*dz))
    AlphamBeta = ((1/(((k0+eye*pml%sigma/c0)**2)*dz**2)) + &
                & (eye*pml%sigma_de)/(2*(omega+eye*pml%sigma)*((k0+eye*pml%sigma/c0)**2)*dz))

    m%a_vect = (m%h10 - kappa*m%h20) - &
                2*(m%h12- kappa*m%h22)/(((k0+eye*pml%sigma/c0)**2)*dz**2)
    m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
    m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


    m%d_vect = (m%h10 + kappa*m%h20) - &
                2*(m%h12 + kappa*m%h22)/(((k0+eye*pml%sigma/c0)**2)*dz**2)
    m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
    m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

    m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
    m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

    m%a_vect(1) = m%a_vect(1) + 2*eye*k0*beta*dz*m%b_vect(1);
    m%d_vect(1) = m%d_vect(1) + 2*eye*k0*beta*dz*m%e_vect(1);


    ! propagation
    !----------------------------------------------------------------------------------------
    ! definition de phi(xn), phi(xn-1) and phi (xn+1) for the calculation of p (cf ostashev)

    ! allocate(phi_m1(nz))
    ! allocate(phi_p1(nz))
    allocate(u(nz))

    phi(1,:) = phi0
    ! phi_p1 = phi(1,:)
    p_over_p0(1,:) = max(abs(sqrt(x_grid(1,:)**2+(z_grid(1,:)-Zs)**2)*phi(1,:) &
                    /sqrt(max(x_grid(1,:),1e-12))*exp(eye*k0*x_grid(1,:))), 1e-10)

    ! x advancement
    !---------------------------------------------------------------------------
    ! print*,"begin x loop ..."

    t_propa = 0 
    t_dl = 0 
    t_interps = 0
    t_interp = 0
    t_interpr = 0
    t_interpp = 0
    CALL CPU_TIME(start)
    do kk = 2, nx
        index = (kk-1)/step+1
        ! update matrix coefficient with flmow grid 
        if(mod(kk-1,step) == 0) then
            m%h10 = h10(index,1:nz)
            m%h12 = h12(index,1:nz)
            m%h20 = h20(index,1:nz)
            m%h22 = h22(index,1:nz)

            m%a_vect = (m%h10 - kappa*m%h20) - &
                        & 2*(m%h12-kappa*m%h22)/(((k0+eye*pml%sigma/c0)**2)*dz**2)
            m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
            m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


            m%d_vect = (m%h10 + kappa*m%h20) - &
                        & 2*(m%h12+kappa*m%h22)/(((k0+eye*pml%sigma/c0)**2)*dz**2)
            m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
            m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

            m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
            m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

            m%a_vect(1) = m%a_vect(1) + 2*eye*k0*beta*dz*m%b_vect(1);
            m%d_vect(1) = m%d_vect(1) + 2*eye*k0*beta*dz*m%e_vect(1);
        endif

        

        ! print*, '--------------------------------------------------------------'
        ! print*, m%h10(1),m%h12(1),m%h20(1),m%h22(1)
        ! print*, m%d_vect(10), m%e_vect(10), m%f_vect(10)
        ! print*, '--------------------------------------------------------------'

        ! phi_p1 = phi(kk-1,:)
        ! solve system 
        CALL CPU_TIME(start_propa)
        u = phi(kk-1,:)
        phi(kk,1) = m%d_vect(1)*u(1) + m%f_vect(1)*u(2)
        phi(kk,2:nz-1) = m%e_vect(2:nz-1)*u(1:nz-2) + m%d_vect(2:nz-1)*u(2:nz-1) &
            & + m%f_vect(2:nz-1)*u(3:nz)
        phi(kk,nz) = m%e_vect(nz)*u(nz-1) + m%d_vect(nz)*u(nz)
        
        call res_trid(phi(kk,1:nz), m%b_vect(2:nz), m%a_vect, m%c_vect(1:nz-1),nz)

        phi(kk,nz) = 0;
        CALL CPU_TIME(finish_propa)
        t_propa = t_propa + finish_propa - start_propa
        ! Compute delta L
        CALL CPU_TIME(start_dl)
        if ((pe_model == 2) .or. (pe_model == 4)) then 
            do ii=1,nz
                ! p_over_p0(kk,ii) = max(abs(sqrt(x_grid(kk,ii)**2+(z_grid(kk,ii)-Zs)**2)*phi(kk,ii) &
                !     /sqrt(x_grid(kk,ii))*exp(eye*k0*x_grid(kk,ii))), 1e-10)
                p_over_p0(kk,ii) = max(abs(sqrt(x_grid(kk,ii)**2+(z_grid(kk,ii)-Zs)**2)*phi(kk,ii) &
                    /sqrt(x_grid(kk,ii))), 1e-10)
            enddo
        else
       
        if (kk>2) then 
            do ii=1,nz
                ! p_over_p0(kk-1,ii) = abs(sqrt(x_grid(kk-1,ii) +  &
                !     (z_grid(kk-1,ii)-Zs)**2/x_grid(kk-1,ii)) * ((1 - mach(index,ii))*phi(kk-1,ii) + &
                !     (eye*mach(index,ii)/(2*k0*dx))*(phi(kk,ii)-phi(kk-2,ii))))
                ! p_over_p0(kk-1,ii) =  sqrt(sqrt((x_grid(kk-1,ii)**2 +(z_grid(kk-1,ii)-Zs)**2))) * &
                p_over_p0(kk-1,ii) =  sqrt((x_grid(kk-1,ii)**2 +(z_grid(kk-1,ii)-Zs)**2))/sqrt(x_grid(kk-1,ii)) * &
                    abs((1 - mach(index,ii))*phi(kk-1,ii) + &
                    (eye*mach(index,ii)/(2*k0*dx))*(phi(kk,ii)-phi(kk-2,ii)))
            enddo
        end if 
        end if 

        !------------------------------------------------------------------
        CALL CPU_TIME(finish_interpr)
        
        CALL CPU_TIME(finish_interp)

        ! t_interp = t_interp + finish_interp - start_interp
        ! t_dl = t_dl + finish_dl - start_dl
        ! t_interps = t_interps + finish_interps- start_interps
        ! t_interpr = t_interpr + finish_interpr- start_interpr
        ! t_interpp = t_interpp + finish_interpp- start_interpp
    enddo
    CALL CPU_TIME(finish)
    ! print*,'done.'
    ! print*, 'x loop time : ', finish-start, 's'
    ! print*, 'resolution: ', t_propa/(finish-start)*100, '%'
    ! print*, 'compute dl  : ', t_dl/(finish-start)*100, '%'
    ! print*, 'interpolation : ', t_interp/(finish-start)*100, '%'
    ! print*, 'interpolation side : ', t_interps/(finish-start)*100, '%'
    ! print*, 'interpolation rec : ', t_interpr/(finish-start)*100, '%'
    ! print*, 'interpolation plane : ', t_interpp/(finish-start)*100, '%'

    deallocate(pml%sigma)
    deallocate(pml%sigma_de)

    deallocate(AlphamBeta)
    deallocate(AlphapBeta)

    ! deallocate(phi_p1)
    ! deallocate(phi_m1)
    deallocate(u)

    deallocate(m%a_vect)
    deallocate(m%b_vect)
    deallocate(m%c_vect)

    deallocate(m%d_vect)
    deallocate(m%e_vect)
    deallocate(m%f_vect)
end subroutine propa_1angle_1freq_1domain
end module simple_wape
