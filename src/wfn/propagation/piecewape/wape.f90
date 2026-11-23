module wape
    implicit none
    public :: wape_1angle_1freq
contains
    subroutine wape_1angle_1freq(var0,tau,f,pe_model,source_height,imp,pml,sol,output)
        use types
        use, intrinsic :: ieee_arithmetic
        use utils
        use interpo
        use fast_interp_1D
        use initialize_flow
        ! use interp_2D_spline
        use read_write
        use hdf5
        use splines
        use impedance_model
        use create_mesh
        implicit none
        ! input parameters
        ! type(domain), intent(inout)                             :: dom
        type(real), intent(in)                                  :: f
        type(real), intent(in)                                  :: tau
        type(impedance), intent(inout)                          :: imp
        type(variables_0), intent(inout)                        :: var0
        type(perfectly_matched_layer),intent(inout)             :: pml
        real, intent(in)                                        :: source_height
        type(output_parameters),intent(in)                      :: output
        integer, intent(in)                                     :: pe_model
        integer                                                 :: nb_receiver
        ! output
        type(solution), intent(inout)                           :: sol
        type(complex), dimension(:,:), allocatable              :: phi0
        type(complex), dimension(:), allocatable                :: phi0_temp, tx,tz
        type(real), dimension(:,:), allocatable                 :: p_over_p0
        type(real), dimension(:), allocatable                   :: x0, z0
        ! local variables
        type(complex),allocatable, dimension(:)                 :: u
        type(solve)                                             :: m
        real                                                    :: omega
        complex                                                 :: kappa
        complex,allocatable,dimension(:)                        :: AlphamBeta, AlphapBeta


        ! interpolated flow
        ! generel variables
        type(frequency)                                         :: freq
        REAL                                                    :: pi = 3.141592653589793, tau_rad
        real                                                    :: ztest,a, xc, zc, ac
        integer                                                 :: nx, nz
        complex                                                 :: eye = (0.0,1.0), phi_interp
        integer                                                 :: ii,jj, kk, index, idom
        real, dimension(:,:), allocatable   :: x_coarse_grid_rot, z_coarse_grid_rot, x_rec_grid_rot, z_rec_grid_rot
        real                                                    :: abs_source_height, t

        real, parameter :: level_threshold_db = -100.0   ! tune to taste
        real            :: phi0_max_abs, phi0_level_db
        integer :: ixend


        ! time computation 
        real                                                    :: start, finish, t_tot
        real                                                    :: start_propa, finish_propa, t_propa
        real                                                    :: start_dl, finish_dl, t_dl
        real                                                    :: start_interp, finish_interp, t_interp
        real                                                    :: start_starter, finish_starter, t_starter
        real                                                    :: start_solve, finish_solve, t_solve
        real                                                    :: start_flow, finish_flow, t_flow

        t_flow = 0
        t_solve = 0
        t_interp =0 
        t_starter = 0

        print*, 'start vector PE propagation...'
        CALL CPU_TIME(start)
        ! init frequency param
        !----------------------------------------------------------------------

        ! init impedance 
        freq%k0 = 2*pi*f/var0%c
        freq%omega  = 2*pi*f
        call compute_beta(imp%mode, f, imp, var0)

        ! loop over the domain 
        !----------------------------------------------------------------------
        CALL CPU_TIME(start_flow)
        ! interpolate flow inside mesh 
        call init_meanflow_freq(pe_model,mesh,output%nbout,mesh%nx,mesh%nz,var0%c,1)
        CALL CPU_TIME(finish_flow)
        t_flow = t_flow + finish_flow - start_flow

        CALL CPU_TIME(start_solve)
        allocate(sol%phi(1, mesh%nz))
        allocate(sol%phi_p1(mesh%nz))
        allocate(sol%phi_m1(mesh%nz))
        allocate(sol%p_over_p0(1, mesh%nz))
        allocate(u(mesh%nz))


        ! Initalize the starter 
        !----------------------------------------------------------------------
        ! salomon start for a point source 
        sol%phi(1,:) = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(mesh%z-source_height))**2)* &
            & (exp(-((freq%k0*(mesh%z-source_height))**2)/3))
    

        ! sol%phi(1,:) = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(mesh%z-source_height))**2)* &
        !       exp(-(1./3.)*(freq%k0*(z_a-source_height))**2)
        sol%phi_m1 = 0.
        sol%phi_p1 = sol%phi(1,:)
        sol%p_over_p0(1,:) = max(abs(sqrt(mesh%x_grid(1,:)**2+(mesh%z_grid(1,:)-mesh%source_height_abs)**2)*sol%phi(1,:) &
            /sqrt(max(mesh%x_grid(1,:),1e-12))*exp(eye*freq%k0*mesh%x_grid(1,:))), 1e-10)


        !totally reflective ground hard code
        ! PML definition
        !-------------------------------------------------------------------------
        allocate(pml%sigma(mesh%nz))
        allocate(pml%sigma_de(mesh%nz))
        pml%sigma = 0.
        pml%sigma_de = 0.

        do ii = mesh%nz - pml%size, mesh%nz
        pml%sigma(ii) = pml%param * (real(ii- mesh%nz + pml%size)/real(pml%size))**pml%n
        pml%sigma_de(ii) = pml%param*pml%n*real((ii- mesh%nz + pml%size))**(pml%n-1)/(mesh%dz*(pml%size**pml%n))
        enddo


        ! init vector for matrix resolution
        !-------------------------------------------------------------------------
        allocate(m%h10(mesh%nz))
        allocate(m%h12(mesh%nz))
        allocate(m%h20(mesh%nz))
        allocate(m%h22(mesh%nz))

        kappa = eye*freq%k0*mesh%dx/2

        m%h10 = flow2D%h10(1,1:mesh%nz)
        m%h12 = flow2D%h12(1,1:mesh%nz)
        m%h20 = flow2D%h20(1,1:mesh%nz)
        m%h22 = flow2D%h22(1,1:mesh%nz)

        allocate(m%a_vect(mesh%nz))
        allocate(m%b_vect(mesh%nz))
        allocate(m%c_vect(mesh%nz))

        allocate(m%d_vect(mesh%nz))
        allocate(m%e_vect(mesh%nz))
        allocate(m%f_vect(mesh%nz))
        allocate(AlphamBeta(mesh%nz))
        allocate(AlphapBeta(mesh%nz))

        AlphapBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz**2)) - &
            & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz))
        AlphamBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz**2)) + &
            & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz))

        m%a_vect = (m%h10 - kappa*m%h20) - &
            2*(m%h12- kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz**2)
        m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
        m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


        m%d_vect = (m%h10 + kappa*m%h20) - &
            2*(m%h12 + kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz**2)
        m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
        m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

        m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
        m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

        m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*mesh%dz*m%b_vect(1);
        m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*mesh%dz*m%e_vect(1);


        ! propagation
        !----------------------------------------------------------------------------------------
        do kk = 2, mesh%nx
        index = (kk-1)/output%nbout+1
        ! update matrix coefficient with flmow grid 
        if(mod(kk-1,output%nbout) == 0) then
            m%h10 = flow2D%h10(index,1:mesh%nz)
            m%h12 = flow2D%h12(index,1:mesh%nz)
            m%h20 = flow2D%h20(index,1:mesh%nz)
            m%h22 = flow2D%h22(index,1:mesh%nz)

            m%a_vect = (m%h10 - kappa*m%h20) - &
                & 2*(m%h12-kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz**2)
            m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
            m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


            m%d_vect = (m%h10 + kappa*m%h20) - &
                & 2*(m%h12+kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*mesh%dz**2)
            m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
            m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

            m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
            m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

            m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*mesh%dz*m%b_vect(1);
            m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*mesh%dz*m%e_vect(1);
        endif

        ! solve system 
        CALL CPU_TIME(start_propa)
        u = sol%phi_p1(:)
        sol%phi_p1(1) = m%d_vect(1)*u(1) + m%f_vect(1)*u(2)
        sol%phi_p1(2:mesh%nz-1) = m%e_vect(2:mesh%nz-1)*u(1:mesh%nz-2) + m%d_vect(2:mesh%nz-1)*u(2:mesh%nz-1) &
            & + m%f_vect(2:mesh%nz-1)*u(3:mesh%nz)
        sol%phi_p1(mesh%nz) = m%e_vect(mesh%nz)*u(mesh%nz-1) + m%d_vect(mesh%nz)*u(mesh%nz)

        call res_trid(sol%phi_p1(1:mesh%nz), m%b_vect(2:mesh%nz), m%a_vect, m%c_vect(1:mesh%nz-1),mesh%nz)

        sol%phi(1,mesh%nz) = 0;
        CALL CPU_TIME(finish_propa)
        t_propa = t_propa + finish_propa - start_propa
        ! Compute delta L
        CALL CPU_TIME(start_dl)


        if ((pe_model == 2) .or. (pe_model == 4)) then 
            do ii=1,mesh%nz
            sol%p_over_p0(1,ii) = max(abs(sqrt(mesh%x_grid(kk,ii)**2+(mesh%z_grid(kk,ii)-source_height)**2)*sol%phi(1,ii) &
                /sqrt(mesh%x_grid(kk,ii))), 1e-10)
            enddo
        else

            if (kk>2) then 
                do ii=1,mesh%nz
                sol%p_over_p0(1,ii) =  sqrt((mesh%x_grid(kk-1,ii)**2 +(mesh%z_grid(kk-1,ii)-source_height)**2))/sqrt(mesh%x_grid(kk-1,ii)) * &
                    abs((1 - flow2D%m(index,ii))*sol%phi(1,ii) + &
                    (eye*flow2D%m(index,ii)/(2*freq%k0*mesh%dx))*(sol%phi_p1(ii)-sol%phi_m1(ii)))
                enddo
            end if 
        end if 


        ! save receiver 
        if (output%side) then
            call linear_interp(sol%p_over_p0(1,:), mesh%z, mesh%z_coarse,sol%deltaL_side_intermediate(kk,:))
        end if 

        call linear_interp(sol%p_over_p0(1,:),mesh%z,output%heights,sol%deltaL_rec_intermediate(kk,:))
        

        sol%phi_m1 = sol%phi(1,:)
        sol%phi(1,:) = sol%phi_p1
        enddo

        deallocate(pml%sigma)
        deallocate(pml%sigma_de)

        deallocate(AlphamBeta)
        deallocate(AlphapBeta)

        deallocate(u)
        deallocate(sol%p_over_p0)
        deallocate(sol%phi)
        deallocate(sol%phi_p1)
        deallocate(sol%phi_m1)
        deallocate(mesh%x_grid)
        deallocate(mesh%z_grid)

        deallocate(m%a_vect)
        deallocate(m%b_vect)
        deallocate(m%c_vect)
        deallocate(m%d_vect)
        deallocate(m%e_vect)
        deallocate(m%f_vect)

        deallocate(flow2D%h10)
        deallocate(flow2D%h12)
        deallocate(flow2D%h20)
        deallocate(flow2D%h22)
        deallocate(flow2D%m)


        CALL CPU_TIME(finish_solve)
        t_solve = t_solve + finish_solve - start_solve

        ! interpolation on the receivers 
        ! ---------------------------------------------------------------------
        CALL CPU_TIME(start_interp)
        if (output%side) then
            do ii=1,mesh%nz_coarse
            ! call interp2(mesh%x,sol%deltaL_side_intermediate(:,ii),mesh%x_coarse,sol%deltaL_side(:,ii))
            call linear_interp(sol%deltaL_side_intermediate(:,ii),mesh%x,mesh%x_coarse,sol%deltaL_side(:,ii))
            enddo
            sol%deltaL_side = 20*log10(max(sol%deltaL_side, 10e-10))
        end if 

        do ii=1,output%nb_receiver
        ! call interp2(mesh%x,sol%deltaL_rec_intermediate(:,ii),mesh%x_coarse,sol%deltaL_rec_intermediate(:,ii))
        call linear_interp(sol%deltaL_rec_intermediate(:,ii),mesh%x,mesh%x_coarse,sol%deltaL_rec(:,ii))
        enddo
        sol%deltaL_rec = 20*log10(max(sol%deltaL_rec, 10e-10))

        CALL CPU_TIME(finish_interp)
        t_interp = t_interp + finish_interp- start_interp


        ! deallocate(sol%phi)
        ! deallocate(sol%p_over_p0)
        CALL CPU_TIME(finish)
        print*,'done.'
        if ((finish - start) > 0.0) then 
            print*, 'x loop time : ', finish-start, 's'
            print*, 'init flow: ', t_flow/(finish-start)*100, '%'
            print*, 'solve : ', t_solve/(finish-start)*100, '%'
            print*, 'interpolation side : ', t_interp/(finish-start)*100, '%'
        else
            print*, 'init flow: N/A (total time too small to measure)'
            print*, 'solve : N/A'
            print*, 'fill starter : N/A'
            print*, 'interpolation side : N/A'
        endif 

    end subroutine
end module
