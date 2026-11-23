module splitstep_pe_flat_mod
    use kind_mod, only: wp
    use constants, only: pi, UNSET, eye
    use input_params, only: inputs_t
    use mesh_mod, only: mesh_t
    use flow_mod, only: flow2D_t
    implicit none
    public :: run_splitstep_pe_flat

contains

    subroutine run_splitstep_pe_flat(freq, cfg, mesh, flow2D, p_out_side, p_out_rec)
        use iim_module, only: iim_stencil
        use pade_mod, only: qp, pade_propagator, pade_starter
        use trisolve_mod, only: apply_rational, tri_step, apply_rational_new
        use interp_1D_mod, only: linear_interp
        implicit none 
        real(wp), intent(in) :: freq
        type(inputs_t), intent(inout) :: cfg
        type(mesh_t), intent(in) :: mesh
        type(flow2D_t), intent(in) :: flow2D
        complex(wp), allocatable, intent(out) :: p_out_side(:,:), p_out_rec(:,:)

        real(wp) :: k0, omega

        real(wp) :: sigma(mesh%nz), sigma_de(mesh%nz)
        complex(wp) :: pml_term(mesh%nz), keff2(mesh%nz)
        complex(wp) :: Xsup(mesh%nz), Xsub(mesh%nz), Xdia(mesh%nz)
        integer :: nb_bot

        complex(wp) :: J(mesh%nx)
        complex(wp) :: c_lo(mesh%nx, 3), c_hi(mesh%nx, 3)

        complex(wp) :: h0, A0(cfg%n_pade), A1(cfg%n_pade), B0(cfg%n_pade), B1(cfg%n_pade)
        
        ! starter 
        integer     :: iz
        real(wp)    :: w, Xs
        complex(wp) :: phi(mesh%nz), rhs(mesh%nz), cp(mesh%nz)
        complex(wp) :: X0sub(mesh%nz), X0dia(mesh%nz), X0sup(mesh%nz)
        complex(wp) :: Tsub(mesh%nz), Tdia(mesh%nz), Tsup(mesh%nz)
        complex(wp) :: zeros(mesh%nz), ones(mesh%nz)

        integer :: ii, kk, j_coarse
        integer :: t0, t1, t2, t11, t12
        real(wp) :: rate


        ! output 
        ! complex(wp), allocatable :: p_out_side0(:,:), p_out_rec0(:,:)

        ! allocate(p_out_rec0(mesh%nx, size(cfg%receiver_heights)))
        ! allocate(p_out_side0(mesh%nx, mesh%nz_coarse))
        call system_clock(t0, rate)
        allocate(p_out_rec(mesh%nx_coarse, size(cfg%receiver_heights)))
        allocate(p_out_side(mesh%nx_coarse, mesh%nz_coarse))
        ! allocate(p_out_rec(mesh%nx, size(cfg%receiver_heights)))
        ! allocate(p_out_side(mesh%nx, mesh%nz))

        omega = 2 * pi * freq
        k0 = omega / cfg%c_ref 
        ! build matrix X 
        ! -----------------------------------------------------------------
        ! PML
        sigma = 0.0_wp
        sigma_de = 0.0_wp
        ! PML top 
        do ii = mesh%nz - cfg%pml_size, mesh%nz
        sigma(ii) = cfg%pml_param * (real(ii- mesh%nz + cfg%pml_size)/real(cfg%pml_size))**cfg%pml_n
        sigma_de(ii) = cfg%pml_param*cfg%pml_n*real((ii- mesh%nz + cfg%pml_size))**(cfg%pml_n-1) &
            /(mesh%dz*(cfg%pml_size**cfg%pml_n))
        enddo

        ! Bottom condition 
        if (.not. any(cfg%seabed_type == [0, 1, 2])) then 
            print*, "WARNING seabed type is not properly set, resuming to default : absorbing"
            cfg%seabed_type = 0 
        endif 
        if (cfg%seabed_type == 1) then 
            print*, "SEABED set to rigid"
            cfg%beta_sb = 0.0
        endif 

        if (cfg%beta_sb == UNSET) then 
            cfg%beta_sb = cfg%rho_water * cfg%c_water / &
                (cfg%rho_sb  * cfg%c_sb)
        endif 

        if (cfg%seabed_type == 0) then 
            nb_bot = int(ceiling(cfg%d_pml_bot * cfg%c_water/freq/mesh%dz))
        else
            nb_bot = 0
        endif 

        if (nb_bot>0) then 
            ! do ii=1, nb_bot
            ! sigma(ii) = cfg%pml_param * (real(nb_bot - ii)/real(nb_bot))**cfg%pml_n
            ! sigma_de(ii) = cfg%pml_param * cfg%pml_n * (real(nb_bot - ii)/real(nb_bot))**(cfg%pml_n-1) &
            !     / (mesh%dz*(cfg%pml_size**cfg%pml_n))
            ! enddo 
            do ii = 1, nb_bot
                sigma(ii)    =  cfg%pml_param * (real(nb_bot - ii + 1, wp)/nb_bot)**cfg%pml_n
                sigma_de(ii) = -cfg%pml_param * cfg%pml_n * real(nb_bot - ii + 1, wp)**(cfg%pml_n - 1) &
                               / (mesh%dz * real(nb_bot, wp)**cfg%pml_n)
            end do
        endif 

        ! d2/dz2 --> (1/(1+sigma/omega) d2/dz2 - i sigma'/omega*(1+isigma/omgea)**3 d/dz
        ! X = epsilon + 1/k0² dz/dz²
        ! basic form 
        keff2 = (k0 + eye*sigma / cfg%c_ref) ** 2
        pml_term = eye * sigma_de / ((omega  + eye * sigma) *2*mesh%dz)
        Xsub = (1/keff2) * (1/mesh%dz**2 + pml_term)
        Xsup = (1/keff2) * (1/mesh%dz**2 - pml_term) 
        Xdia = (-2/(keff2*mesh%dz**2))

        ! seabed impedance
        if ( any(cfg%seabed_type ==[1,2])) then 
            Xdia(1) = Xdia(1) + 2*eye *2 *pi *freq/cfg%c_water * cfg%beta_sb * Xsub(1)
            Xsup(1) = Xsup(1) + Xsub(1)
        endif 

        ! interface 
        ! compute updated coef for second order derivate in z, above and bellow
        ! the interface 
        J = k0**2 * (flow2D%epsilon_eff(mesh%js, :) - flow2D%epsilon_eff(mesh%js+1,:))
        call iim_stencil(mesh%z(mesh%js-1:mesh%js+1), [-1,-1, 1], J, cfg%rho_water, cfg%rho_air, 4, c_lo)
        call iim_stencil(mesh%z(mesh%js  :mesh%js+2), [-1, 1, 1], J, cfg%rho_water, cfg%rho_air, 5, c_hi)
        c_lo = c_lo / k0**2
        c_hi = c_hi / k0**2


        

        ! self starter 
        ! psi(r0) = sqrt(2 pi / k0) e^{i pi/4} (1+X)^(-1/4) e^{i k0 r0 (sqrt(1+X)-1)} delta_zS
        ! evaluated as  g(X) [T(X)^-2 delta],  g = T^2 (1+X)^(-1/4) exp(...),  T = 1 - X/Xs
        ! -----------------------------------------------------------------
        call X_row(1, X0sub, X0dia, X0sup) 
        Xs = maxval(flow2D%epsilon_eff) + 1.0_wp
        
        ! delta function at the source 
        ! spread between 2 points 
        iz = int(floor((cfg%source_height - mesh%z(1))/mesh%dz)) + 1
        w  = (cfg%source_height - mesh%z(iz))/mesh%dz
        phi = (0.0_wp, 0.0_wp)
        phi(iz)   = (1.0_wp - w)/mesh%dz
        phi(iz+1) = w/mesh%dz


        ! T(X) = 1 - X/Xs, applied twice as a solve: T phi_new = phi
        zeros = (0.0_wp, 0.0_wp)
        ones  = (1.0_wp, 0.0_wp)
        Tsub = -X0sub/Xs
        Tdia = 1.0_wp - X0dia/Xs
        Tsup = -X0sup/Xs
        do ii = 1, 2
            call tri_step(phi, Tsub, Tdia, Tsup, zeros, ones, zeros, rhs, cp)
        end do
        
        call pade_starter(real(k0, qp)*real(mesh%x(1), qp), cfg%n_pade, &
                  real(cfg%theta_pade, qp), real(Xs, qp), &
                  h0, A0, A1, B0, B1)
        ! g(X) from the starter Pade factors (start_h0, start_A0... from pade_starter)
        call apply_rational_new(phi, h0, A0, A1, B0, B1, &
                                X0sub, X0dia, X0sup, rhs, cp)
        phi = phi * sqrt(2*pi/k0) * exp(eye*pi/4)



        ! paddé coef 
        ! -----------------------------------------------------------------
        call pade_propagator(real(k0, qp)*real(mesh%dx, qp), cfg%n_pade, &
            real(cfg%theta_pade, qp), &
            h0, A0, A1, B0, B1)
        
        call system_clock(t1)
        t2 = 0
        ! marching schem 
        ! -----------------------------------------------------------------
        do kk=2,mesh%nx
            call X_row(kk, X0sub, X0dia, X0sup)
            call system_clock(t11)
            call apply_rational_new(phi, h0, A0, A1, B0, B1,  &
                    X0sub, X0dia, X0sup, rhs, cp)

            call system_clock(t12)
            t2 = t2 + t12 - t11
            phi(mesh%nz) = 0.0

            ! interpolation to common grid
            if (mod(kk, mesh%nbout) == 0) then
                j_coarse = kk / mesh%nbout
                ! interpolate along z only (phi varies slowly enough in z)
                call linear_interp(phi, mesh%z(1), mesh%dz, mesh%z_coarse, p_out_side(j_coarse,:))
                call linear_interp(phi, mesh%z(1), mesh%dz, cfg%receiver_heights, p_out_rec(j_coarse,:))
                ! far-field factors at the exact range
                p_out_side(j_coarse,:) = p_out_side(j_coarse,:) * exp(eye*k0*mesh%x(kk)) / sqrt(mesh%x(kk))
                p_out_rec(j_coarse,:)  = p_out_rec(j_coarse,:)  * exp(eye*k0*mesh%x(kk)) / sqrt(mesh%x(kk))
            end if
        enddo 
        print *, "pade:", real(t1-t0)/rate, " march:", real(t2)/rate

        
    contains

        ! Operator X at range index kk: regular rows + epsilon, with the two
        ! rows next to the interface replaced by the IIM stencils.
        subroutine X_row(jj, sub, dia, sup)
            integer,     intent(in)  :: jj
            complex(wp), intent(out) :: sub(:), dia(:), sup(:)
            integer :: js
            js = mesh%js

            sub = Xsub
            sup = Xsup
            dia = Xdia + flow2D%epsilon_eff(:, jj)

            ! row js (water side)
            sub(js) = c_lo(jj, 1)
            dia(js) = flow2D%epsilon_eff(js, jj) + c_lo(jj, 2)
            sup(js) = c_lo(jj, 3)

            ! row js+1 (air side)
            sub(js+1) = c_hi(jj, 1)
            dia(js+1) = flow2D%epsilon_eff(js+1, jj) + c_hi(jj, 2)
            sup(js+1) = c_hi(jj, 3)
        end subroutine X_row

    endsubroutine 
end module 
