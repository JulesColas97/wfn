module gtpe

implicit none

public :: propa_1angle_1freq_gtpe
contains


subroutine propa_1angle_1freq_gtpe(var0,tau,f,source_height,imp,pml,flow2D,mesh,sol,output)
    use types
    use utils
    use interpo
    use fast_interp_1D
    use read_write
    use hdf5
    use impedance_model
    implicit none
    ! input parameters
    real, intent(in)                                        :: tau
    type(real), intent(in)                                  :: f
    type(impedance), intent(inout)                          :: imp
    type(variables_0), intent(inout)                        :: var0
    type(perfectly_matched_layer),intent(inout)             :: pml
    real, intent(in)                                        :: source_height
    type(mesh_2D), intent(inout)                            :: mesh
    type(meanflow2D), intent(in)                            :: flow2D
    type(output_parameters),intent(in)                      :: output
    integer                                                 :: nbout, nb_receiver
    ! output
    type(solution), intent(inout)                           :: sol

    ! generel variables
    type(frequency)                                         :: freq
    real                                                    :: pi = 3.1415, ztest
    complex                                                 :: eye = (0,1), phi_interp
    integer                                                 :: ii, kk, index
    real                                                    :: tau_rad
    integer, parameter                                      :: single = kind(1e0), double = kind(1d0)
    ! real, allocatable,dimension(:)                          :: H,Hp,Hpp
    ! local variable for the loop
    real                                                    :: temp_cl
    real                                                    :: H_a,Hp_a,Hpp_a
    real                                                    :: H_b, Hp_b, Hpp_b
    real, allocatable, dimension(:)                         :: z_a, z_b, zb_coarse
    complex,allocatable,dimension(:)                        :: Apml_a,Apml_b,Bpml_a,Bpml_b
    complex,allocatable,dimension(:)                        :: alpha_a, beta_a,chi_a, gamma_a
    complex,allocatable,dimension(:)                        :: alphap_a, betap_a, gammap_a
    complex,allocatable,dimension(:)                        :: alpha_b, beta_b,chi_b, gamma_b
    complex,allocatable,dimension(:)                        :: alphap_b, betap_b, gammap_b
    complex,allocatable,dimension(:)                        :: A_alpha, B_alpha
    complex,allocatable,dimension(:)                        :: A_chi, B_chi
    complex,allocatable,dimension(:)                        :: A_gamma, B_gamma
    complex,allocatable,dimension(:)                        :: B_l,B_u,B_d
    complex,allocatable,dimension(:)                        :: A_l,A_u,A_d
    complex                                                 :: A_a, B_a, C_a
    complex                                                 :: A_b,B_b, C_b
    complex                                                 :: L_CL,M_CL,P_CL,N_CL,K_CL,O_CL
    complex                                                 :: Q_CL,T_CL,R_CL,S_CL,U_CL
    complex                                                 :: phi_abarre
    complex,allocatable,dimension(:)                        :: phi_a,phi_b
    real :: abs_source_height

    real, allocatable, dimension(:)                         :: intermediate1
    complex                                                 ::  one_2ik0, intermediate2, intermediate3
    real :: start0, start1, finish
    real :: Ttotal, Tquantity,Tcl, Tsolve, Twrite

    Ttotal=0
    Tquantity=0
    Tcl=0
    Tsolve=0
    Twrite=0

    tau_rad = tau*pi/180
    nbout = output%nbout
    ! init frequency param
    freq%freq = f
    freq%omega = 2*pi*freq%freq
    freq%lambda = var0%c/freq%freq
    freq%k0 = 2*pi/freq%lambda

    ! imp%Zs = sqrt(4*imp%sigmae/(-eye*freq%omega*var0%gamma*var0%rho)) + var0%c*imp%alphae/(-eye*freq%omega*4*var0%gamma)
    ! imp%beta = 1/imp%Zs
    ! if (imp%rigid) then
    !   !totally reflective
    !   imp%beta = 0
    ! endif
        
    freq%k0 = 2*pi*f/var0%c
    freq%omega  = 2*pi*f
    call compute_beta(imp%mode, f, imp, var0)
    ! ground height definition
    ! allocate(H(mesh%nx))
    ! allocate(Hp(mesh%nx))
    ! allocate(Hpp(mesh%nx))
    ! call  prof(mesh%x,simu%hh,simu%hl,tau_rad,simu%hpos,H,Hp,Hpp)

    ! allocate vector for matrix resolution
    !-------------------------------------------------------------------------
    ! allocate mesh
    allocate(z_a(mesh%nz), z_b(mesh%nz))
    ! allocate pml
    allocate(pml%sigma_a(mesh%nz))
    allocate(pml%sigmap_a(mesh%nz))
    allocate(pml%sigma_b(mesh%nz))
    allocate(pml%sigmap_b(mesh%nz))
    allocate(Apml_a(mesh%nz),Apml_b(mesh%nz),Bpml_a(mesh%nz),Bpml_b(mesh%nz))

    ! allocate coef
    allocate(alpha_a(mesh%nz), beta_a(mesh%nz),chi_a(mesh%nz), gamma_a(mesh%nz),intermediate1(mesh%nz))
    allocate(alphap_a(mesh%nz), betap_a(mesh%nz), gammap_a(mesh%nz))
    allocate(alpha_b(mesh%nz), beta_b(mesh%nz),chi_b(mesh%nz), gamma_b(mesh%nz))
    allocate(alphap_b(mesh%nz), betap_b(mesh%nz), gammap_b(mesh%nz))
    allocate( A_alpha(mesh%nz), B_alpha(mesh%nz), A_chi(mesh%nz), B_chi(mesh%nz),  A_gamma(mesh%nz), B_gamma(mesh%nz))
    allocate(B_l(mesh%nz),B_u(mesh%nz),B_d(mesh%nz))
    allocate(A_l(mesh%nz),A_u(mesh%nz),A_d(mesh%nz))

    print *, 'Module BC-SW'

    ! intermediate 
    one_2ik0 = 1/(2*eye*freq%k0)
    ! init A
    !-----------------------------------------------------------------------
    !init Height
    H_a  = mesh%h(1)
    Hp_a = mesh%hp(1)
    Hpp_a = mesh%hpp(1)
    z_a = mesh%z_grid(1,:)
    ! z_a = (mesh%z*(mesh%Lz - H_a)/mesh%Lz) + H_a

    !init pml
    pml%iza = minloc(abs(z_a-mesh%Lz+pml%size*mesh%dz),1)
    if (z_a(pml%iza)<(mesh%Lz-pml%size*mesh%dz)) pml%iza=pml%iza+1
    pml%sigma_a = 0.
    pml%sigmap_a = 0.
    do ii = pml%iza, mesh%nz
        pml%sigma_a(ii) = pml%param*((z_a(ii)-mesh%Lz+pml%size*mesh%dz)/(pml%size*mesh%dz))**pml%n
        pml%sigmap_a(ii) = pml%param*pml%n*(z_a(ii)- mesh%Lz + pml%size*mesh%dz)**(pml%n-1)/(((pml%size*mesh%dz)**pml%n))
    enddo
    Apml_a = 1./(1.+eye*pml%sigma_a/freq%omega)**2
    Bpml_a =-eye*pml%sigmap_a/freq%omega/(1+eye*pml%sigma_a/freq%omega)**3

    ! init alpha, beta, chi,gamma
    alpha_a = Apml_a*mesh%Lz**2./(mesh%Lz-H_a)**2 + Hp_a**2*(mesh%Lz-mesh%z)**2/(mesh%Lz-H_a)**2
    beta_a = Bpml_a*mesh%Lz/(mesh%Lz-H_a) - Hpp_a*(mesh%Lz-mesh%z)/(mesh%Lz-H_a) &
            - 2*Hp_a**2*(mesh%Lz-mesh%z)/(mesh%Lz-H_a)**2 - 2*eye*freq%k0*Hp_a*(mesh%Lz-mesh%z)/(mesh%Lz-H_a)
    chi_a =  Bpml_a*mesh%Lz/(mesh%Lz-H_a) + Hpp_a*(mesh%Lz-mesh%z)/(mesh%Lz-H_a) - &
                      2*eye*freq%k0*Hp_a*(mesh%Lz-mesh%z)/(mesh%Lz-H_a)
    gamma_a = freq%k0**2*flow2D%epsilon_eff(1,:)

    ! init alphap, betap, gammap
    alphap_a = -alpha_a/(2*eye*freq%k0)
    betap_a = -beta_a/(2*eye*freq%k0) - 2*Hp_a*(mesh%Lz-mesh%z)/(mesh%Lz-H_a)
    gammap_a = 2*eye*freq%k0 - gamma_a/(2*eye*freq%k0)

    ! Boundary condition on the ground v3 (cf Sak et Wes)
    A_a = -Hp_a/(2*mesh%dx)
    B_a = (1+Hp_a**2)/(2*mesh%dz)*mesh%Lz/(mesh%Lz-H_a)
    C_a = (eye*freq%k0*sqrt(1+Hp_a**2)*imp%beta - eye*freq%k0*Hp_a)

    ! Init solution
    !-----------------------------------------------------------------------
    ! allocate(sol%phi(mesh%nz,1,1))
    allocate(phi_a(mesh%nz))
    allocate(phi_b(mesh%nz))
    ! allocate(sol%phi_temp(mesh%nz,1,1))
    ! allocate(sol%deltaL(mesh%nz))
    allocate(sol%p_over_p0(mesh%nz,1))

    ! sol%phi = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(z_a-src%pos_z))**2)* &
    !           exp(-(1./3.)*(freq%k0*(z_a-src%pos_z))**2)
    abs_source_height = source_height + z_a(1)
    phi_b = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(z_a-abs_source_height))**2)* &
              exp(-(1./3.)*(freq%k0*(z_a-abs_source_height))**2)
    phi_abarre = 0.

    ! print*, flow2D%epsilon_eff
    ! x advancement
    !-----------------------------------------------------------------------
    print *, "begin x loop ..."
    CALL CPU_TIME(start0)
    do kk = 2,mesh%nx-1
      index = (kk)/nbout+1

      ! init B
      !-----------------------------------------------------------------------
      !init Height
      H_b  = mesh%h(kk)
      Hp_b = mesh%hp(kk)
      Hpp_b = mesh%hpp(kk)
      ! z_b = (mesh%z*(mesh%Lz - H_b)/mesh%Lz) + H_b
      z_b = mesh%z_grid(kk,:)
      !init pml
      pml%izb = minloc(abs(z_b-mesh%Lz+pml%size*mesh%dz),1)
      if (z_b(pml%izb)<(mesh%Lz-pml%size*mesh%dz)) then
        pml%izb=pml%izb+1
      endif

      pml%sigma_b = 0.
      pml%sigmap_b = 0.
      do ii = pml%izb, mesh%nz
          pml%sigma_b(ii) = pml%param*(abs((z_b(ii)-mesh%Lz+pml%size*mesh%dz)/(pml%size*mesh%dz)))**pml%n
          pml%sigmap_b(ii) = pml%param*pml%n*(abs(z_b(ii)- mesh%Lz + pml%size*mesh%dz))**(pml%n-1)/(((pml%size*mesh%dz)**pml%n))
      enddo
      Apml_b = 1./(1+eye*pml%sigma_b/freq%omega)**2
      Bpml_b =-(eye*pml%sigmap_b/freq%omega)/((1+eye*pml%sigma_b/freq%omega)**3)

      CALL CPU_TIME(start1)
      ! init alpha, beta, chi,gamma
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>> without intermediate
      alpha_b = Apml_b*mesh%Lz**2./(mesh%Lz-H_b)**2 + (Hp_b**2)*(mesh%Lz-mesh%z)**2/((mesh%Lz-H_b)**2)

      beta_b = Bpml_b*mesh%Lz/(mesh%Lz-H_b) - Hpp_b*(mesh%Lz-mesh%z)/(mesh%Lz-H_b) &
              - 2*Hp_b**2*(mesh%Lz-mesh%z)/(mesh%Lz-H_b)**2 - 2*eye*freq%k0*Hp_b*(mesh%Lz-mesh%z)/(mesh%Lz-H_b)
      chi_b =  Bpml_b*mesh%Lz/(mesh%Lz-H_b) + Hpp_b*(mesh%Lz-mesh%z)/(mesh%Lz-H_b) - &
              2*eye*freq%k0*Hp_b*(mesh%Lz-mesh%z)/(mesh%Lz-H_b)

      gamma_b = freq%k0**2*flow2D%epsilon_eff(floor((kk)/nbout + 1.),:)

      ! init alphap, betap, gammap
      alphap_b = -alpha_b*one_2ik0
      betap_b = -beta_b*one_2ik0 - 2*Hp_b*(mesh%Lz-mesh%z)/(mesh%Lz-H_b)
      gammap_b = 2*eye*freq%k0 - gamma_b*one_2ik0

      ! calculate intermediate value
      A_alpha = mesh%dx*(alpha_a/3 + alpha_b/6)
      B_alpha = mesh%dx*(alpha_a/6 + alpha_b/3)

      A_chi = mesh%dx*(chi_a/3 + chi_b/6)
      B_chi = mesh%dx*(chi_a/6 + chi_b/3)

      A_gamma = mesh%dx*(gamma_a/3 + gamma_b/6)
      B_gamma = mesh%dx*(gamma_a/6 + gamma_b/3)

      ! calculate matrix coeficients
      B_l = (alphap_b + B_alpha)/(mesh%dz**2) - (betap_b + B_chi)/(2*mesh%dz)
      B_u = (alphap_b + B_alpha)/(mesh%dz**2) + (betap_b + B_chi)/(2*mesh%dz)
      B_d = -2*(alphap_b + B_alpha)/(mesh%dz**2) + (gammap_b + B_gamma);

      A_l = (alphap_a - A_alpha)/(mesh%dz**2) - (betap_a - A_chi)/(2*mesh%dz)
      A_u = (alphap_a - A_alpha)/(mesh%dz**2) + (betap_a - A_chi)/(2*mesh%dz)
      A_d = -2*(alphap_a - A_alpha)/(mesh%dz**2) + (gammap_a - A_gamma)


      CALL CPU_TIME(finish)
      Tquantity = Tquantity + (finish-start1)

      CALL CPU_TIME(start1)

      ! boundary condition v3 saak et Wes
      ! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      A_b = -Hp_b/(2*mesh%dx)
      B_b = (1+Hp_b**2)/(2*mesh%dz)*mesh%Lz/(mesh%Lz-H_b)
      C_b = (eye*freq%k0*sqrt(1+Hp_b**2)*imp%beta - eye*freq%k0*Hp_b)

      intermediate2 = 1/(3*B_a-C_a)
      intermediate3 = 1/(3*B_b-C_b-3*A_b)

      L_CL = 1 + A_a*intermediate2*(4*A_b)*intermediate3
      M_CL = 4*B_a*intermediate2
      P_CL = -B_a*intermediate2
      N_CL = A_a*intermediate2*4*B_b*intermediate3
      K_CL = -A_a*intermediate2*B_b*intermediate3
      O_CL = A_a*intermediate2*A_b*intermediate3 - A_a*intermediate2
      Q_CL = -4*A_b*intermediate3*4*B_a*intermediate2
      T_CL = 4*A_b*intermediate3*B_a*intermediate2
      R_CL = 4*B_b*intermediate3
      U_CL = -B_b*intermediate3
      S_CL = A_b*intermediate3 + (4*A_b)*intermediate3*A_a*intermediate2


      B_d(2) = L_CL*B_d(2) + R_CL*B_l(2) - N_CL*A_l(2)
      B_u(2) = L_CL*B_u(2) + U_CL*B_l(2) - K_CL*A_l(2)

      A_d(2) = L_CL*A_d(2) + M_CL*A_l(2) - Q_CL*B_l(2)
      A_u(2) = L_CL*A_u(2) + P_CL*A_l(2) - T_CL*B_l(2)

      CALL CPU_TIME(finish)
      Tcl = Tcl + (finish-start1)
      CALL CPU_TIME(start1)

      ! resolution v3 Saak et We
      ! phi_a = sol%phi
      ! ! A*phi(a) 
      ! sol%phi_temp(2) = A_d(2)*sol%phi(2) + A_u(2)*sol%phi(3)
      ! sol%phi_temp(3:mesh%nz-1) = A_l(3:mesh%nz-1)*sol%phi(2:mesh%nz-2) + A_d(3:mesh%nz-1)*sol%phi(3:mesh%nz-1) &
      !     & + A_u(3:mesh%nz-1)*sol%phi(4:mesh%nz)
      ! sol%phi_temp(mesh%nz) = A_l(mesh%nz)*sol%phi(mesh%nz-1) + A_d(mesh%nz)*sol%phi(mesh%nz)
      !
      !
      ! sol%phi_temp(2) = sol%phi_temp(2) + (O_CL*A_l(2) - S_CL*B_l(2))*phi_abarre
      ! call res_trid(sol%phi_temp(2:mesh%nz),B_l(3:mesh%nz), B_d(2:mesh%nz), B_u(2:mesh%nz-1),mesh%nz-1)
      !
      ! sol%phi = sol%phi_temp
      ! phi_a(1) = M_CL/L_CL*phi_a(2) + P_CL/L_CL*phi_a(3) + &
      !   N_CL/L_CL*sol%phi(2) + K_CL/L_CL*sol%phi(3) + O_CL/L_CL*phi_abarre
      !
      ! sol%phi(1) = Q_CL/L_CL*phi_a(2) + T_CL/L_CL*phi_a(3) + &
      !     R_CL/L_CL*sol%phi(2) + U_CL/L_CL*sol%phi(3) + S_CL/L_CL*phi_abarre
      !
      ! sol%phi(mesh%nz) = 0

      ! A*phi(a) 
      phi_a = phi_b
      phi_b(2) = A_d(2)*phi_a(2) + A_u(2)*phi_a(3)
      phi_b(3:mesh%nz-1) = A_l(3:mesh%nz-1)*phi_a(2:mesh%nz-2) + A_d(3:mesh%nz-1)*phi_a(3:mesh%nz-1) &
          & + A_u(3:mesh%nz-1)*phi_a(4:mesh%nz)
      phi_b(mesh%nz) = A_l(mesh%nz)*phi_a(mesh%nz-1) + A_d(mesh%nz)*phi_a(mesh%nz)
      phi_b(2) = phi_b(2) + (O_CL*A_l(2) - S_CL*B_l(2))*phi_abarre

      call res_trid(phi_b(2:mesh%nz),B_l(3:mesh%nz), B_d(2:mesh%nz), B_u(2:mesh%nz-1),mesh%nz-1)

      phi_a(1) = M_CL/L_CL*phi_a(2) + P_CL/L_CL*phi_a(3) + &
        N_CL/L_CL*phi_b(2) + K_CL/L_CL*phi_b(3) + O_CL/L_CL*phi_abarre

      phi_b(1) = Q_CL/L_CL*phi_a(2) + T_CL/L_CL*phi_a(3) + &
          R_CL/L_CL*phi_b(2) + U_CL/L_CL*phi_b(3) + S_CL/L_CL*phi_abarre

      phi_b(mesh%nz) = 0
      CALL CPU_TIME(finish)
      Tsolve = Tsolve + (finish-start1)
      CALL CPU_TIME(start1)


     ! interpolation on receivers 
      ! if (kk>=2) then
      !   ! Interpolation for carto
      !   do ii=1,mesh%nz
      !   sol%deltaL(ii) =  20*log10(MAX(abs(sol%phi(ii)/sqrt(mesh%xi(kk))*exp(eye*freq%k0*mesh%xi(kk)))* &
      !                                       & sqrt(mesh%xi(kk-1)**2 +(z_b(ii)-src%pos_z)**2),1e-10))
      !   enddo
      !   zb_coarse = (mesh%eta_coarse*(mesh%Lz-H_b)/mesh%Lz)+H_b
      !
      !   call interplog(z_b,sol%deltaL,zb_coarse,sol%deltaL_vect(kk,:))
      !   call interplog(z_b,sol%deltaL,output%heights+H_b,sol%receiver(kk,:))
      !
      ! endif


      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ! Compute delta L
        sol%p_over_p0(:,1) =abs(sqrt(mesh%x(kk)**2+(z_b-abs_source_height)**2)*phi_b(:) &
            & /sqrt(mesh%x(kk))*exp(eye*freq%k0*mesh%x(kk)))

        zb_coarse = (mesh%z_coarse*(mesh%Lz_coarse-h_b)/mesh%Lz_coarse)+h_b

        ! print *, sol%p_over_p0(10)
        ! Interpolation for different receiver and planes 
        !------------------------------------------------------------------
        ! interpolation on coarse grid for side view
        if (output%side) then
            ! call interpcomplextoreal(z_b,sol%p_over_p0(:,1,1),zb_coarse,sol%deltaL_vect(kk,:))
            ! call interp2(z_b,sol%p_over_p0(:,1,1),mesh%z_coarse,sol%deltaL_side_intermediate(kk,:))
            call linear_interp(sol%p_over_p0(:,1), z_b, mesh%z_coarse,sol%deltaL_side_intermediate(kk,:))
                    ! sol%deltaL_side_intermediate(kk,:) = 20*log10(MAX(sol%deltaL_side_intermediate(kk,:),1e-10))
        end if 

        ! call interp2(z_b,sol%p_over_p0(:,1,1),output%heights + H_b,sol%deltaL_rec_intermediate(kk,:))
        call linear_interp(sol%p_over_p0(:,1),z_b,output%heights + H_b,sol%deltaL_rec_intermediate(kk,:))
        ! if (kk<mesh%nx)  then
        !     ! interpolation at given point to construct final plane
        !     !--------------------------------------------------------------
        !     ! CALL CPU_TIME(start_interpp)
        !     do ii=1,size(output%xplane)
        !     if (output%xplane(ii)>=0) then
        !         if ((mesh%x(kk)*cos(tau_rad)<=output%xplane(ii)) &
        !             .and. (mesh%x(kk+1)*cos(tau_rad)>=output%xplane(ii))) then 
        !             call interpcomplextoreal(z_b,sol%p_over_p0,zb_coarse,sol%xplane(ii,:))
        !             sol%xplane(ii,:) = 20*log10(MAX(sol%xplane(ii,:),1e-10))
        !
        !             sol%ycoord(ii) = mesh%x(kk)*sin(tau_rad)
        !             sol%xcount(ii) = 1 
        !         endif 
        !     else 
        !         if ((mesh%x(kk)*cos(tau_rad)>=output%xplane(ii)) &
        !             .and. (mesh%x(kk+1)*cos(tau_rad)<=output%xplane(ii))) then 
        !             call interpcomplextoreal(z_b,sol%p_over_p0,zb_coarse,sol%xplane(ii,:))
        !             sol%xplane(ii,:) = 20*log10(MAX(sol%xplane(ii,:),1e-10))
        !             sol%ycoord(ii) = mesh%x(kk)*sin(tau_rad)
        !             sol%xcount(ii) = 1 
        !         endif
        !     endif  
        !     enddo 
        !
        !     do ii=1,size(output%yplane)
        !     if (output%yplane(ii)>0) then 
        !         if ((mesh%x(kk)*sin(tau_rad)<=output%yplane(ii)) &
        !             .and. (mesh%x(kk+1)*sin(tau_rad)>=output%yplane(ii))) then 
        !             call interpcomplextoreal(z_b,sol%p_over_p0,zb_coarse,sol%yplane(ii,:))
        !             sol%yplane(ii,:) = 20*log10(MAX(sol%yplane(ii,:),1e-10))
        !             sol%xcoord(ii) = mesh%x(kk)*cos(tau_rad)
        !             sol%ycount(ii) =  1
        !         endif 
        !     else 
        !         if ((mesh%x(kk)*sin(tau_rad)>=output%yplane(ii)) &
        !             .and. (mesh%x(kk+1)*sin(tau_rad)<=output%yplane(ii))) then 
        !             call interpcomplextoreal(z_b,sol%p_over_p0,zb_coarse,sol%yplane(ii,:))
        !             sol%yplane(ii,:) = 20*log10(MAX(sol%yplane(ii,:),1e-10))
        !             sol%xcoord(ii) = mesh%x(kk)*cos(tau_rad)
        !             sol%ycount(ii) =  1
        !         endif 
        !     endif 
        !     enddo 
        !     ! CALL CPU_TIME(finish_interpp)
        ! endif 

        !------------------------------------------------------------------
        ! CALL CPU_TIME(start_interpr)
        ! call interpcomplextoreal(z_b,sol%p_over_p0,output%heights+H_b,sol%receiver(kk,:)) 
        ! sol%receiver(kk,:) = 20*log10(MAX(sol%receiver(kk,:),1e-10))

        ! CALL CPU_TIME(finish_interpr)
        
        ! CALL CPU_TIME(finish_interp)

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><


      CALL CPU_TIME(finish)
      Twrite = Twrite + (finish-start1)

      H_a  =H_b
      Hp_a = Hp_b
      Hpp_a = Hpp_b
      z_a = z_b
      !init pml
      pml%iza = pml%izb
      pml%sigma_a = pml%sigma_b
      pml%sigmap_a = pml%sigmap_b
      Apml_a = Apml_b
      Bpml_a =Bpml_b
      ! init alpha, beta, chi,gamma
      alpha_a = alpha_b
      beta_a = beta_b
      chi_a = chi_b
      gamma_a = gamma_b
      ! init alphap, betap, gammap
      alphap_a = alphap_b
      betap_a = betap_b
      gammap_a =gammap_b
      ! Boundary condition on the ground
      A_a = A_b
      B_a = B_b
      C_a = C_b
      phi_abarre = phi_a(1)
    enddo



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

    CALL CPU_TIME(finish)
    Ttotal = finish-start0
    print *,'done.'
    print *, 'x loop time : ', Ttotal, 's'
    print *,'Quantity : ', Tquantity/Ttotal*100, '%'
    print *,'CL : ', Tcl/Ttotal*100,'%'
    print *,'Solve : ', Tsolve/Ttotal *100, '%'
    print *,'Interpo :', Twrite/Ttotal*100, '%'

    deallocate(pml%sigma_a)
    deallocate(pml%sigmap_a)
    deallocate(pml%sigma_b)
    deallocate(pml%sigmap_b)
    deallocate(Apml_a,Apml_b,Bpml_a,Bpml_b)

    ! deallocate coef
    deallocate(alpha_a, beta_a,chi_a, gamma_a)
    deallocate(alphap_a, betap_a, gammap_a)
    deallocate(alpha_b, beta_b,chi_b, gamma_b,intermediate1)
    deallocate(alphap_b, betap_b, gammap_b)
    deallocate(A_alpha, B_alpha, A_chi, B_chi,  A_gamma, B_gamma)
    deallocate(B_l,B_u,B_d)
    deallocate(A_l,A_u,A_d)
    ! deallocate(H)
    ! deallocate(Hp)
    ! deallocate(Hpp)
    deallocate(z_a)
    deallocate(z_b)

    !dealocate solution
    ! deallocate(sol%phi)
    deallocate(sol%p_over_p0)
    ! deallocate(sol%phi_temp)
    ! deallocate(sol%deltaL)
    deallocate(phi_a)
    deallocate(phi_b)

end subroutine propa_1angle_1freq_gtpe
end module gtpe


