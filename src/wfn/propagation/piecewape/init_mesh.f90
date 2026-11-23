module create_mesh
    use types

    implicit none
    public :: init_coarse_mesh
    public :: clean_coarse_mesh
    public :: init_mesh
    public :: init_receiver
    public :: compute_x_grid
    public :: compute_z_grid

    type(mesh_2D) :: mesh 

contains

    !! @brief define the parameters and the mesh of the 2D simulation coreseponding to the kk angle
    !!
    !! @details the function first determine the dimension of the 2D mesh from the total dimension of the simulation
    !! and from the angle, The it define the corresponding mesh
    !!
    !!
    !! @param [in]      simu: parameters of the global simulation
    !! @param [in]      tau: angle
    !! @param [in]      dout: coarse deltax
    !! @param[inout]    mesh: the 2D mesh coresponding to tau propagation
    !! @param[inout]    dom: parameters of the 2D mesh coresponding to the kk angles
    !! @return
    subroutine init_coarse_mesh(tau,simu,terrain2D, pml, output, c0)
        use types
        use utils 
        use interpo
        implicit none

        type(simulation), intent(in)                            :: simu
        type(terrain_elevation), intent(in)                     :: terrain2D
        type(perfectly_matched_layer), intent(in)               :: pml
        real, intent(in)                                        :: tau, c0
        type(output_parameters), intent(in) :: output


        integer                                                 :: ii
        real                                                    :: pi = 3.141592653589793, tau_rad
        real                                                    :: z_min, z_max, f_min, L_pml
        real, allocatable, dimension(:) ::  x_temp

        ! conversion in rad
        tau_rad = tau*pi/180.0
        if (tau_rad >= pi) then 
            tau_rad = tau_rad - 2*pi
        end if 

        ! define size 2D domain from 3D dimension and angle
        ! ---------------------------------------------------------------------
        if (tau_rad == 0) then
            mesh%Lx = simu%Lx1
        elseif (tau_rad == -pi) then
            mesh%Lx = simu%Lx2
        elseif (tau_rad == pi/2) then
            mesh%Lx = simu%Ly1     
        elseif (tau_rad == -pi/2) then
            mesh%Lx = simu%Ly2
        elseif (-pi/2<tau_rad .and.  tau_rad<0) then
            mesh%Lx = min(simu%Lx1/abs(cos(tau_rad)),simu%Ly2/abs(sin(tau_rad)))
        elseif (0<tau_rad .and.  tau_rad<pi/2) then
            mesh%Lx = min(simu%Lx1/abs(cos(tau_rad)),simu%Ly1/abs(sin(tau_rad)))
        elseif (pi/2<tau_rad .and.  tau_rad<pi) then
            mesh%Lx = min(simu%Lx2/abs(cos(tau_rad)),simu%Ly1/abs(sin(tau_rad)))
        elseif (-pi<tau_rad .and.  tau_rad<-pi/2) then
            mesh%Lx = min(simu%Lx2/abs(cos(tau_rad)),simu%Ly2/abs(sin(tau_rad)))
        end if

        ! create X vector
        ! ---------------------------------------------------------------------
        mesh%Lx_coarse = mesh%Lx+2*simu%domain_step
        mesh%nx_coarse = floor(mesh%Lx_coarse/output%dout)
        mesh%dx_coarse = output%dout

        if (allocated(mesh%x_coarse)) deallocate(mesh%x_coarse) 
        allocate(mesh%x_coarse(mesh%nx_coarse))
        do ii = 1, mesh%nx_coarse
        mesh%x_coarse(ii) = output%dout*(ii -1)
        enddo

        ! create Terrain elevation 
        ! ---------------------------------------------------------------------
        ! compute terrain height 
        if (terrain2D%ext .eqv. .true.) then 
            print*, "Interpolate terrain elevation from file"
            ! Create a 2D interpolate (must be larger that x_coarse to be able to find the first 
            ! receiver position
            if (allocated(x_temp)) deallocate(x_temp)
            if (allocated(mesh%h_coarse)) deallocate(mesh%h_coarse) 
            allocate(x_temp(floor((mesh%Lx_coarse + 100)/output%dout)))
            do ii=1, floor((mesh%Lx_coarse + 100.0)/output%dout)
                x_temp(ii) =  output%dout * (ii-1) - 50.0
            enddo 
            print*, terrain2D%Xs, terrain2D%Ys
            call interpolate_from_2D_to_1D(terrain2D%h, terrain2D%x, terrain2D%y,  &
                mesh%h_coarse, terrain2D%Xs, terrain2D%Ys, & 
                x_temp, tau_rad)

            if (allocated(mesh%hKnot)) deallocate(mesh%hKnot) 
            if (allocated(mesh%htx)) deallocate(mesh%htx) 
            call create_1D_interpoland(mesh%h_coarse, x_temp, mesh%hKnot, mesh%htx)
            deallocate(x_temp)
            deallocate(mesh%h_coarse)
            !print*, mesh%x_coarse
            ! interpolate h on coarse mesh 
            allocate(mesh%h_coarse(mesh%nx_coarse))
            call interpolate_from_2D_to_1D(terrain2D%h, terrain2D%x, terrain2D%y,  &
                mesh%h_coarse, terrain2D%Xs, terrain2D%Ys, & 
                mesh%x_coarse, tau_rad)
            

        else if (terrain2D%gaussian_hill .eqv. .true.) then 
            if (allocated(mesh%h_coarse)) deallocate(mesh%h_coarse) 
            allocate(mesh%h_coarse(mesh%nx_coarse))
            mesh%h_coarse = prof_h(mesh%x_coarse, terrain2D%hh, terrain2d%hl, tau_rad, terrain2D%hpos)
            if (allocated(mesh%hKnot)) deallocate(mesh%hKnot) 
            if (allocated(mesh%htx)) deallocate(mesh%htx) 
            call create_1D_interpoland(mesh%h_coarse, mesh%x_coarse, mesh%hKnot, mesh%htx)
        else
            if (allocated(mesh%h_coarse)) deallocate(mesh%h_coarse) 
            allocate(mesh%h_coarse(mesh%nx_coarse))
            mesh%h_coarse = 0.0
            if (allocated(mesh%hKnot)) deallocate(mesh%hKnot) 
            if (allocated(mesh%htx)) deallocate(mesh%htx) 
            call create_1D_interpoland(mesh%h_coarse, mesh%x_coarse, mesh%hKnot, mesh%htx)
        endif 


        ! create Z vector 
        ! ---------------------------------------------------------------------
        ! add heigh of the terrain (warning if terrain goes bellow zero this is not accounted for)
        f_min = minval(simu%frequencies) 
        L_pml = pml%size * c0 / (simu%dx_frac * f_min)

        z_min = minval(mesh%h_coarse)
        z_max = maxval(mesh%h_coarse) + simu%Lz + L_pml

        ! the total height of the domain may become higher thet zmax-zmin as the domain are tilted.
        ! the maximum height should be simu%domain step (tilted by 90° which shoul never happened)
        mesh%Lz_coarse = z_max - z_min  + simu%domain_step
        
        mesh%nz_coarse =  floor(mesh%Lz_coarse/output%dout) + 2
        if (allocated(mesh%z_coarse)) deallocate(mesh%z_coarse) 
        allocate(mesh%z_coarse(mesh%nz_coarse))
        do ii = 1, mesh%nz_coarse
            mesh%z_coarse(ii) = output%dout*(ii-1) + z_min
        enddo

        ! Create coarse mesh for Side interpolation
        ! ---------------------------------------------------------------------
        allocate(mesh%x_coarse_grid(mesh%nx_coarse, mesh%nz_coarse))
        allocate(mesh%z_coarse_grid(mesh%nx_coarse, mesh%nz_coarse))

        do ii = 1, mesh%nz_coarse 
            mesh%x_coarse_grid(:,ii) = mesh%x_coarse
        enddo

        do ii = 1, mesh%nx_coarse
            mesh%z_coarse_grid(ii,:) =  mesh%z_coarse
        enddo 
        
        ! Create grid for receiver interpolation
        ! ---------------------------------------------------------------------
        allocate(mesh%x_rec_grid(mesh%nx_coarse, output%nb_receiver))
        allocate(mesh%z_rec_grid(mesh%nx_coarse, output%nb_receiver))

        do ii = 1, output%nb_receiver 
            mesh%x_rec_grid(:,ii) = mesh%x_coarse
        enddo

        do ii = 1, mesh%nx_coarse 
            mesh%z_rec_grid(ii,:) = output%heights + mesh%h_coarse(ii)
        enddo
    end subroutine


    subroutine clean_coarse_mesh(terrain2D_ext)
        logical, intent(in) :: terrain2D_ext

        deallocate(mesh%x_rec_grid) 
        deallocate(mesh%z_rec_grid) 
        deallocate(mesh%x_coarse_grid) 
        deallocate(mesh%z_coarse_grid) 
        deallocate(mesh%x_coarse) 
        deallocate(mesh%z_coarse) 
        deallocate(mesh%h_coarse) 
        if (terrain2D_ext .eqv. .true.) then 
            deallocate(mesh%hKnot)
            deallocate(mesh%htx)
        endif 
    end subroutine 


    function compute_x_grid(ix, iz, idom) result(x)
        integer, intent(in)  :: ix, iz, idom
        real :: x
            x = mesh%xc(idom) + (ix -1) * mesh%dx * cos(mesh%ac(idom)) & 
                - (iz-1)*mesh%dz*sin(mesh%ac(idom))
    end function

    function compute_z_grid(ix, iz, idom) result(z)
        integer, intent(in)  :: ix, iz, idom
        real :: z
            z = mesh%zc(idom) + (ix -1) * mesh%dx * sin(mesh%ac(idom)) &
                +(iz-1)*mesh%dz*cos(mesh%ac(idom))
    end function



    real(kind=8) function local_slope_angle(x, hKnot, htx, dx)
        use types
        use interpoland_1D
        real(kind=8), intent(in)       :: x, hKnot(:), htx(:), dx
        real(kind=8) :: h_fwd, h_bwd

            h_fwd = interpolate(hKnot, htx, x + dx)
            h_bwd = interpolate(hKnot, htx, x - dx)

        local_slope_angle = atan2(h_fwd - h_bwd, 2.0_8 * dx)
    end function local_slope_angle


    subroutine compute_adaptive_xc(terrain2D, xc0, &
                                domain_step)
        use types 
        implicit none
    
        type(terrain_elevation), intent(in)     :: terrain2D
        real(kind=8), intent(in)               :: xc0
        real(kind=8), intent(in)               :: domain_step       ! target max step
    
        ! ── locals ────────────────────────────────────────────────────────────────
        integer,      parameter  :: max_domains = 100000
        real(kind=8), parameter  :: slope_dx    = 0.5_8   ! finite-diff step for slope [m]
        real(kind=8), parameter ::  min_domain_step=1, max_delta_angle=10*3.14/180
        real(kind=8) :: xc_tmp(max_domains+1)
        real(kind=8) :: x_cur, x_next, step
        real(kind=8) :: alpha_cur, alpha_next, delta_alpha, alpha
        integer      :: n, bisect_iter
    
        ! ── helper: local slope angle at position x ───────────────────────────────
        ! (reuse your existing terrain evaluation)
        
        ! initialize x position 
        x_cur    = xc0
        n        = 1
        xc_tmp(1)= xc0
        alpha_cur = local_slope_angle(x_cur, mesh%hKnot, mesh%htx, slope_dx)


        do while (x_cur < xc0 + mesh%Lx - min_domain_step)
    
            ! ── start from target step, bisect until angle change is acceptable ──
            step = min(domain_step, (xc0 + mesh%Lx) - x_cur)
    
            do bisect_iter = 1, 20          ! max 20 halvings → step/2^20 ~ negligible
                x_next     = x_cur + step
                alpha_next = local_slope_angle(x_next, mesh%hKnot, mesh%htx, slope_dx)
                delta_alpha = abs(alpha_next - alpha_cur)
    
                if (delta_alpha <= max_delta_angle) exit   ! step is acceptable
                if (step <= min_domain_step)        exit   ! cannot shrink further
    
                step = max(step * 0.5_8, min_domain_step)
            enddo
    
            ! enforce hard minimum
            step   = max(step, min_domain_step)
            x_cur  = min(x_cur + step, xc0 + mesh%Lx)
    
            n = n + 1
            if (n > max_domains) then
                print*, "WARNING: max_domains reached in adaptive xc placement"
                exit
            endif
            xc_tmp(n) = x_cur
            alpha_cur = alpha_next
        enddo
    
        ! ── ensure last point reaches end of domain ───────────────────────────────
        if (xc_tmp(n) < xc0 + mesh%Lx) then
            n = n + 1
            xc_tmp(n) = xc0 + mesh%Lx
        endif
    
        mesh%ndomain = n - 1
        allocate(mesh%xc(n))
        mesh%xc(1:n) = xc_tmp(1:n)
    
        print*, "Adaptive xc: ndomain =", mesh%ndomain, &
                "  mean step =", mesh%Lx/mesh%ndomain, &
                "  min step  =", min_domain_step, &
                "  max delta angle =", max_delta_angle
    
    end subroutine compute_adaptive_xc


    !! @brief Create the computing mesh for a given frequency and angle
    !!
    !! @details define a dx given the frequency and 'dx_frac', create the coresponding mesh
    !!
    !!
    !! @param [in]      freq: frequency of the calculation
    !! @param [in]      c0 : sound sped
    !! @param [in]      simu: parameters of the global simulation
    !! @param[inout]    mesh: the 2D mesh coresponding to the computed angle
    !! @param[inout]    mesh: parameters of the 2D mesh coresponding to the angle computed
    !! @return
    subroutine init_mesh(tau,freq,c0,simu,terrain2D,pml,dout)
        use types
        use utils
        ! use interpo
        use interpoland_1D
        use geometry_utils
        implicit none

        type(simulation), intent(in)                            :: simu
        type(terrain_elevation), intent(in)                     :: terrain2D
        real, intent(in)                                        :: tau,freq,c0
        type(perfectly_matched_layer), intent(in)               :: pml
        real, intent(in)                                        :: dout

        integer                                                 :: ii, jj, kk
        real                                                    :: pi = 3.141592653589793, tau_rad, ax, az, bx, bz, xc0
        real :: Lzc_raw, Lz_overlap, Lzc_max
        logical                                                 :: flag

        ! conversion in rad
        tau_rad = tau*pi/180
        if (tau_rad >= pi) then 
            tau_rad = tau_rad - 2*pi
        end if 

        ! define solving mesh quantities
        mesh%dx = c0/(simu%dx_frac*freq)
        if (mesh%dx> dout/5) then
            mesh%dx = dout/5
        endif
        mesh%dz = mesh%dx

        ! add pml on top of the domain
        mesh%Lz = simu%Lz + pml%size*mesh%dz

        if (simu%multi_domain) then 
            print*, "Initialize with multi domain"
            ! compute number of domain 
            mesh%ndomain = floor(mesh%Lx/simu%domain_step) + 1 

            ! compute characteristic of each domain 
            allocate(mesh%xc(mesh%ndomain+1)) ! to compute the last angle 
            allocate(mesh%zc(mesh%ndomain+1))
            allocate(mesh%Lxc(mesh%ndomain))
            allocate(mesh%Lzc(mesh%ndomain))
            allocate(mesh%nxc(mesh%ndomain))
            allocate(mesh%nzc(mesh%ndomain))
            allocate(mesh%ac(mesh%ndomain))

            mesh%source_height = simu%source_height
            mesh%source_height_abs = mesh%source_height + interpolate(mesh%hknot, mesh%htx, 0.0)
            call find_first_xc(mesh%x_coarse, mesh%hKnot, mesh%htx, mesh%nx_coarse, simu%source_height, xc0, flag)

            if (.not. flag) then 
                print*, "did not find first position for multi domain, set to 0.0"
                xc0 = 0.0
            endif 
            ! xc0 = 0.0
            print*, "first position at ", xc0

            ! set th initial  x coordinate of each sub domain 
            do ii = 1, mesh%ndomain+1
            mesh%xc(ii) = simu%domain_step * (ii-1) + xc0
            enddo

            ! set the initial  z coordinate  of each subdomain 
            if (terrain2D%ext .eqv. .true.) then 
                mesh%zc =  interpolate(mesh%hKnot, mesh%htx, mesh%xc)
            else if (terrain2D%gaussian_hill .eqv. .true.) then 
                mesh%zc = prof_h(mesh%xc,terrain2D%hh,terrain2D%hl,tau_rad,terrain2D%hpos)
            else
                mesh%zc=0.0
            end if

            ! compute angle and length of each domain 
            do ii = 1, mesh%ndomain
            mesh%ac(ii) = atan((mesh%zc(ii+1) - mesh%zc(ii))/simu%domain_step)
            mesh%Lxc(ii) = simu%domain_step/cos(mesh%ac(ii))
            enddo
            

            ! computation of domain height and 
            ! correction of domain length for overlapping
            ! see thesis from arthur guibard 
            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ! ii=1
            ! mesh%Lzc(ii) = mesh%Lz ! - mesh%zc(ii)
            ! if (mesh%ac(ii+1)<mesh%ac(ii)) then 
            !     mesh%Lxc(ii) = mesh%Lxc(ii) + tan(mesh%ac(ii) - mesh%ac(ii+1))*mesh%Lzc(ii)
            ! end if 
            ! mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            ! mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1
            !
            ! do ii = 2,mesh%ndomain-1
            !
            ! ! new method 
            !  if (tan(mesh%ac(ii) - mesh%ac(ii-1))<(mesh%Lxc(ii-1)/mesh%Lzc(ii-1))) then 
            !      mesh%Lzc(ii) = min(mesh%Lzc(ii-1)/cos(mesh%ac(ii-1) - mesh%ac(ii)),200.0)
            !  else
            !     ! angle to large need to find domain height from intersetion with previous sub domain 
            !     flag = .false.
            !     jj = ii - 1
            !     do while ((.not. flag) .and. (jj>0))
            !         ax = mesh%xc(jj) - (mesh%Lzc(jj)-pml%size*mesh%dz) * sin(mesh%ac(jj))
            !         az = mesh%zc(jj) + (mesh%Lzc(jj)-pml%size*mesh%dz) * cos(mesh%ac(jj))
            !         bx = mesh%xc(jj) + mesh%Lxc(jj) * cos(mesh%ac(jj)) - (mesh%Lzc(jj)-pml%size*mesh%dz) * sin(mesh%ac(jj))
            !         bz = mesh%zc(jj) + mesh%Lxc(jj) * sin(mesh%ac(jj)) + (mesh%Lzc(jj)-pml%size*mesh%dz) * cos(mesh%ac(jj))
            !         call find_height_to_intersect(mesh%xc(ii), mesh%zc(ii), &
            !             ax, az, bx, bz, mesh%ac(ii) + 0.5*pi, mesh%Lzc(ii), flag)
            !         jj = jj -1
            !     enddo
            !     mesh%Lzc(ii) = min(200.0, mesh%Lzc(ii) + pml%size * mesh%dz)
            !
            ! endif
            !
            ! if (mesh%ac(ii+1)<=mesh%ac(ii)) then 
            !     mesh%Lxc(ii) = mesh%Lxc(ii) + tan(mesh%ac(ii) - mesh%ac(ii+1))*mesh%Lzc(ii)
            ! endif 




            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            ! ═══════════════════════════════════════════════════════════════════
            ! PASS 1 – compute all Lzc first (with clamping), no Lxc overlap yet
            ! ═══════════════════════════════════════════════════════════════════
            Lzc_max = 200
            ii = 1
            mesh%Lzc(ii) = mesh%Lz
            
            do ii = 2, mesh%ndomain-1
                if (tan(mesh%ac(ii) - mesh%ac(ii-1)) < (mesh%Lxc(ii-1)/mesh%Lzc(ii-1))) then
                    Lzc_raw = mesh%Lzc(ii-1) / cos(mesh%ac(ii-1) - mesh%ac(ii))
                else
                    flag    = .false.
                    jj      = ii - 1
                    Lzc_raw = Lzc_max
                    do while ((.not. flag) .and. (jj > 0))
                        ax = mesh%xc(jj) - (mesh%Lzc(jj) - pml%size*mesh%dz) * sin(mesh%ac(jj))
                        az = mesh%zc(jj) + (mesh%Lzc(jj) - pml%size*mesh%dz) * cos(mesh%ac(jj))
                        bx = mesh%xc(jj) + mesh%Lxc(jj)*cos(mesh%ac(jj)) &
                                         - (mesh%Lzc(jj) - pml%size*mesh%dz)*sin(mesh%ac(jj))
                        bz = mesh%zc(jj) + mesh%Lxc(jj)*sin(mesh%ac(jj)) &
                                         + (mesh%Lzc(jj) - pml%size*mesh%dz)*cos(mesh%ac(jj))
                        call find_height_to_intersect(mesh%xc(ii), mesh%zc(ii), &
                            ax, az, bx, bz, mesh%ac(ii) + 0.5_8*pi, Lzc_raw, flag)
                        jj = jj - 1
                    enddo
                    Lzc_raw = Lzc_raw + pml%size * mesh%dz
                endif
                mesh%Lzc(ii) = min(Lzc_raw, Lzc_max)   ! clamp here
            enddo
            
            ii = mesh%ndomain
            Lzc_raw       = mesh%Lzc(ii-1) / cos(mesh%ac(ii-1) - mesh%ac(ii))
            mesh%Lzc(ii)  = min(Lzc_raw, Lzc_max)
            
            ! ═══════════════════════════════════════════════════════════════════
            ! PASS 2 – now compute Lxc using both Lzc(ii) AND Lzc(ii+1)
            ! ═══════════════════════════════════════════════════════════════════
            ii = 1
            if (mesh%ac(2) <= mesh%ac(1)) then
                Lz_overlap    = min(mesh%Lzc(1), mesh%Lzc(2))          ! ← key fix
                mesh%Lxc(1)   = mesh%Lxc(1) + tan(mesh%ac(1) - mesh%ac(2)) * Lz_overlap
            endif
            mesh%nxc(1) = floor(mesh%Lxc(1)/mesh%dx) + 3
            mesh%nzc(1) = floor(mesh%Lzc(1)/mesh%dz) + 1
            
            do ii = 2, mesh%ndomain-1
                ! forward overlap: only cover up to what the NEXT domain actually needs
                if (mesh%ac(ii+1) <= mesh%ac(ii)) then
                    Lz_overlap  = min(mesh%Lzc(ii), mesh%Lzc(ii+1))    ! ← key fix
                    mesh%Lxc(ii) = mesh%Lxc(ii) + tan(mesh%ac(ii) - mesh%ac(ii+1)) * Lz_overlap
                endif
                mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
                mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1
            enddo
            
            ii = mesh%ndomain
            mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1




            ! mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            ! mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1
            ! enddo 
            
            ! print*, mesh%nzc
            ! print*, '---------'
            ii = mesh%ndomain     
            mesh%Lzc(ii) = mesh%Lzc(ii-1)/cos(mesh%ac(ii-1) - mesh%ac(ii))
            mesh%Lxc(ii) = mesh%Lxc(ii)
            mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1

            print*, "initialized multi domain"
            print*, "max Lxc=", maxval(mesh%Lxc)
            print*, "max Lzc=", maxval(mesh%Lzc)
            mesh%nxmax = maxval(mesh%nxc)
            mesh%nzmax = maxval(mesh%nzc)
            mesh%nx = mesh%nxmax
            mesh%nz = mesh%nzmax

            ! allocate(mesh%x_grid(mesh%nxmax, mesh%nzmax, mesh%ndomain))
            ! allocate(mesh%z_grid(mesh%nxmax, mesh%nzmax,mesh%ndomain))
            allocate(mesh%x(mesh%nxmax))
            allocate(mesh%z(mesh%nzmax))

            ! local coordinate
            do ii = 1, mesh%nxmax
            mesh%x(ii) = mesh%dx*(ii -1)
            enddo

            do ii = 1, mesh%nzmax
            mesh%z(ii) = mesh%dz*(ii -1)
            enddo

            ! global coordinate
            ! do ii =1,mesh%ndomain 
            ! do kk=1,mesh%nzmax
            ! do jj=1,mesh%nxmax
            ! mesh%x_grid(jj,kk,ii)=mesh%xc(ii) &
            !     + (jj-1)*mesh%dx*cos(mesh%ac(ii)) - (kk-1)*mesh%dz*sin(mesh%ac(ii))
            !
            ! mesh%z_grid(jj,kk,ii)=mesh%zc(ii) &
            !     + (jj-1)*mesh%dx*sin(mesh%ac(ii)) + (kk-1)*mesh%dz*cos(mesh%ac(ii))
            ! enddo
            ! enddo 
            ! enddo 

            ! set the connectivty matrix 
            allocate(mesh%connectivity(mesh%ndomain, mesh%ndomain))
            mesh%connectivity = .false.
            ! loop over the domains
            do ii=1,mesh%ndomain-1
            do jj=ii+1, mesh%ndomain
            ! check if start of domain jj interect domain ii
            ! if (mesh%x_grid(1,mesh%nzc(jj),jj) < max(mesh%x_grid(mesh%nxc(ii),1,ii), &
            !     mesh%x_grid(mesh%nxc(ii),mesh%nzc(ii),ii))) then
            !     mesh%connectivity(ii,jj) = .true. 
            ! end if 
            if (compute_x_grid(1,mesh%nzc(jj),jj) < max(compute_x_grid(mesh%nxc(ii),1,ii), &
                compute_x_grid(mesh%nxc(ii),mesh%nzc(ii),ii))) then
                mesh%connectivity(ii,jj) = .true. 
            end if 
            enddo
            enddo

        elseif (simu%pe_model == 5) then 
            print*, "Initialize with curvilinear mesh"
            mesh%Lz = mesh%Lz_coarse 
            mesh%ndomain = 1 
            mesh%nx = floor(mesh%Lx/mesh%dx)+1
            mesh%nz = floor(mesh%Lz/mesh%dz)+1

            ! init mesh
            allocate(mesh%x(mesh%nx))
            allocate(mesh%z(mesh%nz))
            do ii = 1, mesh%nx
            mesh%x(ii) = mesh%dx*(ii -1)
            enddo

            do ii = 1, mesh%nz
            mesh%z(ii) = mesh%dz*(ii -1)
            enddo

            ! Curvilinear mesh for GTPE 
            allocate(mesh%h(mesh%nx))
            allocate(mesh%hp(mesh%nx))
            allocate(mesh%hpp(mesh%nx))

            if (terrain2D%ext .eqv. .true.) then 
                mesh%h = interpolate(mesh%hKnot, mesh%htx, mesh%x)
                do ii= 2, mesh%nx-1
                mesh%hp(ii) = (mesh%h(ii+1) - mesh%h(ii))/mesh%dx
                mesh%hpp(ii) = (mesh%h(ii+1) -2*mesh%h(ii) +mesh%h(ii-1))/(mesh%dx**2)
                enddo 
                mesh%hp(1) = (mesh%h(2) - mesh%h(1))/mesh%dx
                mesh%hp(mesh%nx) = (mesh%h(mesh%nx) - mesh%h(mesh%nx-1))/mesh%dx

                mesh%hpp(1) = (mesh%h(3) -2*mesh%h(2) +mesh%h(1))/(mesh%dx**2)
                mesh%hpp(mesh%nx) = (mesh%h(mesh%nx) - 2*mesh%h(mesh%nx-1) + mesh%h(mesh%nx-2))/(mesh%dx**2)

            else if (terrain2D%gaussian_hill .eqv. .true.) then 
                call prof(mesh%x,terrain2D%hh,terrain2D%hl,tau_rad,terrain2D%hpos, &
                    mesh%h,mesh%hp,mesh%hpp)
            else
                mesh%h = 0.0
                mesh%hp = 0.0
                mesh%hpp = 0.0
            endif

            allocate(mesh%x_grid(mesh%nx, mesh%nz))
            allocate(mesh%z_grid(mesh%nx, mesh%nz))
            do ii = 1, mesh%nz 
            mesh%x_grid(:, ii) = mesh%x
            enddo
            do ii = 1, mesh%nx
            mesh%z_grid(ii,:) = mesh%h(ii) + mesh%z*(mesh%Lz - (mesh%h(ii) - minval(mesh%h)))/mesh%Lz
            enddo 
        else 
            mesh%ndomain = 1 

            ! compute characteristic of each domain 
            allocate(mesh%xc(mesh%ndomain)) ! to compute the last angle 
            allocate(mesh%zc(mesh%ndomain))
            allocate(mesh%Lxc(mesh%ndomain))
            allocate(mesh%Lzc(mesh%ndomain))
            allocate(mesh%nxc(mesh%ndomain))
            allocate(mesh%nzc(mesh%ndomain))
            allocate(mesh%ac(mesh%ndomain))

            ! set th initial  x coordinate of each sub domain 
            mesh%xc = 0.0
            mesh%zc = 0.0
            mesh%ac = 0
            mesh%Lxc = mesh%Lx
            mesh%Lzc = mesh%Lz
            
            mesh%nxc(1) = floor(mesh%Lxc(1)/mesh%dx) + 1 
            mesh%nzc(1) = floor(mesh%Lzc(1)/mesh%dz) + 1


            mesh%nxmax = maxval(mesh%nxc)
            mesh%nzmax = maxval(mesh%nzc)
            mesh%nx = mesh%nxmax
            mesh%nz = mesh%nzmax

            ! allocate(mesh%x_grid(mesh%nxmax, mesh%nzmax))
            ! allocate(mesh%z_grid(mesh%nxmax, mesh%nzmax))
            allocate(mesh%x(mesh%nxmax))
            allocate(mesh%z(mesh%nzmax))

            ! local coordinate
            do ii = 1, mesh%nxmax
                mesh%x(ii) = mesh%dx*(ii -1)
            enddo

            do ii = 1, mesh%nzmax
                mesh%z(ii) = mesh%dz*(ii -1)
                ! mesh%x_grid(:,ii) = mesh%x
            enddo

            do ii = 1, mesh%nxmax
                ! mesh%z_grid(ii,:) = mesh%z
            enddo

            ! set the connectivty matrix 
            allocate(mesh%connectivity(mesh%ndomain, mesh%ndomain))
            mesh%connectivity = .false.
        end if 
    end subroutine init_mesh

    subroutine init_mesh2(tau,freq,c0,simu,terrain2D,pml,dout)
        use types
        use utils
        ! use interpo
        use interpoland_1D
        use geometry_utils
        implicit none

        type(simulation), intent(in)                            :: simu
        type(terrain_elevation), intent(in)                     :: terrain2D
        real, intent(in)                                        :: tau,freq,c0
        type(perfectly_matched_layer), intent(in)               :: pml
        real, intent(in)                                        :: dout

        integer                                                 :: ii, jj, kk
        real                                                    :: pi = 3.141592653589793, tau_rad, ax, az, bx, bz, xc0
        real :: Lzc_raw, Lz_overlap, Lzc_max
        logical                                                 :: flag

        ! conversion in rad
        tau_rad = tau*pi/180
        if (tau_rad >= pi) then 
            tau_rad = tau_rad - 2*pi
        end if 

        ! define solving mesh quantities
        mesh%dx = c0/(simu%dx_frac*freq)
        if (mesh%dx> dout/5) then
            mesh%dx = dout/5
        endif
        mesh%dz = mesh%dx

        ! add pml on top of the domain
        mesh%Lz = simu%Lz + pml%size*mesh%dz

        if (simu%multi_domain) then 
            print*, "Initialize with multi domain"
            ! compute number of domain 
            mesh%ndomain = floor(mesh%Lx/simu%domain_step) + 1 

            mesh%source_height = simu%source_height
            mesh%source_height_abs = mesh%source_height + interpolate(mesh%hknot, mesh%htx, 0.0)
            call find_first_xc(mesh%x_coarse, mesh%hKnot, mesh%htx, mesh%nx_coarse, simu%source_height, xc0, flag)

            if (.not. flag) then 
                print*, "did not find first position for multi domain, set to 0.0"
                xc0 = 0.0
            endif 
            print*, "first position at ", xc0
            
            call compute_adaptive_xc(terrain2D, xc0, simu%domain_step)

            allocate(mesh%zc(mesh%ndomain+1))
            allocate(mesh%Lxc(mesh%ndomain))
            allocate(mesh%Lzc(mesh%ndomain))
            allocate(mesh%nxc(mesh%ndomain))
            allocate(mesh%nzc(mesh%ndomain))
            allocate(mesh%ac(mesh%ndomain))

            ! set the initial  z coordinate  of each subdomain 
            if (terrain2D%ext .eqv. .true.) then 
                mesh%zc =  interpolate(mesh%hKnot, mesh%htx, mesh%xc)
            else if (terrain2D%gaussian_hill .eqv. .true.) then 
                mesh%zc = prof_h(mesh%xc,terrain2D%hh,terrain2D%hl,tau_rad,terrain2D%hpos)
            else
                mesh%zc=0.0
            end if

            ! compute angle and length of each domain 
            do ii = 1, mesh%ndomain
            mesh%ac(ii) = atan((mesh%zc(ii+1) - mesh%zc(ii))/(mesh%xc(ii+1) - mesh%xc(ii)))
            mesh%Lxc(ii) = simu%domain_step/cos(mesh%ac(ii))
            enddo
            

            ! computation of domain height and 
            ! correction of domain length for overlapping
            ! see thesis from arthur guibard 
            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            Lzc_max = 300
            ii=1
            mesh%Lzc(ii) = mesh%Lz ! - mesh%zc(ii)
            if (mesh%ac(ii+1)<mesh%ac(ii)) then 
                mesh%Lxc(ii) = mesh%Lxc(ii) + tan(mesh%ac(ii) - mesh%ac(ii+1))*mesh%Lzc(ii)
            end if 
            mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1

            do ii = 2,mesh%ndomain-1

            ! new method 
             if (tan(mesh%ac(ii) - mesh%ac(ii-1))<(mesh%Lxc(ii-1)/mesh%Lzc(ii-1))) then 
                 mesh%Lzc(ii) = min(Lzc_max,mesh%Lzc(ii-1)/cos(mesh%ac(ii-1) - mesh%ac(ii)))
             else
                ! angle to large need to find domain height from intersetion with previous sub domain 
                flag = .false.
                jj = ii - 1
                do while ((.not. flag) .and. (jj>0))
                    ax = mesh%xc(jj) - (mesh%Lzc(jj)-pml%size*mesh%dz) * sin(mesh%ac(jj))
                    az = mesh%zc(jj) + (mesh%Lzc(jj)-pml%size*mesh%dz) * cos(mesh%ac(jj))
                    bx = mesh%xc(jj) + mesh%Lxc(jj) * cos(mesh%ac(jj)) - (mesh%Lzc(jj)-pml%size*mesh%dz) * sin(mesh%ac(jj))
                    bz = mesh%zc(jj) + mesh%Lxc(jj) * sin(mesh%ac(jj)) + (mesh%Lzc(jj)-pml%size*mesh%dz) * cos(mesh%ac(jj))
                    call find_height_to_intersect(mesh%xc(ii), mesh%zc(ii), &
                        ax, az, bx, bz, mesh%ac(ii) + 0.5*pi, mesh%Lzc(ii), flag)
                    jj = jj -1
                enddo
                ! mesh%Lzc(ii) = min(200.0, mesh%Lzc(ii) + pml%size * mesh%dz)
                mesh%Lzc(ii) = min(Lzc_max,mesh%Lzc(ii) + pml%size * mesh%dz)

            endif

            if (mesh%ac(ii+1)<=mesh%ac(ii)) then 
                mesh%Lxc(ii) = mesh%Lxc(ii) + tan(mesh%ac(ii) - mesh%ac(ii+1))*mesh%Lzc(ii)
            endif 
            mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1
            enddo 
            ii = mesh%ndomain     
            mesh%Lzc(ii) = mesh%Lzc(ii-1)/cos(mesh%ac(ii-1) - mesh%ac(ii))
            mesh%Lxc(ii) = mesh%Lxc(ii)
            mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1




            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            ! ═══════════════════════════════════════════════════════════════════
            ! PASS 1 – compute all Lzc first (with clamping), no Lxc overlap yet
            ! ═══════════════════════════════════════════════════════════════════
            ! Lzc_max = 300
            ! ii = 1
            ! mesh%Lzc(ii) = mesh%Lz
            !
            ! do ii = 2, mesh%ndomain-1
            !     if (tan(mesh%ac(ii) - mesh%ac(ii-1)) < (mesh%Lxc(ii-1)/mesh%Lzc(ii-1))) then
            !         Lzc_raw = mesh%Lzc(ii-1) / cos(mesh%ac(ii-1) - mesh%ac(ii))
            !     else
            !         flag    = .false.
            !         jj      = ii - 1
            !         Lzc_raw = Lzc_max
            !         do while ((.not. flag) .and. (jj > 0))
            !             ax = mesh%xc(jj) - (mesh%Lzc(jj) - pml%size*mesh%dz) * sin(mesh%ac(jj))
            !             az = mesh%zc(jj) + (mesh%Lzc(jj) - pml%size*mesh%dz) * cos(mesh%ac(jj))
            !             bx = mesh%xc(jj) + mesh%Lxc(jj)*cos(mesh%ac(jj)) &
            !                              - (mesh%Lzc(jj) - pml%size*mesh%dz)*sin(mesh%ac(jj))
            !             bz = mesh%zc(jj) + mesh%Lxc(jj)*sin(mesh%ac(jj)) &
            !                              + (mesh%Lzc(jj) - pml%size*mesh%dz)*cos(mesh%ac(jj))
            !             call find_height_to_intersect(mesh%xc(ii), mesh%zc(ii), &
            !                 ax, az, bx, bz, mesh%ac(ii) + 0.5_8*pi, Lzc_raw, flag)
            !             jj = jj - 1
            !         enddo
            !         Lzc_raw = Lzc_raw + pml%size * mesh%dz
            !     endif
            !     mesh%Lzc(ii) = min(Lzc_raw, Lzc_max)   ! clamp here
            ! enddo
            !
            ! ii = mesh%ndomain
            ! Lzc_raw       = mesh%Lzc(ii-1) / cos(mesh%ac(ii-1) - mesh%ac(ii))
            ! mesh%Lzc(ii)  = min(Lzc_raw, Lzc_max)
            !
            ! ! ═══════════════════════════════════════════════════════════════════
            ! ! PASS 2 – now compute Lxc using both Lzc(ii) AND Lzc(ii+1)
            ! ! ═══════════════════════════════════════════════════════════════════
            ! ii = 1
            ! if (mesh%ac(2) <= mesh%ac(1)) then
            !     Lz_overlap    = min(mesh%Lzc(1), mesh%Lzc(2))          ! ← key fix
            !     mesh%Lxc(1)   = mesh%Lxc(1) + tan(mesh%ac(1) - mesh%ac(2)) * Lz_overlap
            ! endif
            ! mesh%nxc(1) = floor(mesh%Lxc(1)/mesh%dx) + 3
            ! mesh%nzc(1) = floor(mesh%Lzc(1)/mesh%dz) + 1
            !
            ! do ii = 2, mesh%ndomain-1
            !     ! forward overlap: only cover up to what the NEXT domain actually needs
            !     if (mesh%ac(ii+1) <= mesh%ac(ii)) then
            !         Lz_overlap  = min(mesh%Lzc(ii), mesh%Lzc(ii+1))    ! ← key fix
            !         mesh%Lxc(ii) = mesh%Lxc(ii) + tan(mesh%ac(ii) - mesh%ac(ii+1)) * Lz_overlap
            !     endif
            !     mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            !     mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1
            ! enddo
            !
            ! ii = mesh%ndomain
            ! mesh%nxc(ii) = floor(mesh%Lxc(ii)/mesh%dx) + 3
            ! mesh%nzc(ii) = floor(mesh%Lzc(ii)/mesh%dz) + 1


            
            print*, "initialized multi domain"
            print*, "max Lxc=", maxval(mesh%Lxc)
            print*, "max Lzc=", maxval(mesh%Lzc)

            mesh%nxmax = maxval(mesh%nxc)
            mesh%nzmax = maxval(mesh%nzc)
            mesh%nx = mesh%nxmax
            mesh%nz = mesh%nzmax

            ! allocate(mesh%x_grid(mesh%nxmax, mesh%nzmax, mesh%ndomain))
            ! allocate(mesh%z_grid(mesh%nxmax, mesh%nzmax,mesh%ndomain))
            allocate(mesh%x(mesh%nxmax))
            allocate(mesh%z(mesh%nzmax))

            ! local coordinate
            do ii = 1, mesh%nxmax
            mesh%x(ii) = mesh%dx*(ii -1)
            enddo

            do ii = 1, mesh%nzmax
            mesh%z(ii) = mesh%dz*(ii -1)
            enddo

            ! global coordinate
            ! do ii =1,mesh%ndomain 
            ! do kk=1,mesh%nzmax
            ! do jj=1,mesh%nxmax
            ! mesh%x_grid(jj,kk,ii)=mesh%xc(ii) &
            !     + (jj-1)*mesh%dx*cos(mesh%ac(ii)) - (kk-1)*mesh%dz*sin(mesh%ac(ii))
            !
            ! mesh%z_grid(jj,kk,ii)=mesh%zc(ii) &
            !     + (jj-1)*mesh%dx*sin(mesh%ac(ii)) + (kk-1)*mesh%dz*cos(mesh%ac(ii))
            ! enddo
            ! enddo 
            ! enddo 

            ! set the connectivty matrix 
            allocate(mesh%connectivity(mesh%ndomain, mesh%ndomain))
            mesh%connectivity = .false.
            ! loop over the domains
            do ii=1,mesh%ndomain-1
            do jj=ii+1, mesh%ndomain
            ! check if start of domain jj interect domain ii
            ! if (mesh%x_grid(1,mesh%nzc(jj),jj) < max(mesh%x_grid(mesh%nxc(ii),1,ii), &
            !     mesh%x_grid(mesh%nxc(ii),mesh%nzc(ii),ii))) then
            !     mesh%connectivity(ii,jj) = .true. 
            ! end if 
            if (compute_x_grid(1,mesh%nzc(jj),jj) < max(compute_x_grid(mesh%nxc(ii),1,ii), &
                compute_x_grid(mesh%nxc(ii),mesh%nzc(ii),ii))) then
                mesh%connectivity(ii,jj) = .true. 
            end if 
            enddo
            enddo

        elseif (simu%pe_model == 5) then 
            print*, "Initialize with curvilinear mesh"
            mesh%ndomain = 1 
            mesh%Lz = mesh%Lz_coarse
            mesh%nx = floor(mesh%Lx/mesh%dx)+1
            mesh%nz = floor(mesh%Lz/mesh%dz)+1

            ! init mesh
            allocate(mesh%x(mesh%nx))
            allocate(mesh%z(mesh%nz))
            do ii = 1, mesh%nx
            mesh%x(ii) = mesh%dx*(ii -1)
            enddo

            do ii = 1, mesh%nz
            mesh%z(ii) = mesh%dz*(ii -1)
            enddo

            ! Curvilinear mesh for GTPE 
            allocate(mesh%h(mesh%nx))
            allocate(mesh%hp(mesh%nx))
            allocate(mesh%hpp(mesh%nx))

            if (terrain2D%ext .eqv. .true.) then 
                mesh%h = interpolate(mesh%hKnot, mesh%htx, mesh%x)
                do ii= 2, mesh%nx-1
                mesh%hp(ii) = (mesh%h(ii+1) - mesh%h(ii))/mesh%dx
                mesh%hpp(ii) = (mesh%h(ii+1) -2*mesh%h(ii) +mesh%h(ii-1))/(mesh%dx**2)
                enddo 
                mesh%hp(1) = (mesh%h(2) - mesh%h(1))/mesh%dx
                mesh%hp(mesh%nx) = (mesh%h(mesh%nx) - mesh%h(mesh%nx-1))/mesh%dx

                mesh%hpp(1) = (mesh%h(3) -2*mesh%h(2) +mesh%h(1))/(mesh%dx**2)
                mesh%hpp(mesh%nx) = (mesh%h(mesh%nx) - 2*mesh%h(mesh%nx-1) + mesh%h(mesh%nx-2))/(mesh%dx**2)

            else if (terrain2D%gaussian_hill .eqv. .true.) then 
                call prof(mesh%x,terrain2D%hh,terrain2D%hl,tau_rad,terrain2D%hpos, &
                    mesh%h,mesh%hp,mesh%hpp)
            else
                mesh%h = 0.0
                mesh%hp = 0.0
                mesh%hpp = 0.0
            endif

            allocate(mesh%x_grid(mesh%nx, mesh%nz))
            allocate(mesh%z_grid(mesh%nx, mesh%nz))
            do ii = 1, mesh%nz 
            mesh%x_grid(:, ii) = mesh%x
            enddo
            do ii = 1, mesh%nx
            ! mesh%z_grid(ii,:) = mesh%h(ii) + mesh%z*(mesh%Lz - mesh%h(ii))/mesh%Lz
            mesh%z_grid(ii,:) = mesh%h(ii) + mesh%z*(mesh%Lz +minval(mesh%h)- (mesh%h(ii)))/(mesh%Lz + minval(mesh%h))
            enddo 
        else 
            mesh%ndomain = 1 

            ! compute characteristic of each domain 
            allocate(mesh%xc(mesh%ndomain)) ! to compute the last angle 
            allocate(mesh%zc(mesh%ndomain))
            allocate(mesh%Lxc(mesh%ndomain))
            allocate(mesh%Lzc(mesh%ndomain))
            allocate(mesh%nxc(mesh%ndomain))
            allocate(mesh%nzc(mesh%ndomain))
            allocate(mesh%ac(mesh%ndomain))

            ! set th initial  x coordinate of each sub domain 
            mesh%xc = 0.0
            mesh%zc = 0.0
            mesh%ac = 0
            mesh%Lxc = mesh%Lx
            mesh%Lzc = mesh%Lz
            
            mesh%nxc(1) = floor(mesh%Lxc(1)/mesh%dx) + 1 
            mesh%nzc(1) = floor(mesh%Lzc(1)/mesh%dz) + 1


            mesh%nxmax = maxval(mesh%nxc)
            mesh%nzmax = maxval(mesh%nzc)
            mesh%nx = mesh%nxmax
            mesh%nz = mesh%nzmax

            allocate(mesh%x_grid(mesh%nxmax, mesh%nzmax))
            allocate(mesh%z_grid(mesh%nxmax, mesh%nzmax))
            allocate(mesh%x(mesh%nxmax))
            allocate(mesh%z(mesh%nzmax))

            ! local coordinate
            do ii = 1, mesh%nxmax
                mesh%x(ii) = mesh%dx*(ii -1)
            enddo

            do ii = 1, mesh%nzmax
                mesh%z(ii) = mesh%dz*(ii -1)
                mesh%x_grid(:,ii) = mesh%x
            enddo

            do ii = 1, mesh%nxmax
                mesh%z_grid(ii,:) = mesh%z
            enddo

            ! set the connectivty matrix 
            allocate(mesh%connectivity(mesh%ndomain, mesh%ndomain))
            mesh%connectivity = .false.
        end if 
    end subroutine init_mesh2



    subroutine clean_mesh(simu_pe_model)
        integer, intent(in) :: simu_pe_model
        deallocate(mesh%x)
        deallocate(mesh%z)

        if (simu_pe_model == 5) then 
            deallocate(mesh%h)
            deallocate(mesh%hp)
            deallocate(mesh%hpp)
            deallocate(mesh%x_grid)
            deallocate(mesh%z_grid)
        else 
            deallocate(mesh%connectivity)
            deallocate(mesh%xc)
            deallocate(mesh%zc)
            deallocate(mesh%ac)
            deallocate(mesh%Lxc)
            deallocate(mesh%Lzc)
            deallocate(mesh%nxc)
            deallocate(mesh%nzc)
            if (allocated(mesh%x_grid)) deallocate(mesh%x_grid)
            if (allocated(mesh%z_grid)) deallocate(mesh%z_grid)
        end if 
    end subroutine 

    subroutine init_receiver(tau, output, sol) 
        use types
        implicit none
        ! input variables
        real, intent(in)                                    :: tau
        type(output_parameters), intent(inout)              :: output
        type(solution), intent(inout)                       :: sol 
        integer                                             :: jj


        ! Init vector for solution
        !-----------------------------------------------------------------------
        ! to save the full plane solution 
        ! check if plane coresponds to propagation angle
        do jj=1,size(output%xplane)
        if ((output%xplane(jj)==0) .and. ((tau==90) .or. (tau==-90) .or. (tau==270))) then 
            output%side = .true.
            print*, 'output%side set to true'
        endif 
        enddo

        do jj=1,size(output%yplane)
        if ((output%yplane(jj)==0) .and. ((tau==0) .or. (tau==180) .or. (tau==-180))) then 
            output%side = .true.
            print*, 'output%side set to true'
        endif 
        enddo 

        ! to save solution at a given position to reconstruct planes
        if (size(output%xplane)>0) then 
            allocate(sol%xplane(size(output%xplane),mesh%nz_coarse))
            allocate(sol%xcount(size(output%xplane)))
            allocate(sol%ycoord(size(output%xplane)))
            sol%xcount = 0
        endif

        if (size(output%yplane)>0) then 
            allocate(sol%yplane(size(output%yplane),mesh%nz_coarse))
            allocate(sol%ycount(size(output%yplane)))
            allocate(sol%xcoord(size(output%yplane)))
            sol%ycount=0
        endif 

        if (output%side) then
            allocate(sol%deltaL_side_intermediate(mesh%nx,mesh%nz_coarse))
            allocate(sol%deltaL_side(mesh%nx_coarse,mesh%nz_coarse))
            sol%deltaL_side = 2000000
            sol%deltaL_side_intermediate = 0
            ! sol%deltaL_int = 0
        endif

        ! to save solution on a line of receiver a z=cst 
        allocate(sol%deltaL_rec_intermediate(mesh%nx,output%nb_receiver))
        allocate(sol%deltaL_rec(mesh%nx_coarse,output%nb_receiver))
        sol%deltaL_rec_intermediate = 0
        sol%deltaL_rec = 0
    end subroutine init_receiver 

end module
