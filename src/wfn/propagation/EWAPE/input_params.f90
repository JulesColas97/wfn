module input_params
    use kind_mod, only: wp
    use constants, only: PI, UNSET
    implicit none
    private
    public :: inputs_t, read_inputs

    integer, parameter :: MAX_N = 1000
    ! real(wp), parameter :: UNSET = -huge(1.0_wp)   ! sentinel

    type inputs_t
        character(len=256) :: case_name

        real(wp), allocatable :: frequencies(:)
        real(wp), allocatable :: tau(:)
        real(wp) :: taumin, taumax
        integer :: nb_tau
        real(wp) :: source_height , H
        real(wp) :: Lx1, Lx2, Ly1, Ly2, Lz

        ! physics 
        real(wp) :: c_air=340, rho_air=1.2, eta_air=0.0
        real(wp) :: c_water=1500, rho_water=1025, eta_water=0.0
        ! real(wp) :: rho_sb=1900, c_sb = 1650, beta_sb = 0
        real(wp) :: rho_sb, cp_sb, cs_sb, alpha_p, alpha_s, theta_sb
        complex(wp) :: beta_sb
        integer :: seabed_type
        real(wp) :: d_pml_bot

        ! numerical 
        real(wp) :: c_ref=850, dx_frac=2, dz_frac=20
        integer :: n_pade=10, n_pade_start = 10
        real(wp) :: theta_pade=PI/8
        integer :: pml_size = 30
        real(wp) :: pml_param = 5e4, pml_n = 2.5

        ! flow field 
        logical :: external_flow, uniform, logarithmic, power
        real(wp) :: u_ref, z_ref, alpha, z_roughness
        character(len=256) :: fname_flow
        real(wp) :: Xs_flow, Ys_flow


        ! surface 
        logical :: external_wave, sinus
        real(wp) :: A_wave, lambda_wave, phi_wave
        character(len=256) :: fname_wave
        real(wp) :: XS_wave, Ys_wave

        ! atmospheric absorption
        logical :: abs_air, abs_water
        real(wp) :: T0, p0, rh 

        ! output 
        real(wp) :: dx_out, dz_out 
        logical :: side
        real, allocatable :: receiver_heights(:) 
    end type


contains 


    subroutine allocate_and_fill_float(dummy_x, x, default)
        real(wp),dimension(:), intent(in) :: dummy_x
        real(wp), ALLOCATABLE, dimension(:), intent(out) :: x
        real(wp), intent(in) :: default
        integer :: n, ii 


        n = 0
        do ii = 1, MAX_N
        if (dummy_x(ii) == UNSET) exit
        n = ii
        end do

        if (n == 0) then                       ! nothing given: use defaults
            n = 1
            allocate(x(1))
            x(1) = default

        else 
            allocate(x(n))
            x = dummy_x(1:n)
        end if
    end subroutine 


    subroutine read_inputs(cfg)
        type(inputs_t), intent(out) :: cfg
        character(len=256) ::  filename="input.nml"
        integer :: num_args, ii
        character(len=256) :: arg
        logical :: file_exists

        ! local copy of the cfg variable for the name list 
        ! arrays are set of fix sized MAX_N 
        character(len=64) :: case_name = "c0"

        real(wp) :: frequencies(MAX_N)
        real(wp) :: tau(MAX_N)
        real(wp) :: taumin=UNSET, taumax=UNSET, dtau
        integer :: nb_tau=UNSET
        real(wp) :: source_height=119 , H=100
        real(wp) :: Lx1=300, Lx2=300, Ly1=300, Ly2=300, Lz=100

        ! physics 
        real(wp) :: c_air=340, rho_air=1.2, eta_air=0.0
        real(wp) :: c_water=1500, rho_water=1000, eta_water=0.0
        integer  :: seabed_type = 0
        real(wp) :: d_pml_bot = 2
        ! real(wp) :: rho_sb=1900, c_sb = 1650, beta_sb = UNSET
        real(wp) :: rho_sb=1900, cp_sb=1700, cs_sb=0.0, alpha_p=0.0, alpha_s=0.0,theta_sb=0.0
        complex(wp) :: beta_sb=UNSET

        ! numerical 
        real(wp) :: c_ref=850, dx_frac=2, dz_frac=20
        integer  :: n_pade=10, n_pade_start = 10
        real(wp) :: theta_pade=PI/8
        integer  :: pml_size = 30
        real(wp) :: pml_param = 5e4, pml_n = 2.5


        ! flow field 
        logical :: external_flow=.false.
        logical :: uniform=.false., logarithmic=.false., power=.false.
        real(wp) :: u_ref=0.0, z_ref=90, alpha=0.2, z_roughness=0.01

        character(len=256) :: fname_flow="./u.h5"
        real(wp) :: Xs_flow=0.0, Ys_flow=0.0
        ! 
        ! surface 
        logical :: external_wave=.false., sinus=.false.
        real(wp) :: A_wave=2, lambda_wave=80, phi_wave=0.0
        character(len=256) :: fname_wave="./wave.h5"
        real(wp) :: XS_wave=0.0, Ys_wave=0.0

        ! absorption 
        logical :: abs_air=.false., abs_water=.false.
        real(wp) :: T0 = 290, p0=101325 , rh=80

        ! output 
        real(wp) :: dx_out=1.0, dz_out=1.0
        logical :: side = .true.
        real(wp) :: receiver_heights(MAX_N)


        namelist/inputs/ case_name, frequencies, tau, taumin, taumax, nb_tau, source_height, H, &
            Lx1, Lx2, Ly1, Ly2, Lz,  & 
            c_air, rho_air, eta_air, c_water, rho_water, eta_water, &
            seabed_type, rho_sb, cp_sb,cs_sb,alpha_p, alpha_s, beta_sb, theta_sb,d_pml_bot, &
            c_ref, dx_frac, dz_frac, n_pade, n_pade_start, theta_pade, pml_size, pml_param, pml_n, &
            external_flow, uniform, logarithmic, power, u_ref, z_ref, alpha, z_roughness, &
            fname_flow, Xs_flow, Ys_flow, &
            external_wave, sinus, A_wave, lambda_wave, phi_wave, &
            fname_wave, Xs_wave, Ys_wave, &
            abs_air, abs_water, T0, p0, rh, &
            dx_out, dz_out, receiver_heights, side


        ! Default angle an frequencies 
        ! ---------------------------------------------------------------------
        tau = UNSET
        frequencies = UNSET
        receiver_heights = UNSET
        ! read file name from -f <arg>, default is input.nml
        ! ---------------------------------------------------------------------
        num_args = command_argument_count()
        do ii = 1, num_args
        call get_command_argument(ii, arg)

        ! Check for the -f flag
        if (trim(arg) == '-f') then
            if (ii < num_args) then
                call get_command_argument(ii + 1, filename)
            else
                print *, "Error: -f flag provided without a filename."
                stop
            end if
        end if
        end do

        print*, "Using input file: ", trim(filename)
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            print*, "File does not exists... stopping"
            stop
        endif

        open(unit=1,file=trim(filename))
        read(unit=1,nml=inputs)
        close(1)

        ! fill cfg with default vaues or read values 
        ! ---------------------------------------------------------------------
        cfg%case_name = case_name
        cfg%taumin = taumin 
        cfg%taumax = taumax 
        cfg%nb_tau = nb_tau 
        cfg%source_height = source_height
        cfg%H = H
        cfg%Lx1=Lx1
        cfg%Lx2 = Lx2
        cfg%Ly1 = Ly1
        cfg%Ly2 = Ly2
        cfg%Lz = Lz

        ! physics 
        cfg%c_air = c_air
        cfg%rho_air = rho_air
        cfg%eta_air = eta_air 
        cfg%c_water = c_water
        cfg%rho_water = rho_water
        cfg%eta_water = eta_water
        cfg%seabed_type = seabed_type
        cfg%d_pml_bot = d_pml_bot
        cfg%cp_sb = cp_sb 
        cfg%cs_sb = cs_sb 
        cfg%alpha_p = alpha_p 
        cfg%alpha_s = alpha_s
        cfg%theta_sb = theta_sb

        cfg%rho_sb = rho_sb
        cfg%beta_sb = beta_sb
        
        ! numerical 
        cfg%c_ref = c_ref 
        cfg%dx_frac = dx_frac
        cfg%dz_frac = dz_frac
        cfg%theta_pade = theta_pade
        cfg%n_pade = n_pade
        cfg%n_pade_start = n_pade_start
        cfg%pml_size = pml_size
        cfg%pml_param = pml_param
        cfg%pml_n = pml_n

        ! flow 
        cfg%external_flow = external_flow
        cfg%uniform = uniform
        cfg%logarithmic = logarithmic
        cfg%power = power
        cfg%u_ref = u_ref
        cfg%z_ref = z_ref
        cfg%alpha = alpha
        cfg%z_roughness = z_roughness
        cfg%fname_flow = fname_flow
        cfg%Xs_flow = Xs_flow
        cfg%Ys_flow = Ys_flow

        ! surface 
        cfg%external_wave = external_wave
        cfg%sinus = sinus
        cfg%A_wave = A_wave
        cfg%lambda_wave = lambda_wave
        cfg%phi_wave = phi_wave
        cfg%fname_wave = fname_wave
        cfg%Xs_wave = Xs_wave
        cfg%Ys_wave = Ys_wave


        ! atmospheric absorption
        cfg%abs_air = abs_air
        cfg%abs_water = abs_water
        cfg%T0 = T0
        cfg%p0 = p0
        cfg%rh = rh

        ! output 
        cfg%dx_out = dx_out 
        cfg%dz_out = dz_out 
        cfg%side = side

        ! Dealing with array 
        ! ---------------------------------------------------------------------
        call allocate_and_fill_float(frequencies, cfg%frequencies, default=100.0_wp)



        if (taumin/=UNSET .and. taumax/=Unset .and. nb_tau/=UNSET) then
            ! CASE A: Build from min/max/count
            print *, "Generating tau vector from taumin/taumax..."
            cfg%nb_tau = nb_tau
            allocate(cfg%tau(cfg%nb_tau))
            if (nb_tau == 1) then
                cfg%tau(1) = taumin
            else
                dtau = (taumax - taumin) / real(nb_tau - 1)
                do ii = 1, nb_tau
                cfg%tau(ii) = taumin + real(ii-1) * dtau
                end do
            end if
            print*, cfg%tau
        else
            call allocate_and_fill_float(tau, cfg%tau, default=0.0_wp)
        endif
        call allocate_and_fill_float(receiver_heights, cfg%receiver_heights, default=2.0_wp)
    end subroutine

end module input_params 
