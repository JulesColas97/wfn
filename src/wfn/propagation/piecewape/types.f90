
!! Types module for definition of all the customs types of the progra=
module types

    type ragged_vector
        real,allocatable                    ::v(:)
    end type ragged_vector
    
    type ragged_array
        real,allocatable                    ::v(:,:)
    end type ragged_array

    ! Simulation global parameters
    type simulation
        character(len=40)                   :: case_name
        integer                             :: nb_tau                         ! number of angles
        integer                             :: nb_freq                          ! number of frequency
        real, allocatable,dimension (:)     :: tau                            ! angles
        real, allocatable,dimension (:)     :: frequencies                      ! frequencies
        integer                             :: pe_model
        logical                             :: continuation                     ! if set to .true. does not override existing frequency
        logical                             :: multi_domain
        real                                :: domain_step
        ! hill parameters (height, length, position) 
        real                                :: hh, hl, hpos

        real                                :: Lx1,Lx2,Ly,Ly1,Ly2,Lz                    ! 3D domain dimension
        real                                :: source_height
        real                                :: dx_frac                          ! number of points per wavelength (dx = lambda/dx_frac)

        logical                             :: external_flow                    ! read from h5
        logical                             :: interpolation                    ! interpolated from 3D field
        logical                             :: uniform                          ! uniform flow
        logical                             :: logarithmic                      ! logarithmic profil
        real                                :: u0                               ! unifrom value
    end type simulation

    ! 2D mesh parameters
    ! type domain
    !     integer                             :: nx                               ! number of point in x
    !     integer                             :: nz                               ! number of point in z
    !     real                                :: dx                               ! grid step
    !     real                                :: dz                               ! grid step
    !     real                                :: Lx                               ! x dimension
    !     real                                :: Lz                               ! z dimension
    !     real                                :: x0                               ! source position
    !     real                                :: z0                               ! source position
    !     ! parameter for recording mesh
    !     integer                             :: nx_coarse
    !     integer                             :: nz_coarse
    !     real                                :: Lx_coarse
    !     real                                :: Lz_coarse
    ! end type domain

    ! 2D mesh vector
    ! type mesh_2D
    !     real, allocatable, dimension(:)         :: x, x_coarse
    !     real, allocatable, dimension(:)         :: z, z_coarse
    !     real, allocatable, dimension(:,:)       :: x_coarse_grid
    !     real, allocatable, dimension(:,:)       :: z_coarse_grid
    !     real, allocatable, dimension(:,:)       :: x_rec_grid
    !     real, allocatable, dimension(:,:)       :: z_rec_grid
    !
    !     real, allocatable, dimension(:)         :: h, hp, hpp
    !     real, allocatable, dimension(:)         :: h_coarse
    !     real, allocatable,dimension(:)          :: hKnot, htx
    !
    !     real, allocatable, dimension(:,:,:)     :: x_grid
    !     real, allocatable, dimension(:,:,:)     :: z_grid
    !     real, allocatable, dimension(:)         :: xc,zc,ac,Lxc,Lzc
    !     logical, allocatable, dimension(:,:)    :: connectivity
    !     integer, allocatable, dimension(:)      :: nxc,nzc
    !     integer                                 :: nxmax, nzmax, ndomain
    ! end type mesh_2D
    
    type mesh_2D
        integer                             :: nx, nz, nx_coarse, nz_coarse
        real                                :: dx, dz, Lx, Lz, Lx_coarse, Lz_coarse, x0, z0, dx_coarse
        real :: source_height, source_height_abs

        real, allocatable, dimension(:)         :: x, x_coarse
        real, allocatable, dimension(:)         :: z, z_coarse
        real, allocatable, dimension(:,:)       :: x_coarse_grid
        real, allocatable, dimension(:,:)       :: z_coarse_grid
        real, allocatable, dimension(:,:)       :: x_rec_grid
        real, allocatable, dimension(:,:)       :: z_rec_grid
        
        ! terrain elevation 
        real, allocatable, dimension(:)         :: h, hp, hpp
        real, allocatable, dimension(:)         :: h_coarse
        real, allocatable,dimension(:)          :: hKnot, htx
        

        real, allocatable, dimension(:,:)     :: x_grid
        real, allocatable, dimension(:,:)     :: z_grid

        ! multi domain parameters 
        real, allocatable, dimension(:)         :: xc,zc,ac,Lxc,Lzc
        logical, allocatable, dimension(:,:)    :: connectivity
        integer, allocatable, dimension(:)      :: nxc,nzc
        integer                                 :: nxmax, nzmax, ndomain
    end type mesh_2D


    !impedance parameters
    type impedance
        integer                             :: mode
        real                                :: sigmae, alphae                         ! effective porosity
        real                                :: Rs, Om, q2, e                 ! effective resistivity
        complex                             :: Zs                               ! impedance
        complex                             :: beta                             ! admitance
        logical                             :: rigid                            ! if true beta=0
    end type impedance

    ! 3D fields for the mean flow
    type meanflow3D
        logical                             :: ext                              ! external file
        character(len=100)                  :: fname                            ! file name
        real, allocatable,dimension(:,:,:)  :: u, v,w, u_projected, c             ! 3D fields
        real, allocatable,dimension(:)      :: x, y, z                          ! coordinates
        integer                            :: nx, ny, nz                        ! coordinates
        real                                :: Xs, Ys, Zs, tau                ! source position in flow
    end type meanflow3D

    ! terrain elevation
    type terrain_elevation
        logical                             :: ext                              ! external file
        character(len=100)                  :: fname                            ! file name
        real, allocatable,dimension(:,:)    :: h               
        real                                :: Xs, Ys                          ! source position in the terrain
        real, allocatable,dimension(:)      :: x, y 
        integer                             :: nx, ny
        logical                             :: gaussian_hill
        real                                :: hh, hl,hpos
    end type terrain_elevation



    ! 2D fields for the mean flow
    type meanflow2D
        logical                             :: ext, uniform, logarithmic
        character(len=100)                  :: fdir                             ! directory for 2D data
        character(len=100)                  :: fname                            ! namle of h5 file
        real, allocatable, dimension (:)    :: utx, uty, ctx, cty, wtx, wty
        real, allocatable, dimension(:,:)   :: uKnot, wKnot
        real, allocatable, dimension(:,:)   :: cKnot
        real, allocatable,dimension(:,:)    :: u, v, w, c, epsilon, epsilon_eff, gamma, tau, m, n, n_eff ! flow quantities
        real, allocatable,dimension(:,:)    :: h10,h12,h20,h22 ! flow quantities

        real, allocatable,dimension(:,:)    :: x, z                             ! coarse grid
        real                                :: Xs, Ys, Zs                ! WT position
        real                                :: u0
    end type meanflow2D

    ! state varibales
    type variables_0
        real                                :: rho, T, P,c, gamma                          !masse volumique
    end type variables_0


    ! frequency quantities
    type frequency
        real                                :: freq                             ! frequency
        real                                :: omega                            ! pulsation
        real                                :: lambda                           ! wave length
        real                                :: k0                               ! wave number
    end type frequency

    ! surce term
    type source
        real                                :: pos_x                            ! position of acoustic source
        real                                :: pos_z
        real                                :: amp                              ! amplitude
    end type source

    ! pml parameters
    type perfectly_matched_layer
        integer                             :: size                             ! nb of points of the layer
        real                                :: param                            ! strength
        real                                :: n                                ! power
        real, allocatable, dimension(:)     :: sigma                            ! coeficient of the layer (along z)
        real, allocatable, dimension(:)     :: sigma_de                         ! first derivative
        
        ! for GTPE 
        integer                             :: iza,izb
        real, allocatable, dimension(:)     :: sigma_a, sigma_b
        real, allocatable, dimension(:)     :: sigmap_a,sigmap_b

    end type perfectly_matched_layer

    ! vector and matrix for the resolution
    type solve
        real, allocatable,dimension(:)      :: h10, h12, h20, h22             ! vector for arbitrati M

        ! vector for classic WAPE
        complex                             :: nu
        complex                             :: mu
        complex, allocatable,dimension(:)   :: nu_vect
        complex, allocatable,dimension(:)   :: mu_vect

        ! vector for tridiagonal matrix split step
        complex, allocatable,dimension(:)   :: a_vect
        complex, allocatable,dimension(:)   :: b_vect
        complex, allocatable,dimension(:)   :: c_vect

        complex, allocatable,dimension(:)   :: d_vect
        complex, allocatable,dimension(:)   :: e_vect
        complex, allocatable,dimension(:)   :: f_vect
    end type solve

    type solution
        complex, allocatable,dimension(:)       :: u
        complex, allocatable,dimension(:)       :: phi_temp,phi_m1, phi_p1
        complex, allocatable,dimension(:,:)   :: phi_vect, phi
        real, allocatable,dimension(:,:)      :: p_over_p0
        complex, allocatable,dimension(:,:)     :: p_vect
        real, allocatable,dimension(:,:)        :: deltaL_side, deltaL_side_intermediate
        real, allocatable,dimension(:,:)        :: deltaL_rec, deltaL_rec_intermediate


        ! real, allocatable,dimension(:,:)    :: deltaL_vect
        ! real, allocatable,dimension(:,:)    :: deltaL_int
        ! real,allocatable,dimension(:,:)     :: receiver
        ! real,allocatable,dimension(:,:)     :: receiver_int
        integer                             :: nbout

        ! varibale for plane x= cst 
        ! type(ragged_vector),allocatable     :: ycoord(:)
        ! type(ragged_array),allocatable      :: xplane(:)
        ! integer,allocatable,dimension(:)    :: xcount

        ! ! variable for plane y = cst 
        ! type(ragged_vector),allocatable     :: xcoord(:)
        ! type(ragged_array),allocatable      :: yplane(:)
        ! integer,allocatable,dimension(:)    :: ycount
        real, allocatable, DIMENSION(:,:)       :: xplane
        real, allocatable, dimension(:)         :: ycoord
        real, allocatable, dimension(:)         :: xcount

        real, allocatable, DIMENSION(:,:)       :: yplane
        real, allocatable, dimension(:)         :: xcoord
        real, allocatable, dimension(:)         :: ycount

    end type solution

    type output_parameters
        real                                :: dout
        integer                             :: nbout
        integer                             :: nb_receiver
        real, allocatable,dimension(:)      :: heights
        logical                             :: side
        logical                             :: write_flow
        logical                             :: gen_side
        logical                             :: top
        real,allocatable,dimension(:)       :: xplane    
        real,allocatable,dimension(:)       :: yplane    
    end type output_parameters


end module types
