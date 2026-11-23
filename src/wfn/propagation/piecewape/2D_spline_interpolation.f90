module interp_2D_spline
    implicit none

    ! Define the generic interface
    interface spline2D_interp
        module procedure interp_real_1d, interp_real_2d
        module procedure interp_cmplx_1d, interp_cmplx_2d
        module procedure interp_cmplx_from_int_1d, interp_cmplx_from_int_2D
    end interface

    public :: create_2D_complex_interpoland
    public :: interpolate_2D_complex

    private :: interp_core_real, interp_core_complex

contains
    ! --- 2D WRAPPERS ---
    subroutine interp_real_2d(values, x1, y1, x, y, res)
        real, intent(in)  :: values(:,:),x1(:), y1(:), x(:,:), y(:,:)
        real, intent(out) :: res(:,:)
        call interp_core_real(values, x1, y1, size(values,1), size(values,2), &
                               x, y, size(x), res)
    end subroutine

    subroutine interp_cmplx_2d(values, x1, y1, x, y, res)
        real, intent(in)     :: x1(:), y1(:), x(:,:), y(:,:)
        complex, intent(in)  :: values(:,:)
        complex, intent(out) :: res(:,:)
        call interp_core_complex(values, x1, y1, size(values,1), size(values,2), &
                               x, y, size(x), res)
    end subroutine

    ! --- 1D WRAPPERS ---
    subroutine interp_real_1d(values, x1, y1, x, y, res)
        real, intent(in)  :: values(:,:), x1(:), y1(:), x(:), y(:)
        real, intent(out) :: res(:)
        call interp_core_real(values, x1, y1, size(values,1), size(values,2), &
                               x, y, size(x), res)
    end subroutine

    subroutine interp_cmplx_1d(values, x1, y1, x, y, res)
        real, intent(in)     :: x1(:), y1(:), x(:), y(:)
        complex, intent(in)  :: values(:,:)
        complex, intent(out) :: res(:)
        call interp_core_complex(values, x1, y1, size(values,1), size(values,2), &
                               x, y, size(x), res)
    end subroutine

    subroutine interp_core_real(values, x1, y1, nx, ny,  x2, y2, n_pts, res)

        use bspline_sub_module
        integer, intent(in) :: nx, ny, n_pts
        real, intent(in) :: values(nx, ny)
        real, intent(in) :: x1(nx), y1(ny)
        real, intent(in) :: x2(n_pts), y2(n_pts)
        real, intent(out) :: res(n_pts)

        integer :: ii
        !         Declaration of locale variabales used for the interpolation
        !-------------------------------------------------------------------------------
        integer,parameter :: kx=3, ky=3 , iknot=0
        real, allocatable, dimension(:) :: tx,ty
        real,allocatable, dimension (:,:) :: fKnot

        real,dimension(ky) :: w1
        real,dimension(3*max(kx,ky)) :: w2
        integer :: inbvx,inbvy,iloy, idx,idy, iflag


        ! if (.not. allocated(res)) allocate(res(n_pts))
        allocate(tx(nx+kx))
        allocate(ty(ny+ky))
        allocate(fKnot(nx,ny))

        !have to set these before the first evaluate call:
        inbvx = 1
        inbvy = 1
        iloy  = 1
        idx = 0
        idy = 0
        ! create knots for interpolation
        call db2ink(x1,nx,y1,ny,real(values), &
            & kx,ky,iknot,tx(:),ty(:),fKnot,iflag)

        ! loop over the 2D points of the domain
        do ii = 1,n_pts
        ! interpolation of c
        if (x1(1)<x2(ii) .and. x2(ii)<x1(nx) &
            .and. y1(1)<y2(ii)  .and. y2(ii)<y1(ny)) then 
            call db2val(x2(ii),y2(ii),idx,idy,&
                tx(:),ty(:),nx,ny,kx,ky,&
                fKnot,res(ii) ,iflag,&
                inbvy,inbvy,iloy,&
                w1(:),w2(:),.true.)
        endif 
        enddo
    end subroutine

    ! --- Version for COMPLEX (Double Precision) ---
    subroutine interp_core_complex(values, x1, y1, nx, ny,  x2, y2, n_pts, res)
        use bspline_sub_module
        integer, intent(in) :: nx, ny, n_pts
        complex, intent(in) :: values(nx, ny)
        real, intent(in) :: x1(nx), y1(ny)
        real, intent(in) :: x2(n_pts), y2(n_pts)
        complex, intent(out) :: res(n_pts)

        integer :: ii
        !         Declaration of locale variabales used for the interpolation
        !-------------------------------------------------------------------------------
        integer,parameter :: kx=3, ky=3 , iknot=0
        real, allocatable, dimension(:) :: txi,tyi, txr, tyr
        real,allocatable, dimension (:,:) :: fKnoti, fKnotr
        ! real, allocatable,  dimension (:) :: f2i, f2r
        real :: f2i, f2r

        real,dimension(ky) :: w1
        real,dimension(3*max(kx,ky)) :: w2
        real,parameter ::  eps = 1.0e-12
        integer :: inbvx,inbvy,iloy, idx,idy, iflag


        ! if (.not. allocated(f2i)) allocate(f2i(n_pts))
        ! if (.not. allocated(f2r)) allocate(f2r(n_pts))
        allocate(txi(nx+kx))
        allocate(tyi(ny+ky))
        allocate(fKnoti(nx,ny))
        allocate(txr(nx+kx))
        allocate(tyr(ny+ky))
        allocate(fKnotr(nx,ny))

        !have to set these before the first evaluate call:
        inbvx = 1
        inbvy = 1
        iloy  = 1
        idx = 0
        idy = 0
        ! create knots for interpolation
        call db2ink(x1,nx,y1,ny,real(values), &
            & kx,ky,iknot,txr(:),tyr(:),fKnotr,iflag)
        
        ! create knots for interpolation
        call db2ink(x1,nx,y1,ny,aimag(values), &
            & kx,ky,iknot,txi(:),tyi(:),fKnoti,iflag)

        ! loop over the 2D points of the domain
        do ii = 1,n_pts
            if (x2(ii)<(x1(1)-eps) .or. x2(ii)>(x1(nx)+eps) &
                .or. y2(ii)<(y1(1)-eps)  .or. y2(ii)>(y1(ny)+eps)) then 
                cycle
            endif
            call db2val(x2(ii),y2(ii),idx,idy,&
                txr(:),tyr(:),nx,ny,kx,ky,&
                fKnotr,f2r,iflag,&
                inbvy,inbvy,iloy,&
                w1(:),w2(:),.true.)
            call db2val(x2(ii),y2(ii),idx,idy,&
                txi(:),tyi(:),nx,ny,kx,ky,&
                fKnoti,f2i,iflag,&
                inbvy,inbvy,iloy,&
                w1(:),w2(:),.true.)

            res(ii) = cmplx(f2r, f2i)
        enddo
    end subroutine


    subroutine create_2D_complex_interpoland(u1, x1, y1, uKnot,tx,ty)
        use bspline_sub_module
    
        implicit none
        ! input variables
        !-------------------------------------------------------------------------------
        complex, dimension (:,:),intent(in)                      :: u1
        real, dimension (:), intent(in)                       :: x1
        real, dimension (:), intent(in)                       :: y1
    
        complex,allocatable,  dimension (:,:),intent(inout)      :: uKnot
        complex,allocatable,  dimension (:),intent(inout)        :: tx,ty
    
        ! local variables
        integer                                 :: nx1,ny1
    
        !         Declaration of locale variabales used for the interpolation
        !-------------------------------------------------------------------------------
        integer,parameter                       :: kx    = 3     !! order in x
        integer,parameter                       :: ky    = 3     !! order in y
        integer,parameter                       :: iknot = 0
        INTEGER                                 :: iflag
        complex :: eye = (0.0,1.0)
        
        real, allocatable :: uKnot_temp(:,:), tx_temp(:), ty_temp(:)

        ! initialisation
        !-------------------------------------------------------------------------------
        nx1 = size(x1)
        ny1 = size(y1)
        
        allocate(uKnot_temp(nx1,ny1))
        allocate(tx_temp(nx1+kx))
        allocate(ty_temp(ny1+ky))
    
        allocate(uKnot(nx1,ny1))
        allocate(tx(nx1+kx))
        allocate(ty(ny1+ky))
    
        ! create knots for interpolation
        call db2ink(x1,nx1,y1,ny1,real(u1),kx,ky,iknot,tx_temp,ty_temp,uKnot_temp,iflag)
    
        uKnot =  uKnot_temp
        tx =  tx_temp
        ty =  ty_temp
        
        call db2ink(x1,nx1,y1,ny1,aimag(u1),kx,ky,iknot,tx_temp,ty_temp,uKnot_temp,iflag)
        uKnot = uKnot + eye * uKnot_temp
        tx = tx + eye * tx_temp
        ty = ty + eye * ty_temp
    end subroutine create_2D_complex_interpoland
    
    subroutine interp_cmplx_from_int_2d(uKnot,tx, ty, x, y, res)
        real, intent(in)  :: x(:,:), y(:,:)
        complex, intent(in)  :: tx(:), ty(:), uKnot(:,:)
        complex, intent(out) :: res(:,:)
        call interpolate_2D_complex(uKnot, tx, ty, size(uKnot,1), size(uKnot,2),  &
                               x, y, size(x), res)
    end subroutine

    subroutine interp_cmplx_from_int_1d(uKnot, tx, ty, x, y, res)
        real, intent(in)     :: x(:), y(:)
        complex, intent(in)  :: uKnot(:,:), tx(:), ty(:)
        complex, intent(out) :: res(:)
        call interpolate_2D_complex(uKnot, tx, ty, size(uKnot,1), size(uKnot,2), &
                               x, y, size(x), res)
    end subroutine
    
    subroutine interpolate_2D_complex(uKnot,tx, ty, nx1, ny1, x2, y2, n_pts, uval)
        use bspline_sub_module
    
        implicit none
        ! input variables
        !-------------------------------------------------------------------------------
    
        complex, intent(in)       :: uKnot(nx1,ny1)
        complex, intent(in)         :: tx(nx1+3),ty(ny1+3)
    
        integer, intent(in) :: n_pts
        real, intent(in) :: x2(n_pts)
        real, intent(in) :: y2(n_pts)
        complex,intent(out) :: uval(n_pts)
    
        ! local variables
        integer                                             :: nx1,ny1
        integer                                             :: nx2,ny2
        integer                                             :: ii,jj
    
        !         Declaration of locale variabales used for the interpolation
        !-------------------------------------------------------------------------------
        integer,parameter                       :: kx    = 3     !! order in x
        integer,parameter                       :: ky    = 3     !! order in y
        integer,parameter                       :: iknot = 0
        real,dimension(ky)                      :: w1
        real,dimension(3*max(kx,ky))            :: w0
        integer                                 :: inbvx,inbvy,inbvz,iloy, idx,idy,iflag
        real :: f2r, f2i
        real :: uKnotr(nx1,ny1), uKnotf(nx1,ny1)
    
        !have to set these before the first evaluate call:
        inbvx = 1
        inbvy = 1
        iloy  = 1
        idx = 0
        idy = 0
    
        uKnotr = real(uKnot)
        uKnotf = aimag(uKnot)
        ! evalutae function 
        do ii = 1,n_pts
        call db2val(x2(ii),y2(ii),idx,idy,real(tx),real(ty),nx1,ny1,kx,ky,uKnotr,f2r,iflag,inbvx,inbvy,iloy,w1,w0,.false.)
        call db2val(x2(ii),y2(ii),idx,idy,aimag(tx),aimag(ty),nx1,ny1,kx,ky,uKnotf,f2i,iflag,inbvx,inbvy,iloy,w1,w0,.false.)
        uval(ii) = cmplx(f2r, f2i)
        enddo 
    end subroutine interpolate_2D_complex

end module interp_2D_spline
