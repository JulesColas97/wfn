module interpoland_1D
    implicit none

    ! Define the generic interface
    interface create_interpoland
        module procedure create_1D_interpoland
        module procedure create_2D_interpoland
    end interface
    
    interface interpolate
        module procedure interpolate_1D_vector, interpolate_1D_value
        module procedure interpolate_2D_vector2D, interpolate_2D_vector1D, interpolate_2D_value
    end interface

    private :: create_1D_interpoland
    private :: interpolate_1D_vector
    private :: interpolate_1D_value

contains

subroutine create_1D_interpoland(val, x1, Knot,tx)
    use bspline_sub_module
    implicit none

    ! inout variables
    real, dimension(:), intent(inout)                         :: val
    real, dimension(:), intent(inout)                         :: x1
    real,allocatable,  dimension (:),intent(inout)            :: Knot
    real,allocatable,  dimension (:),intent(inout)            :: tx


    ! local variables
    integer                                                   :: kx    !! order in x
    integer                                                   :: iknot
    integer                                                   :: iflag, nx1
    integer                                                   :: inbvx, idx
    kx = 3
    iknot =0

    nx1 = size(x1)

    if (.not. allocated(Knot)) allocate(Knot(nx1))
    if (.not. allocated(tx)) allocate(tx(nx1+4))

    idx = 0
    inbvx = 1

    ! real part interpolation
    call db1ink(x1,nx1,val,kx,iknot,tx,Knot,iflag)
end subroutine create_1D_interpoland

subroutine create_2D_interpoland(u1, x1, y1, uKnot,tx,ty)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------
    real, dimension (:,:),intent(in)                      :: u1
    real, dimension (:), intent(in)                       :: x1
    real, dimension (:), intent(in)                       :: y1

    real,allocatable,  dimension (:,:),intent(inout)      :: uKnot
    real,allocatable,  dimension (:),intent(inout)        :: tx,ty

    ! local variables
    integer                                 :: nx1,ny1

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                       :: kx=3 !! order in x
    integer,parameter                       :: ky=3 !! order in y
    integer,parameter                       :: iknot=0
    INTEGER                                 :: iflag
    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(x1)
    ny1 = size(y1)

    allocate(uKnot(nx1,ny1))
    allocate(tx(nx1+kx))
    allocate(ty(ny1+ky))

    ! create knots for interpolation
    call db2ink(x1,nx1,y1,ny1,u1,kx,ky,iknot,tx,ty,uKnot,iflag)
end subroutine create_2D_interpoland



function interpolate_1D_vector(Knot,tx, x2) result(val2)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real, dimension (:),intent(in)        :: Knot
    real,  dimension (:),intent(in)       :: tx

    real, dimension (:),intent(in)        :: x2
    real,allocatable,  dimension (:)      :: val2

    ! local variables
    integer                                             :: nx1
    integer                                             :: nx2
    integer                                             :: ii

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                                   :: kx=3 !! order in x
    real,dimension(3*4)                                 :: w1
    integer                                             :: inbvx, idx, iflag

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(Knot)
    nx2 = size(x2)

    if(.not. allocated(val2)) allocate(val2(nx2))

    !have to set these before the first evaluate call:
    idx = 0
    inbvx = 1

    ! evalutae function 
    do ii = 1,nx2
        call db1val(x2(ii),idx,tx,nx1,kx,Knot,val2(ii),iflag,inbvx,w1)
    enddo
end function

function interpolate_1D_value(Knot,tx, x2) result(val2)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real, dimension (:),intent(in)         :: Knot
    real,  dimension (:),intent(in)         :: tx

    real, intent(in)                        :: x2
    real                                    :: val2

    ! local variables
    integer                                             :: nx1
    integer                                             :: nx2
    integer                                             :: ii

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                                   :: kx=3 !! order in x
    real,dimension(3*4)                                 :: w1
    integer                                             :: inbvx, idx, iflag

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(Knot)

    !have to set these before the first evaluate call:
    idx = 0
    inbvx = 1

    ! evalutae function 
    call db1val(x2,idx,tx,nx1,kx,Knot,val2,iflag,inbvx,w1)
end function

function interpolate_2D_vector2D(uKnot,tx,ty, x2, y2) result(uval)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real,allocatable,  dimension (:,:),intent(in)       :: uKnot
    real,allocatable,  dimension (:),intent(in)         :: tx,ty

    real, dimension (:,:),intent(in)                    :: x2
    real, dimension (:,:),intent(in)                    :: y2
    real,allocatable,  dimension (:,:)                  :: uval

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
    
    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(uKnot,1)
    ny1 = size(uKnot,2)

    nx2 = size(x2,1)
    ny2 = size(x2,2)

    if(.not. allocated(uval)) allocate(uval(nx2,ny2))

    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    iloy  = 1
    idx = 0
    idy = 0

    ! evalutae function 
    do ii = 1,nx2
        do jj = 1,ny2
          call db2val(x2(ii,jj),y2(ii,jj),idx,idy,tx,ty,nx1,ny1,kx,ky,uKnot,uval(ii,jj),iflag,inbvx,inbvy,iloy,w1,w0,.true.)
        enddo 
    enddo 
end function

function interpolate_2D_vector1D(uKnot,tx,ty, x2, y2) result( uval)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real,allocatable,  dimension (:,:),intent(in)       :: uKnot
    real,allocatable,  dimension (:),intent(in)         :: tx,ty

    real, dimension (:),intent(in)                    :: x2
    real, dimension (:),intent(in)                    :: y2
    real,allocatable,  dimension (:)                  :: uval

    ! local variables
    integer                                             :: nx1,ny1
    integer                                             :: nx2,ny2
    integer                                             :: ii,jj

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                       :: kx=3    !! order in x
    integer,parameter                       :: ky=3    !! order in y
    integer,parameter                       :: iknot=0
    real,dimension(ky)                      :: w1
    real,dimension(3*max(kx,ky))            :: w0
    integer                                 :: inbvx,inbvy,inbvz,iloy, idx,idy,iflag

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(uKnot,1)
    ny1 = size(uKnot,2)

    nx2 = size(x2)

    if(.not. allocated(uval)) allocate(uval(nx2))

    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    iloy  = 1
    idx = 0
    idy = 0

    ! evalutae function 
    do ii = 1,nx2
          call db2val(x2(ii),y2(ii),idx,idy,tx,ty,nx1,ny1,kx,ky,uKnot,uval(ii),iflag,inbvx,inbvy,iloy,w1,w0,.true.)
    enddo 
end function

function interpolate_2D_value(uKnot,tx,ty, x2, y2) result(uval)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real,allocatable,  dimension (:,:),intent(in)       :: uKnot
    real,allocatable,  dimension (:),intent(in)         :: tx,ty

    real, intent(in) :: x2, y2
    real :: uval

    ! local variables
    integer                                             :: nx1,ny1
    integer                                             :: nx2,ny2
    integer                                             :: ii,jj

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                       :: kx=3       !! order in x
    integer,parameter                       :: ky=3     !! order in y
    integer,parameter                       :: iknot=0
    real,dimension(ky)                      :: w1
    real,dimension(3*max(kx,ky))            :: w0
    integer                                 :: inbvx,inbvy,inbvz,iloy, idx,idy,iflag

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(uKnot,1)
    ny1 = size(uKnot,2)


    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    iloy  = 1
    idx = 0
    idy = 0

    call db2val(x2,y2,idx,idy,tx,ty,nx1,ny1,kx,ky,uKnot,uval,iflag,inbvx,inbvy,iloy,w1,w0,.true.)
end function
end module 
