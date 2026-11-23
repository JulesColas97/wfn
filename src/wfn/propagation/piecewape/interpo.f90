module interpo

implicit none

public :: create_2D_interpoland
public :: interpolate_2D
public :: create_1D_interpoland
public :: interpolate_1D
public :: interpolate_from_3D_to_2D
public :: interpolate_from_2D_to_1D
public :: interp2

contains

subroutine interp2(xData, yData, xVal, yVal)
  ! Inputs: xData = a vector of the x-values of the data to be interpolated
  !         yData = a vector of the y-values of the data to be interpolated
  !         xVal  = a vector of the x-values where interpolation should be performed
  ! Output: yVal  = a vector of the resulting interpolated values

  implicit none

  real, intent(in) :: xData(:), yData(:), xVal(:)
  real, intent(out) :: yVal(:)
  integer :: inputIndex, dataIndex
  real :: weight

  ! Possible checks on inputs could go here
  ! Things you may want to check:
  !   monotonically increasing xData
  !   size(xData) == size(yData)
  !   size(xVal) == size(yVal)

  do inputIndex = 1, size(xVal)
    ! print*, inputIndex, 'over',  size(xVal)
      dataIndex = 1
      do while ((xData(dataIndex+1)<xVal(inputIndex)) .or. (dataIndex+1<size(xData)))
        dataIndex = dataIndex + 1
      enddo
      weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
      yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1)
  end do
end subroutine

subroutine interpolate_from_2D_to_1D(h1, x1, y1, h2, Xs, Ys, x2, tau)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------
    real, dimension (:,:),intent(in)                    :: h1
    real, dimension (:), intent(in)                     :: x1
    real, dimension (:), intent(in)                     :: y1

    real,allocatable,  dimension (:),intent(inout)      :: h2
    real, dimension (:),intent(in)                      :: x2

    real, intent(in)                                    :: Xs
    real, intent(in)                                    :: Ys
    real, intent(in)                                    :: tau


    real,allocatable, dimension (:,:)                 :: hKnot
    real,allocatable, dimension (:,:)                 :: h1_new
    real,allocatable,dimension(:)                     :: x1_new
    real,allocatable,dimension(:)                     :: y1_new


    ! local variables
    integer                                 :: nx1,ny1, nx2
    integer                                 :: nx1_new,ny1_new
    real                                    :: dx1, dy1, dx2
    real                                    :: xmax, xmin, xmaxOverx,xmaxOvery,xminOverx,xminOvery
    integer                                 :: nxAddPos,nxAddNeg,nyAddPos,nyAddNeg
    integer                                 :: ixi, ieta,ii

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                                   :: kx    = 3     !! order in x
    integer,parameter                                   :: ky    = 3     !! order in z
    integer,parameter                                   :: iknot = 0
    real, allocatable, dimension(:)                   :: tx,ty


    real,dimension(ky)                                :: w1
    real,dimension(3*max(kx,ky))                      :: w2
    integer                                           :: iflag, inbvx,inbvy,iloy, idx,idy

    real                                              :: pi = 3.141592653589793

    real :: x_world_start, y_world_start, x_world_end, y_world_end
    real :: min_x, max_x, min_y, max_y


    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(x1)
    ny1 = size(y1)
    dx1 = x1(2) - x1(1)
    dy1 = y1(2) - y1(1)

    nx2 = size(x2)
    dx2 = x2(2) - x2(1)

    if (.not. allocated(h2)) allocate(h2(nx2))

    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    iloy  = 1
    idx = 0
    idy = 0

    !-------------------------------------------------------------------------------------------
    !       calculation of  nxi max if the 2D mesh is larger than the 3D grid in x direction
    !-------------------------------------------------------------------------------------------
    xmax = x2(nx2)
    xmin = x2(1)

    ! ! projection of 2D plane dimension to the 3D axis
    ! xmaxOverx = xmax*cos(tau)
    ! xmaxOvery = xmax*sin(tau)
    ! xminOverx = xmin*cos(tau + pi)
    ! xminOvery = xmin*sin(tau + pi)

    ! ! different cases according to the size of the projection
    ! ! add points  at the beginning or ath end of the array
    ! if (xmaxOverx >= 0) then
    !   nxAddPos = floor(max (xmaxOverx + Xs - x1(nx1) , 0.)/dx1)
    !   nxAddNeg =  floor(max(xminOverx - Xs, 0.)/dx1)
    ! else
    !   nxAddNeg =  floor(max(xmaxOverx - Xs, 0.)/dx1)
    !   nxAddPos = floor(max (xminOverx - x1(nx1) + Xs, 0.)/dx1)
    ! endif

    ! if (xmaxOvery >= 0) then
    !   nyAddPos = floor(max (xmaxOvery - y1(ny1) + Ys, 0.)/dy1) 
    !   nyAddNeg =  floor(max(xminOvery - Ys, 0.)/dy1)
    ! else
    !   nyAddNeg =  floor(max(xmaxOvery - Ys, 0.)/dy1)
    !   nyAddPos = floor(max (xminOvery - y1(ny1) + Ys, 0.)/dy1)
    ! endif
    ! nxAddNeg = 0
    x_world_start = Xs + x2(1)   * cos(tau)
    y_world_start = Ys + x2(1)   * sin(tau)
    x_world_end   = Xs + x2(nx2) * cos(tau)
    y_world_end   = Ys + x2(nx2) * sin(tau)
    
    ! 2. Find the extremes of the segment
    min_x = min(x_world_start, x_world_end)
    max_x = max(x_world_start, x_world_end)
    min_y = min(y_world_start, y_world_end)
    max_y = max(y_world_start, y_world_end)
    
    ! 3. Calculate required padding (in terms of index counts)
    ! We use ceiling to ensure we encompass the whole range
    nxAddNeg = ceiling(max(x1(1) - min_x, 0.0) / dx1)
    nxAddPos = ceiling(max(max_x - x1(nx1), 0.0) / dx1)
    nyAddNeg = ceiling(max(y1(1) - min_y, 0.0) / dy1)
    nyAddPos = ceiling(max(max_y - y1(ny1), 0.0) / dy1)

    
    print*, 'added points', nxAddNeg, nxAddPos, nyAddNeg, nyAddPos

    ! set new dimensions
    nx1_new = nx1+nxAddNeg+nxAddPos
    ny1_new = ny1+nyAddNeg+nyAddPos

    allocate(h1_new(nx1_new,ny1_new))
    allocate(hKnot(nx1_new,ny1_new))
    allocate(tx(nx1_new+kx))
    allocate(ty(ny1_new+ky))
    allocate(x1_new(nx1_new))
    allocate(y1_new(ny1_new))

    ! create the x,y,z vector
    do ii=1,nx1_new
      x1_new(ii) = dx1 * (ii-1 - nxAddNeg) + x1(1) 
    end do

    do ii=1,ny1_new
      y1_new(ii) = dy1 * (ii-1 - nyAddNeg) + y1(1)
    end do

    ! copy of the orginal values in the new array
    h1_new(nxAddneg+1:nx1_new-nxAddPos,nyAddNeg+1:ny1_new-nyAddPos) = h1(:,:)
    ! recopy of the first or last value if the domain extend further
    do ii = 1,nxAddNeg
      h1_new(ii,:) = h1_new(nxAddNeg+1,:)
    end do

    do ii = nx1_new-nxAddPos+1,nx1_new
      h1_new(ii,:) = h1_new(nx1_new-nxAddPos,:)
    end do

    do ii =1,NyAddNeg
      h1_new(:,ii) = h1_new(:,nyAddNeg+1)
    end do

    do ii = ny1_new-nyAddPos+1,ny1_new
      h1_new(:,ii) = h1_new(:,ny1_new-nyAddPos)
    end do

    ! create knots for interpolation
    call db2ink(x1_new-Xs,nx1_new,y1_new-Ys,ny1_new,h1_new, &
            & kx,ky,iknot,tx(:),ty(:),hKnot,iflag)

    do ixi = 1,nx2
       ! interpolation of v
      call db2val(cos(tau)*x2(ixi),sin(tau)*x2(ixi),idx,idy,&
                                          tx(:),ty(:),nx1_new,ny1_new,kx,ky,&
                                          hKnot,h2(ixi) ,iflag,&
                                          inbvx,inbvy,iloy,&
                                          w1(:),w2(:),.true.)
    enddo

    deallocate(h1_new)
    deallocate(hKnot)
    deallocate(tx)
    deallocate(ty)
    deallocate(x1_new)
    deallocate(y1_new)
end subroutine interpolate_from_2D_to_1D

subroutine interpolate_from_3D_to_2D(u1, x1, y1, z1, u2, Xs, Ys, x2, z2, tau)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------
    real, dimension (:,:,:),intent(in)   :: u1
    ! real, dimension (:,:,:),intent(in)   :: c1
    real, dimension (:), intent(in)      :: x1
    real, dimension (:), intent(in)      :: y1
    real, dimension (:), intent(in)      :: z1

    real,allocatable,  dimension (:,:),intent(inout)     :: u2
    ! real,allocatable,  dimension (:,:),intent(inout)     :: c2

    real, intent(in)                        :: Xs
    real, intent(in)                        :: Ys
    real, intent(in)                        :: tau

    real, dimension (:,:),intent(in)        :: x2
    real, dimension (:,:),intent(in)        :: z2

    real,allocatable, dimension (:,:,:)                 :: uKnot
    real,allocatable, dimension (:,:,:)                 :: u1_new
    real,allocatable,dimension(:)                       :: x1_new
    real,allocatable,dimension(:)                       :: y1_new
    real,allocatable,dimension(:)                       :: z1_new


    ! local variables
    integer                                 :: nx1,ny1, nz1, nx2, nz2
    integer                                 :: nx1_new,ny1_new, nz1_new
    real                                    :: dx1, dy1, dx2
    real                                    :: xmax, xmin, xmaxOverx,xmaxOvery,xminOverx,xminOvery
    integer                                 :: nxAddPos,nxAddNeg,nyAddPos,nyAddNeg
    integer                                 :: ixi, ieta,ii

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                       :: kx    = 3     !! order in x
    integer,parameter                       :: ky    = 3     !! order in y
    integer,parameter                       :: kz    = 3     !! order in z
    integer,parameter                       :: iknot = 0
    real, allocatable, dimension(:)       :: tx,ty,tz
    real,dimension(ky,kz)                 :: w1_3d
    real,dimension(kz)                    :: w2_3d
    real,dimension(3*max(kx,ky,kz))       :: w3_3d
    integer, dimension(3)                   :: inbvx,inbvy,inbvz,iloy,iloz, idx,idy,idz

    integer,dimension(6,2)                  :: iflag
    real                                    :: pi = 3.141592653589793
    real :: x_world_start, y_world_start, x_world_end, y_world_end
    real :: min_x, max_x, min_y, max_y

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(x1)
    ny1 = size(y1)
    nz1 = size(z1)
    dx1 = x1(2) - x1(1)
    dy1 = y1(2) - y1(1)

    nx2 = size(x2, 1)
    nz2 = size(x2, 2)
    dx2 = x2(2,1) - x2(1,1)

    if (.not. allocated(u2)) allocate(u2(nx2,nz2))
    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    inbvz = 1
    iloy  = 1
    iloz  = 1
    idx = 0
    idy = 0
    idz = 0

    !-------------------------------------------------------------------------------------------
    !       calculation of  nxi max if the 2D mesh is larger than the 3D grid in x direction
    !-------------------------------------------------------------------------------------------

    xmax = x2(nx2,1)
    xmin = x2(1,1)

    ! ! projection of 2D plane dimension to the 3D axis
    ! xmaxOverx = xmax*cos(tau)
    ! xmaxOvery = xmax*sin(tau)
    ! xminOverx = xmin*cos(tau + pi)
    ! xminOvery = xmin*sin(tau + pi)

    ! ! different cases according to the size of the projection
    ! ! add points  at the beginning or ath end of the array
    ! if (xmaxOverx >= 0) then
    !   nxAddPos = floor(max (xmaxOverx + Xs - x1(nx1) , 0.)/dx1)
    !   nxAddNeg =  floor(max(xminOverx - Xs, 0.)/dx1)
    ! else
    !   nxAddNeg =  floor(max(xmaxOverx - Xs, 0.)/dx1)
    !   nxAddPos = floor(max (xminOverx - x1(nx1) + Xs, 0.)/dx1)
    ! endif

    ! if (xmaxOvery >= 0) then
    !   nyAddPos = floor(max (xmaxOvery - y1(ny1) + Ys, 0.)/dy1)
    !   nyAddNeg =  floor(max(xminOvery - Ys, 0.)/dy1)
    ! else
    !   nyAddNeg =  floor(max(xmaxOvery - Ys, 0.)/dy1)
    !   nyAddPos = floor(max (xminOvery - y1(ny1) + Ys, 0.)/dy1)
    ! endif
    
    x_world_start = Xs + x2(1,1)   * cos(tau)
    y_world_start = Ys + x2(1,1)   * sin(tau)
    x_world_end   = Xs + x2(nx2,1) * cos(tau)
    y_world_end   = Ys + x2(nx2,1) * sin(tau)
    
    ! 2. Find the extremes of the segment
    min_x = min(x_world_start, x_world_end)
    max_x = max(x_world_start, x_world_end)
    min_y = min(y_world_start, y_world_end)
    max_y = max(y_world_start, y_world_end)
    
    ! 3. Calculate required padding (in terms of index counts)
    ! We use ceiling to ensure we encompass the whole range
    nxAddNeg = ceiling(max(x1(1) - min_x, 0.0) / dx1)
    nxAddPos = ceiling(max(max_x - x1(nx1), 0.0) / dx1)
    nyAddNeg = ceiling(max(y1(1) - min_y, 0.0) / dy1)
    nyAddPos = ceiling(max(max_y - y1(ny1), 0.0) / dy1)
    print*, 'added points', nxAddNeg, nxAddPos, nyAddNeg, nyAddPos

    ! set new dimensions
    nx1_new = nx1+nxAddNeg+nxAddPos
    ny1_new = ny1+nyAddNeg+nyAddPos
    nz1_new = nz1


    allocate(u1_new(nx1_new,ny1_new,nz1_new))
    allocate(uKnot(nx1_new,ny1_new,nz1_new))
    allocate(tx(nx1_new+kx))
    allocate(ty(ny1_new+ky))
    allocate(tz(nz1_new+kz))
    allocate(x1_new(nx1_new))
    allocate(y1_new(ny1_new))
    allocate(z1_new(nz1_new))


    ! create the x,y,z vector
    do ii=1,nx1_new
      x1_new(ii) = dx1 * (ii-1 - nxAddNeg) + x1(1)
    end do

    do ii=1,ny1_new
      y1_new(ii) = dy1 * (ii-1 - nyAddNeg) + y1(1)
    end do

    z1_new = z1

    ! copy of the orginal values in the new array
    u1_new(nxAddneg+1:nx1_new-nxAddPos,nyAddNeg+1:ny1_new-nyAddPos,:) = u1
    ! recopy of the first or last value if the domain extend further
    do ii = 1,nxAddNeg
      u1_new(ii,:,:) = u1_new(nxAddNeg+1,:,:)
    end do

    do ii = nx1_new-nxAddPos+1,nx1_new
      u1_new(ii,:,:) = u1_new(nx1_new-nxAddPos,:,:)
    end do

    do ii =1,NyAddNeg
      u1_new(:,ii,:) = u1_new(:,nyAddNeg+1,:)
    end do

    do ii = ny1_new-nyAddPos+1,ny1_new
      u1_new(:,ii,:) = u1_new(:,ny1_new-nyAddPos,:)
    end do

    x1_new = x1_new - Xs
    y1_new = y1_new - Ys
    ! create knots for interpolation
    call db3ink(x1_new,nx1_new,y1_new,ny1_new,z1_new,nz1_new,u1_new, &
            & kx,ky,kz,iknot,tx(:),ty(:),tz(:),uKnot,iflag(3,1))

    ! loop over the 2D points of the domain
    do ieta = 1,nz2
        do ixi = 1,nx2
            ! interpolation of u
            call db3val(cos(tau)*x2(ixi,ieta),sin(tau)*x2(ixi,ieta),z2(ixi,ieta),idx(1),idy(1),idz(1),&
                                                tx(:),ty(:),tz(:),nx1_new,ny1_new,nz1_new,kx,ky,kz,&
                                                uKnot, u2(ixi,ieta), iflag(3,1),&
                                                inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1),&
                                                w1_3d(:,:),w2_3d(:),w3_3d(:),.true.)

        enddo
    enddo
    
    deallocate(u1_new)
    deallocate(uKnot)
    deallocate(tx)
    deallocate(ty)
    deallocate(tz)
    deallocate(x1_new)
    deallocate(y1_new)
    deallocate(z1_new)
end subroutine interpolate_from_3D_to_2D

subroutine interpolate_from_2D_to_2D(u1,  x1, z1, u2,  Xs, x2, z2, tau)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------
    real, dimension (:,:),intent(in)                    :: u1
    ! real, dimension (:,:),intent(in)                    :: c1
    real, dimension (:), intent(in)                     :: x1
    real, dimension (:), intent(in)                     :: z1

    real,allocatable,  dimension (:,:),intent(inout)    :: u2
    ! real,allocatable,  dimension (:,:),intent(inout)    :: c2

    real, intent(in)                                    :: Xs
    real, intent(in)                                    :: tau

    real, dimension (:,:),intent(in)                    :: x2
    real, dimension (:,:),intent(in)                    :: z2

    real,allocatable, dimension (:,:)                   :: uKnot
    real,allocatable, dimension (:,:)                   :: cKnot
    real,allocatable, dimension (:,:)                   :: u1_new
    real,allocatable, dimension (:,:)                   :: c1_new
    real,allocatable,dimension(:)                       :: x1_new
    real,allocatable,dimension(:)                       :: z1_new


    ! local variables
    integer                                             :: nx1, nz1, nx2, nz2
    integer                                             :: nx1_new, nz1_new
    real                                                :: dx1, dx2
    real                                                :: xmax, xmin, xmaxOverx,xminOverx
    integer                                             :: nxAddPos,nxAddNeg
    integer                                             :: ixi, ieta,ii

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                                   :: kx    = 3     !! order in x
    integer,parameter                                   :: kz    = 3     !! order in z
    integer,parameter                                   :: iknot = 0
    real, allocatable, dimension(:,:)                   :: tx,tz


    real,dimension(kz,2)                                :: w1
    real,dimension(3*max(kx,kz),2)                      :: w2
    integer, dimension(3)                               :: inbvx,inbvz,iloz, idx,idz

    integer,dimension(6,2)                              :: iflag
    real                                                :: pi = 3.141592653589793

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(x1)
    nz1 = size(z1)
    dx1 = x1(2) - x1(1)

    nx2 = size(x2, 1)
    nz2 = size(x2, 2)
    dx2 = x2(2,1) - x2(1,1)

    if (.not. allocated(u2)) allocate(u2(nx2,nz2))
    ! if (.not. allocated(c2)) allocate(c2(nx2,nz2))

    !have to set these before the first evaluate call:
    inbvx = 1
    inbvz = 1
    iloz  = 1
    idx = 0
    idz = 0

    !-------------------------------------------------------------------------------------------
    !       calculation of  nxi max if the 2D mesh is larger than the 3D grid in x direction
    !-------------------------------------------------------------------------------------------

    xmax = x2(nx2,1)
    xmin = x2(1,1)

    ! projection of 2D plane dimension to the 3D axis
    xmaxOverx = xmax*cos(tau)
    xminOverx = xmin*cos(tau + pi)

    ! different cases according to the size of the projection
    ! add points  at the beginning or ath end of the array
    if (xmaxOverx >= 0) then
      nxAddPos = floor(max (xmaxOverx + Xs - x1(nx1) , 0.)/dx1)
      nxAddNeg =  floor(max(xminOverx - Xs, 0.)/dx1)
    else
      nxAddNeg =  floor(max(xmaxOverx - Xs, 0.)/dx1)
      nxAddPos = floor(max (xminOverx - x1(nx1) + Xs, 0.)/dx1)
    endif

    
    print*, 'added ', nxAddNeg, NxAddPos

    ! set new dimensions
    nx1_new = nx1+nxAddNeg+nxAddPos
    nz1_new = nz1


    allocate(u1_new(nx1_new,nz1_new))
    ! allocate(c1_new(nx1_new,nz1_new))
    allocate(uKnot(nx1_new,nz1_new))
    ! allocate(cKnot(nx1_new,nz1_new))
    allocate(tx(nx1_new+kx,2))
    allocate(tz(nz1_new+kz,2))
    allocate(x1_new(nx1_new))
    allocate(z1_new(nz1_new))


    ! create the x,y,z vector
    do ii=1,nx1_new
      x1_new(ii) = dx1 * (ii-1 - nxAddNeg) + x1(1)
    end do

    z1_new = z1

    ! copy of the orginal values in the new array
    ! and reduce dimension 
    u1_new(nxAddneg+1:nx1_new-nxAddPos,:) = u1(:,:)
    ! c1_new(nxAddneg+1:nx1_new-nxAddPos,:) = c1(:,:)
    ! recopy of the first or last value if the domain extend further
    do ii = 1,nxAddNeg
      u1_new(ii,:) = u1_new(nxAddNeg,:)
      ! c1_new(ii,:) = c1_new(nxAddNeg,:)
    end do

    do ii = nx1_new-nxAddPos+1,nx1_new
      u1_new(ii,:) = u1_new(nx1_new-nxAddPos,:)
      ! c1_new(ii,:) = c1_new(nx1_new-nxAddPos,:)
    end do

    ! create knots for interpolation
    call db2ink(x1_new-Xs,nx1_new,z1_new,nz1_new,u1_new, &
            & kx,kz,iknot,tx(:,1),tz(:,1),uKnot,iflag(3,1))
    ! call db2ink(x1_new-Xs,nx1_new,z1_new,nz1_new,c1_new,&
    !         & kx,kz,iknot,tx(:,2),tz(:,2),cKnot,iflag(3,2))

    ! loop over the 2D points of the domain
    do ieta = 1,nz2
        do ixi = 1,nx2
            ! interpolation of c
           ! call db2val(cos(tau)*x2(ixi,ieta),z2(ixi,ieta),idx(1),idz(1),&
           !                                     tx(:,1),tz(:,1),nx1_new,nz1_new,kx,kz,&
           !                                     cKnot,c2(ixi,ieta) ,iflag(1,1),&
           !                                     inbvx(1),inbvz(1),iloz(1),&
           !                                     w1(:,1),w2(:,1),.true.)
            
            ! interpolation of u
           call db2val(cos(tau)*x2(ixi,ieta),z2(ixi,ieta),idx(2),idz(2),&
                                               tx(:,2),tz(:,2),nx1_new,nz1_new,kx,kz,&
                                               uKnot,u2(ixi,ieta) ,iflag(1,2),&
                                               inbvx(2),inbvz(2),iloz(2),&
                                               w1(:,2),w2(:,2),.true.)

        enddo
    enddo
    
    deallocate(u1_new)
    ! deallocate(c1_new)
    deallocate(uKnot)
    ! deallocate(cKnot)
    deallocate(tx)
    deallocate(tz)
    deallocate(x1_new)
    deallocate(z1_new)
end subroutine interpolate_from_2D_to_2D


subroutine interpolate_from_1D_to_2D(u1,z1,u2,x2,z2,tau)
  use bspline_sub_module
  implicit none
    
  real, dimension (:),intent(in)                            :: u1
  real, dimension (:), intent(in)                           :: z1

  real,allocatable,  dimension (:,:),intent(inout)          :: u2

  real, intent(in)                                          :: tau

  real, dimension (:,:),intent(in)                          :: x2
  real, dimension (:,:),intent(in)                          :: z2

  ! local variables
  integer                                                   :: nz1,nx2,nz2
  integer                                                   :: kz = 3     !! order in x
  integer                                                   :: iknot = 0
  real, dimension(:),allocatable                            :: tz
  real,dimension(3*3)                                       :: w1_1d
  integer                                                   :: iflag
  integer                                                   :: inbvz, idz
  real, dimension(:), allocatable                           :: knot
  integer                                                   :: ixi, ieta

  idz = 0
  inbvz = 1

  nz1 = size(z1)
  nx2 = size(x2, 1)
  nz2 = size(x2, 2)
  if (.not. allocated(u2)) allocate(u2(nx2,nz2))


  nx2 = size(x2, 1)
  nz2 = size(x2, 2)
  allocate(knot(nz1))
  allocate(tz(nz1+3))
  ! real part interpolation

  call db1ink(z1,nz1,u1,kz,iknot,tz,knot,iflag)
  do ieta = 1,nz2
    do ixi = 1,nx2
        call db1val(z2(ixi, ieta),idz,tz,nz1,kz,knot,u2(ixi,ieta),iflag,inbvz,w1_1d)
    enddo 
  enddo


  ! idz = 0
  ! inbvz = 1
  !
  ! call db1ink(z1,nz1,c1,kz,iknot,tz,knot,iflag)
  !
  ! do ieta = 1,nz2
  !   do ixi = 1,nx2
  !       call db1val(z2(ixi, ieta),idz,tz,nz1,kz,knot,c2(ixi,ieta),iflag,inbvz,w1_1d)
  !   enddo 
  ! enddo

end subroutine interpolate_from_1D_to_2D


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
    integer,parameter                       :: kx    = 3     !! order in x
    integer,parameter                       :: ky    = 3     !! order in y
    integer,parameter                       :: iknot = 0
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


subroutine interpolate_2D(uKnot,tx,ty, x2, y2, uval)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real,allocatable,  dimension (:,:),intent(in)       :: uKnot
    real,allocatable,  dimension (:),intent(in)         :: tx,ty

    real, dimension (:,:),intent(in)                    :: x2
    real, dimension (:,:),intent(in)                    :: y2
    real,allocatable,  dimension (:,:),intent(inout)    :: uval

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
end subroutine interpolate_2D


subroutine create_1D_interpoland(h1, x1, hKnot,tx)
    use bspline_sub_module
    implicit none

    ! inout variables
    real, dimension(:), intent(inout)                         :: h1
    real, dimension(:), intent(inout)                         :: x1
    real,allocatable,  dimension (:),intent(inout)            :: hKnot
    real,allocatable,  dimension (:),intent(inout)            :: tx


    ! local variables
    integer                                                   :: kx = 3     !! order in x
    integer                                                   :: iknot = 0
    integer                                                   :: iflag, nx1
    integer                                                   :: inbvx, idx


    nx1 = size(x1)

    if (.not. allocated(hKnot)) allocate(hKnot(nx1))
    if (.not. allocated(tx)) allocate(tx(nx1+4))

    idx = 0
    inbvx = 1

    ! real part interpolation
    call db1ink(x1,nx1,h1,kx,iknot,tx,hKnot,iflag)

end subroutine create_1D_interpoland


subroutine interpolate_1D(hKnot,tx, x2, hval)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real, dimension (:),intent(in)         :: hKnot
    real,  dimension (:),intent(in)         :: tx

    real, dimension (:),intent(in)                      :: x2
    real,allocatable,  dimension (:),intent(inout)      :: hval

    ! local variables
    integer                                             :: nx1
    integer                                             :: nx2
    integer                                             :: ii

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                                   :: kx=3     !! order in x
    real,dimension(3*4)                                 :: w1
    integer                                             :: inbvx, idx, iflag

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(hKnot)
    nx2 = size(x2)

    if(.not. allocated(hval)) allocate(hval(nx2))

    !have to set these before the first evaluate call:
    idx = 0
    inbvx = 1

    ! evalutae function 
    do ii = 1,nx2
        call db1val(x2(ii),idx,tx,nx1,kx,hKnot,hval(ii),iflag,inbvx,w1)
    enddo 
end subroutine interpolate_1D


end module interpo
