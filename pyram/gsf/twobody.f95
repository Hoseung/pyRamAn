subroutine star_potential(md,xd,yd,zd,mg,xg,yg,zg,ms,xs,ys,zs,eps,nd,ng,ns,pots)
use omp_lib
implicit none

integer*4, intent(in) :: nd,ng,ns
real*8, intent(in), dimension(nd) :: md,xd,yd,zd
real*8, intent(in), dimension(ng) :: mg,xg,yg,zg
real*8, intent(in), dimension(ns) :: ms,xs,ys,zs
real*8, intent(in) :: eps
real*8, intent(out), dimension(ns) :: pots
integer*4 :: i,j
real*8 :: dist
integer :: num_threads

real*8, dimension(ns) :: ms2,xs2,ys2,zs2

do i=1,ns
  pots(i) = 0.0
  ms2(i) = ms(i)
  xs2(i) = xs(i)
  ys2(i) = ys(i)
  zs2(i) = zs(i)
enddo

!$omp parallel do private(i,j,dist) shared(pots,xs,ys,zs,ms,xd,yd,zd,md)
do i=1,ns
  do j=1,nd
    dist = sqrt( (xs(i)-xd(j))**2 + (ys(i)-yd(j))**2 + (zs(i)-zd(j))**2 + eps**2)
    pots(i) = pots(i)+md(j)/dist
  enddo
enddo
!$omp end parallel do

!$omp parallel do private(i,j,dist) shared(pots,xs,ys,zs,ms,xg,yg,zg,mg)
do i=1,ns
  do j=1,ng
    dist = sqrt( (xs(i)-xg(j))**2 + (ys(i)-yg(j))**2 + (zs(i)-zg(j))**2 + eps**2)
    pots(i) = pots(i)+mg(j)/dist
  enddo
enddo
!$omp end parallel do

!$omp parallel do private(i,j,dist) shared(pots,xs,ys,zs,ms,xs2,ys2,zs2,ms2)
do i=1,ns
  do j=1,ns
    dist = sqrt( (xs(i)-xs2(j))**2 + (ys(i)-ys2(j))**2 + (zs(i)-zs2(j))**2 + eps**2)
    pots(i) = pots(i)+ms2(j)/dist
  enddo
enddo
!$omp end parallel do

end subroutine star_potential


subroutine midplane_potential(md,xd,yd,zd,mg,xg,yg,zg,ms,xs,ys,zs,xin,yin,zin,eps,nd,ng,ns,ni,pote)
implicit none

integer*4, intent(in) :: nd,ng,ns,ni
real*8, intent(in), dimension(nd) :: md,xd,yd,zd
real*8, intent(in), dimension(ng) :: mg,xg,yg,zg
real*8, intent(in), dimension(ns) :: ms,xs,ys,zs
real*8, intent(in), dimension(ni) :: xin,yin,zin
real*8, intent(in) :: eps
real*8, intent(out), dimension(ni) :: pote
integer*4 :: i,j
real*8 :: dist

do i=1,ni
  pote(i) = 0.0
enddo

do i=1,ni
  do j=1,nd
    dist = sqrt( (xin(i)-xd(j))**2 + (yin(i)-yd(j))**2 + (zin(i)-zd(j))**2 + eps**2)
    pote(i) = pote(i)+md(j)/dist
  enddo
enddo

do i=1,ni
  do j=1,ng
    dist = sqrt( (xin(i)-xg(j))**2 + (yin(i)-yg(j))**2 + (zin(i)-zg(j))**2 + eps**2)
    pote(i) = pote(i)+mg(j)/dist
  enddo
enddo

do i=1,ni
  do j=1,ns
    dist = sqrt( (xin(i)-xs(j))**2 + (yin(i)-ys(j))**2 + (zin(i)-zs(j))**2 + eps**2)
    pote(i) = pote(i)+ms(j)/dist
  enddo
enddo

end subroutine midplane_potential


subroutine midplane_vcirc2(md,xd,yd,zd,mg,xg,yg,zg,ms,xs,ys,zs,xin,yin,zin,eps,nd,ng,ns,ni,vcirc2)
implicit none

integer*4, intent(in) :: nd,ng,ns,ni
real*8, intent(in), dimension(nd) :: md,xd,yd,zd
real*8, intent(in), dimension(ng) :: mg,xg,yg,zg
real*8, intent(in), dimension(ns) :: ms,xs,ys,zs
real*8, intent(in), dimension(ni) :: xin,yin,zin
real*8, intent(in) :: eps
real*8, intent(out), dimension(ni) :: vcirc2
integer*4 :: i,j
real*8 :: dist32
real*8, dimension(ni) :: acc_x, acc_y, acc_z

do i=1,ni
  acc_x(i) = 0.0
  acc_y(i) = 0.0
!   acc_z(i) = 0.0
enddo

do i=1,ni
  do j=1,nd
    dist32 = (sqrt( (xin(i)-xd(j))**2 + (yin(i)-yd(j))**2 + (zin(i)-zd(j))**2 + eps**2))**3
    acc_x(i) = acc_x(i)+md(j)/dist32*(xin(i)-xd(j))
    acc_y(i) = acc_y(i)+md(j)/dist32*(yin(i)-yd(j))
!     acc_z(i) = acc_z(i)+md(j)/dist32*(zin(i)-zd(j))
  enddo
enddo

do i=1,ni
  do j=1,ng
    dist32 = (sqrt( (xin(i)-xg(j))**2 + (yin(i)-yg(j))**2 + (zin(i)-zg(j))**2 + eps**2))**3
    acc_x(i) = acc_x(i)+mg(j)/dist32*(xin(i)-xg(j))
    acc_y(i) = acc_y(i)+mg(j)/dist32*(yin(i)-yg(j))
!     acc_z(i) = acc_z(i)+mg(j)/dist32*(zin(i)-zg(j))
  enddo
enddo

do i=1,ni
  do j=1,ns
    dist32 = (sqrt( (xin(i)-xs(j))**2 + (yin(i)-ys(j))**2 + (zin(i)-zs(j))**2 + eps**2))**3
    acc_x(i) = acc_x(i)+ms(j)/dist32*(xin(i)-xs(j))
    acc_y(i) = acc_y(i)+ms(j)/dist32*(yin(i)-ys(j))
!     acc_z(i) = acc_z(i)+ms(j)/dist32*(zin(i)-zs(j))
  enddo
enddo

do i=1,ni
  vcirc2(i) = acc_x(i)*xin(i) + acc_y(i)*yin(i)
enddo    

end subroutine midplane_vcirc2