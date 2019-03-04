subroutine potential(md,xd,yd,zd,pots)
use omp_lib
implicit none

integer*4, intent(in) :: nd
real*8, intent(in), dimension(nd) :: md,xd,yd,zd
!real*8, intent(in), dimension(ng) :: mg,xg,yg,zg
!real*8, intent(in), dimension(ns) :: ms,xs,ys,zs
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
!$omp end parallel do

end subroutine potential
