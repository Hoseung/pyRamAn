program aaa
write(*,*)"aa"
call amr2cell("/home/hoseung/Work/data/036370/snapshots/output_00035/", 0.0, 0.6, 0.0, 0.6, 0.0, 0.6, 8, 100, 100, 100)
end program aaa

subroutine amr2cell(repository, xmin, xmax, ymin, ymax, zmin, zmax, lmax, nx_sample, ny_sample, nz_sample)!, bit_length, npoint)
  implicit none

!  integer     ,INTENT(IN)                     ::bit_length,npoint
!  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
!  real(kind=8),INTENT(OUT),dimension(1:npoint)::order
  character(LEN=*),INTENT(IN)::repository!,outfich,filetype='bin'
  character(LEN=128)::nomfich,outfich="test.dat"
  real(KIND=8),INTENT(IN)::xmin,xmax,ymin,ymax,zmin,zmax
!  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  integer     ,INTENT(IN)::lmax
! default value from python side is 0.
!  integer::nx_sample=0,ny_sample=0,nz_sample=0
  integer     ,INTENT(IN)::nx_sample,ny_sample,nz_sample

  integer::ndim,n,i,j,k,twotondim,ncoarse,type=0,ii
  integer::ivar,nvar,ncpu,ncpuh,nboundary,ngrid_current
  integer::nx,ny,nz,ilevel,idim,jdim,kdim,icell
  integer::nlevelmax,ilevel1,ngrid1,ngridtot
  integer::nlevelmaxs,nlevel,iout
  integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer::ngridmax,nstep_coarse,icpu,ncpu_read
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  real::smallr,smallc
  real::boxlen,boxlen2
  real::t,aexp,hexp,t2,aexp2,hexp2
  real::omega_m,omega_l,omega_k,omega_b
  real::scale_l,scale_d,scale_t
  real::omega_m2,omega_l2,omega_k2,omega_b2

  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
!  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin,levelmin
  integer::ix,iy,iz,ndom,impi,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,dummy
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dx2
  real(KIND=8),dimension(:,:),allocatable::x,xg
  real(KIND=8),dimension(:,:,:),allocatable::var
!  real(KIND=8),dimension(:,:,:),allocatable::varr
  real(KIND=8),dimension(:,:),allocatable::xarr,varr
  real(KIND=8),dimension(:),allocatable::dxarr

!  real(KIND=4),dimension(:,:,:),allocatable::toto
  real(KIND=8),dimension(:)  ,allocatable::rho
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:)  ,allocatable::isp
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  logical::ok,ok_part,ok_cell,aa
  real(KIND=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
!  integer,dimension(:),allocatable::cpu_list
  integer,dimension(1)::cpu_list=(/1/)!,2,3/)
  character(LEN=1)::proj='z'

  type level
     integer::ilevel
     integer::ngrid
     real(KIND=4),dimension(:,:,:),pointer::cube
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level

  type(level),dimension(1:100)::grid

  ! Temporary space for reading labels from the info file.
  character(LEN=128)::temp_label

  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  write(*,*)repository
  ipos=INDEX(repository,'output_')
  write(*,*)"ipos",ipos
  nchar=repository(ipos+7:ipos+13)
  write(*,*)"nchar",nchar
  nomfich="/home/hoseung/Work/data/036370/snapshots/output_00035/hydro_00035"//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  nomfich="/home/hoseung/Work/data/036370/snapshots/output_00035/amr_00035"//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich="/home/hoseung/Work/data/036370/snapshots/output_00035/amr_00035"//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)ngrid_current
  read(10)boxlen
  close(10)
  twotondim=2**ndim
  write(*,*)"2**ndim",twotondim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------
  nomfich=TRIM(outfich)
!  write(*,*)'Writing file '//TRIM(nomfich)
!  open(unit=20,file=nomfich,form='formatted')

  ! Loop over processor files
  ncpu_read=1
  do k=1,ncpu_read
     Write(*,*)"Loop over cpu"
     icpu=cpu_list(k)
!     write(ncharcpu, "(A5,I05)") icpu
     write(ncharcpu, "(I5.5)") icpu
     ! Open AMR file and skip header
     nomfich="/home/hoseung/Work/data/036370/snapshots/output_00035/amr_00035.out00001"
     open(unit=10,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     do i=1,21
        read(10)
     end do
     ! Read grid numbers
     read(10)ngridlevel
     ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
     read(10)
     if(nboundary>0)then
        do i=1,2
           read(10)
        end do
        read(10)ngridbound
        ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
     endif
     read(10)
! ROM: comment the single follwing line for old stuff
     read(10)
     if(TRIM(ordering).eq.'bisection')then
        do i=1,5
           read(10)
        end do
     else
        read(10)
     endif
     read(10)
     read(10)
     read(10)

     ! Open HYDRO file and skip header
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
!     nomfich="/home/hoseung/Work/data/036370/snapshots/output_00035/hydro_00035.out00001"
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh
     read(11)
     read(11)
     read(11)
     read(11)


     ngridtot = 0
     do ilevel=1, lmax
        ngridtot = ngridtot + ngridfile(icpu,ilevel)
     enddo
     write(*,*)"ngridtot", ngridtot
     ngridtot=2
     allocate(varr(1:ngridtot,1:nvarh))
     allocate(xarr  (1:ngridtot,1:ndim))
     allocate(dxarr(1:ngridtot))
     ii = 0
     ! Loop over levels
     do ilevel=1,lmax
        write(*,*)"ilevel",ilevel, lmax

        ! Geometry
        dx=0.5**ilevel
        dx2=0.5*dx
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(xg(1:ngrida,1:ndim))
           allocate(son(1:ngrida,1:twotondim))
           allocate(var(1:ngrida,1:twotondim,1:nvarh))
           allocate(x  (1:ngrida,1:ndim))
           allocate(rho(1:ngrida))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu
!           write(*,*)"over domains",j
           ! Read AMR data
           if(ngridfile(j,ilevel)>0)then
              read(10) ! Skip grid index
              read(10) ! Skip next index
              read(10) ! Skip prev index
              ! Read grid center
              do idim=1,ndim
                 if(j.eq.icpu)then
                    read(10)xg(:,idim)
                 else
                    read(10)
                 endif
              end do
              read(10) ! Skip father index
              do ind=1,2*ndim
                 read(10) ! Skip nbor index
              end do
              ! Read son index
              do ind=1,twotondim
                 if(j.eq.icpu)then
                    read(10)son(:,ind)
                 else
                    read(10)
                 end if
              end do
              ! Skip cpu map
              do ind=1,twotondim
                 read(10)
              end do
              ! Skip refinement map
              do ind=1,twotondim
                 read(10)
              end do
           endif

           ! Read HYDRO data
           read(11)
           read(11)
           if(ngridfile(j,ilevel)>0)then
              ! Read hydro variables
              do ind=1,twotondim
                 do ivar=1,nvarh
                    if(j.eq.icpu)then
                       read(11)var(:,ind,ivar)
                    else
                       read(11)
                    end if
                 end do
              end do
           end if
        end do
	write(*,*)"Read done"
        ! Compute map
        write(*,*)"ngrida",ngrida
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim
!              write(*,*)"Check1"

              ! Compute cell center
              do i=1,ngrida
                 x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                 x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                 x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
              end do
              ! Check if cell is refined
              do i=1,ngrida
                 ref(i)=son(i,ind)>0.and.ilevel<lmax
              end do
              ! Store data cube
              do i=1,ngrida
                 ok_cell= .not.ref(i).and. &
                      & (x(i,1)+dx2)>=xmin.and.&
                      & (x(i,2)+dx2)>=ymin.and.&
                      & (x(i,3)+dx2)>=zmin.and.&
                      & (x(i,1)-dx2)<=xmax.and.&
                      & (x(i,2)-dx2)<=ymax.and.&
                      & (x(i,3)-dx2)<=zmax
                 if(ok_cell)then
!                    write(*,*)ii
		    xarr(ii,1:3)=x(i,1:3)
		    dxarr(ii)=x(i,1)
		    varr(ii,ivar)=var(i,ind,ivar)
!		    write(*,*)xarr(ii,1:3)
!		    write(*,*)dxarr(ii)
!		    write(*,*)varr(ii,1)
!                    write(*,*)"done"
!                    write(20,999)x(i,1),x(i,2),x(i,3),dx,icpu,ilevel,&
!                         & (var(i,ind,ivar),ivar=1,nvarh)
                    ii = ii + 1
                 end if
              end do

           end do
           ! End loop over cell
           write(*,*)"Save done; ioct, ilevel",ind, ilevel
           aa = allocated(xg)
	   write(*,*)aa
           deallocate(xg)! if ngrida, deallocate arrays.
           write(*,*)"1"
           deallocate(son)
           write(*,*)"2"
           deallocate(var)
           write(*,*)"3"
           deallocate(ref)
           write(*,*)"4"
           deallocate(rho)
           write(*,*)"5"
           deallocate(x)
           write(*,*)"6"
           aa = allocated(xg)
	   write(*,*)aa
        endif! if ngrida > 0
	write(*,*)"All done",ilevel

     end do
     ! End loop over levels
     write(*,*)"levels over"
     close(10)
     close(11)

!     deallocate(xarr,varr,dxarr)! one array per output file
  end do
  ! End loop over cpus
!  close(20)
!999 format(4(1pe12.5,1x),2(i6,1x),10(e12.5,1x))

end subroutine  
