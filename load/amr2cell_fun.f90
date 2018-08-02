subroutine a2c_count(ngridtot, nvarh, repository, xmin, xmax, ymin, &
                  & ymax, zmin, zmax, lmax, cpu_list)
  !--------------------------------------------------------------------------
  ! Ce programme calcule le cube cartesien pour les
  ! variables hydro d'une simulation RAMSES.
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none

  character(LEN=*),INTENT(IN)::repository
  character(LEN=128)::nomfich
  real(KIND=8),INTENT(IN)::xmin,xmax,ymin,ymax,zmin,zmax
  integer     ,INTENT(IN)::lmax
  integer,dimension(:),INTENT(IN)::cpu_list
  integer     ,INTENT(OUT)::ngridtot
! default value from python side is 0.

  integer::ndim,i,j,k,twotondim
  integer::ivar,ncpu,nboundary
  integer::nx,ny,nz,ilevel,idim
  integer::nlevelmax,ngrida
  integer::ind,ipos
  integer::ngridmax,icpu,ncpu_read
  real::boxlen

  integer::imin,imax,jmin,jmax,kmin,kmax
  integer,INTENT(OUT)::nvarh
  integer::ix,iy,iz
  integer::nx_full,ny_full,nz_full
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dx2
  real(KIND=8),dimension(:,:),allocatable::x,xg

  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  logical::ok,ok_cell

  type level
     integer::ilevel
     integer::ngrid
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level


  type(level),dimension(1:50)::grid

  ncpu_read = size(cpu_list)

  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)
  read(10)boxlen
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1
     kmax=int(zzmax*dble(nz_full))+1
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin
     grid(ilevel)%kmax=kmax
  end do


  !-----------------------------------------------
  ! Compute array size
  !----------------------------------------------
  ! Loop over processor files
  ngridtot = 0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=10,file=nomfich,status='old',form='unformatted')
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
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh
     read(11)
     read(11)
     read(11)
     read(11)

     ! Loop over levels
     do ilevel=1,lmax

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
           allocate(x  (1:ngrida,1:ndim))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

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
                    read(11)
                 end do
              end do
           end if
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

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
                   ngridtot = ngridtot + 1
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(ref,x,son,xg)
        endif

     end do
     ! End loop over levels

     close(10)
     close(11)
  end do

end subroutine a2c_count


subroutine a2c_load(xarr, dxarr, varr, cpuarr, refarr, repository, &
                    & xmin, xmax, ymin, ymax, zmin, zmax, lmax, &
                    & ngridtot, nvarh, cpu_list)
  implicit none

  integer,dimension(:),INTENT(IN)::cpu_list
  character(len=*),intent(IN)::repository
  character(len=128)::nomfich
  real(kind=8),intent(IN)::xmin,xmax,ymin,ymax,zmin,zmax
  integer     ,intent(IN)::lmax,ngridtot
  !default value from python side is 0.

  real(kind=8),intent(out),dimension(ngridtot,3)::xarr
  real(kind=8),intent(out),dimension(ngridtot,nvarh)::varr
  real(kind=8),intent(out),dimension(ngridtot)::dxarr
  integer,intent(out),dimension(ngridtot)::cpuarr
  logical,intent(out),dimension(ngridtot)::refarr

  integer::ndim,i,j,k,twotondim
  integer::ivar,ncpu,nboundary
  integer::nx,ny,nz,ilevel,idim
  integer::nlevelmax
  integer::ind,ipos,ngrida,icnt
  integer::ngridmax,icpu,ncpu_read
  !real::boxlen,t

  integer::imin,imax,jmin,jmax,kmin,kmax
  integer::nvarh, nvarh_org
  integer::ix,iy,iz
  integer::nx_full,ny_full,nz_full
  real(kind=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dx2
  real(kind=8),dimension(:,:),allocatable::x,xg
  real(kind=8),dimension(:,:,:),allocatable::var

  !real(kind=8),dimension(:)  ,allocatable::rho
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(kind=8),dimension(1:8,1:3)::xc
  real(kind=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(len=5)::nchar,ncharcpu
  character(len=80)::ordering
  logical::ok,ok_cell

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

  type(level),dimension(1:50)::grid

  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
!  read(10)ngrid_current
!  read(10)boxlen
  close(10)
  write(*,*) "HERE"
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1
     kmax=int(zzmax*dble(nz_full))+1
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin
     grid(ilevel)%kmax=kmax
  end do

  !------------------------------------------------
  ! Allocate arrays
  !------------------------------------------------
  icnt=0

  ! Loop over processor files
  ncpu_read = size(cpu_list)
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
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
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh_org
     read(11)
     read(11)
     read(11)
     read(11)

     ! Loop over levels
     do ilevel=1,lmax

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
           !allocate(rho(1:ngrida))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

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
                 do ivar=1,nvarh_org
                    if(j.eq.icpu .and. ivar .le. nvarh)then
                       read(11)var(:,ind,ivar)
                    else
                       read(11)
                    end if
                 end do
              end do
           end if
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

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
                    icnt = icnt + 1
                    xarr(icnt,1:3)=x(i,1:3)
                    dxarr(icnt)=dx
                    cpuarr(icnt) = icpu
                    varr(icnt,1:nvarh)=var(i,ind,1:nvarh)
                    refarr(icnt) = ref(i)
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,var,ref,x)
        endif

     end do
     ! End loop over levels

     close(10)
     close(11)

  end do
  ! End loop over cpus

end subroutine


subroutine a2c_load_level(xarr, dxarr, varr, cpuarr, refarr, repository, &
                          & xmin, xmax, ymin, ymax, zmin, zmax, lmax, &
                          & ngridtot, nvarh, cpu_list)
  implicit none

  integer,dimension(:),INTENT(IN)::cpu_list
  character(len=*),intent(in)::repository
  character(len=128)::nomfich
  real(kind=8),intent(in)::xmin,xmax,ymin,ymax,zmin,zmax
  integer     ,intent(in)::lmax,ngridtot
! default value from python side is 0.

  real(kind=8),intent(out),dimension(ngridtot,3)::xarr
  real(kind=8),intent(out),dimension(ngridtot,nvarh)::varr
  real(kind=8),intent(out),dimension(ngridtot)::dxarr
  integer,intent(out),dimension(ngridtot)::cpuarr
  logical,intent(out),dimension(ngridtot)::refarr

  integer::ndim,i,j,k,twotondim
  integer::ivar,ncpu,nboundary
  integer::nx,ny,nz,ilevel,idim
  integer::nlevelmax
  integer::ind,ipos,ngrida,icnt
  integer::ngridmax,icpu,ncpu_read
!  real::boxlen,t

  integer::imin,imax,jmin,jmax,kmin,kmax
  integer::nvarh, nvarh_org
  integer::ix,iy,iz
  integer::nx_full,ny_full,nz_full
  real(kind=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dx2
  real(kind=8),dimension(:,:),allocatable::x,xg
  real(kind=8),dimension(:,:,:),allocatable::var

  !real(kind=8),dimension(:)  ,allocatable::rho
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(kind=8),dimension(1:8,1:3)::xc
  real(kind=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(len=5)::nchar,ncharcpu
  character(len=80)::ordering
  logical::ok,ok_cell

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

  type(level),dimension(1:50)::grid

  ! Temporary space for reading labels from the info file.
!  character(LEN=128)::temp_label

!  call read_params

  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
!  read(10)ngrid_current
!  read(10)boxlen
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

!  if(ndim==2)then
!     write(*,*)'Output file contains 2D data'
!     write(*,*)'Aborting'
!     stop
!  endif

  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1
     kmax=int(zzmax*dble(nz_full))+1
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin
     grid(ilevel)%kmax=kmax
  end do

  !------------------------------------------------
  ! Allocate arrays
  !------------------------------------------------
  icnt=0

  ! Loop over processor files
  ncpu_read = size(cpu_list)
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=10,file=nomfich,status='old',form='unformatted')
!     write(*,*)'Processing file '//TRIM(nomfich)
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
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh_org
     read(11)
     read(11)
     read(11)
     read(11)

     ! Loop over levels
     do ilevel=1,lmax

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
           !allocate(rho(1:ngrida))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

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
                 do ivar=1,nvarh_org
                    if(j.eq.icpu .and. ivar .le. nvarh)then
                       read(11)var(:,ind,ivar)
                    else
                       read(11)
                    end if
                 end do
              end do
           end if
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

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
                    icnt = icnt + 1
                    xarr(icnt,1:3)=x(i,1:3)
                    dxarr(icnt)=dx
                    cpuarr(icnt) = icpu
                    varr(icnt,1:nvarh)=var(i,ind,1:nvarh)
                    refarr(icnt) = ref(i)
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,var,ref,x)
        endif

     end do
     ! End loop over levels

     close(10)
     close(11)

  end do
  ! End loop over cpus

end subroutine



!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order
!  real(kind=8),INTENT(OUT)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================
