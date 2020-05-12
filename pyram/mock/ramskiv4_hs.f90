!! ========================================================================= !!
!! === RAMSKI - RAMSES wearing SKIRT ======================================= !!
!! === This code makes the files required for SKIRT simulation. ============ !!
!! === Written by Jongwon Park (Yonsei University) ========================= !!
!! ========================================================================= !!





!! ========================================================================= !!
!! == MODULES ============================================================== !!

module ramski_mods

  implicit none

  !! define structures

  !! catalog
  type catalog
    integer::np,id,level,host,sub,nsub,nextsub
    real*8::m,mvir,r,rvir,tvir,cvel,ax,ay,az,sp
    real*8,dimension(1:3)::xx,vv
  end type catalog

  !! information
  type info_struct
    real*8::unit_l,unit_d,unit_t,unit_m,t,aexp,H0,boxlen
    real*8::omega_m,omega_l,omega_k,omega_b
  end type info_struct

  !! cell
  type cell_struct
    real*8,dimension(1:3)::xx
    real*8::rr,tt,zz    !! density, temperature and metallicity
    integer::level      !! AMR level
    logical::refinement !! refinement
    real*8::dx          !! cell size
  end type cell_struct

  !! part
  type part_struct
    real*8,dimension(1:3)::xp,vp
    real*8::mp,ap,zp
  end type part_struct

  !! morton
  type morton_struct
    real*8,dimension(:),allocatable::dd !! density
    logical,dimension(:),allocatable::rr,cc !! refinement
  end type morton_struct

  contains

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine rd_info(repository,info)

    implicit none

    integer::ipos
    character(len=128)::nomfich,repository
    character(len=5)::nchar
    logical::ok
    character(len=80)::GMGM
    real*8::t,aexp,h0,unit_l,unit_d,unit_t,boxlen
    real*8::omega_m,omega_l,omega_k,omega_b

    type(info_struct),intent(out)::info

    ipos=INDEX(repository,'output_')
    nchar=repository(ipos+7:ipos+13)
    nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
    inquire(file=nomfich, exist=ok) ! verify input file 
    if ( .not. ok ) then
       print *,TRIM(nomfich)//' not found.'
       stop
    endif

    open(unit=10,file=nomfich,form='formatted',status='old')
    read(10,'(A13,I11)') !! GMGM,ncpu
    read(10,'(A13,I11)') !! GMGM,ndim
    read(10,'(A13,I11)') !! GMGM,levelmin
    read(10,'(A13,I11)') !! GMGM,levelmax
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,'(A13,E23.15)')GMGM,boxlen
    read(10,'(A13,E23.15)')GMGM,t
    read(10,'(A13,E23.15)')GMGM,aexp
    read(10,'(A13,E23.15)')GMGM,h0
    read(10,'(A13,E23.15)')GMGM,omega_m
    read(10,'(A13,E23.15)')GMGM,omega_l
    read(10,'(A13,E23.15)')GMGM,omega_k
    read(10,'(A13,E23.15)')GMGM,omega_b
    read(10,'(A13,E23.15)')GMGM,unit_l
    read(10,'(A13,E23.15)')GMGM,unit_d
    read(10,'(A13,E23.15)')GMGM,unit_t
    close(10)

    info%t=t ; info%aexp=aexp ; info%H0=H0
    info%unit_d=unit_d ; info%unit_l=unit_l ; info%unit_t=unit_t
    info%unit_m=unit_d*unit_l*unit_l*unit_l

    info%omega_m=omega_m ; info%omega_l=omega_l
    info%omega_k=omega_k ; info%omega_b=omega_b

  end subroutine rd_info

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine amr2cell(repository,arr_range,cell,ntot)
    !! This reads hydro data from RAMSES snapshot.
    !! It is just a simple modification of amr2cell.f90.

    implicit none
    integer::ndim,n,i,j,k,twotondim,ncoarse,type=0,domax=0
    integer::ivar,nvar,ncpu,ncpuh,lmax=0,nboundary,ngrid_current
    integer::nx,ny,nz,ilevel,idim,jdim,kdim,icell
    integer::nlevelmax,ilevel1,ngrid1
    integer::nlevelmaxs,nlevel,iout
    integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
    integer::ngridmax,nstep_coarse,icpu,ncpu_read
    integer::nhx,nhy,ihx,ihy,ivar1,ivar2
    real::gamma,smallr,smallc,gammah
    real::boxlen,boxlen2
    real::t,aexp,hexp,t2,aexp2,hexp2
    real::omega_m,omega_l,omega_k,omega_b
    real::scale_l,scale_d,scale_t
    real::omega_m2,omega_l2,omega_k2,omega_b2

    integer::nx_sample=0,ny_sample=0,nz_sample=0
    integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
    integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
    integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
    integer::nx_full,ny_full,nz_full,lmin,levelmin
    integer::ix,iy,iz,ndom,impi,bit_length,maxdom
    integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
    real(KIND=8),dimension(1:8)::bounding_min,bounding_max
    real(KIND=8)::dkey,order_min,dmax,dummy
    real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
    real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dx2
    real(KIND=8),dimension(:,:),allocatable::x,xg
    real(KIND=8),dimension(:,:,:),allocatable::var
    real(KIND=4),dimension(:,:,:),allocatable::toto
    real(KIND=8),dimension(:)  ,allocatable::rho
    logical,dimension(:)  ,allocatable::ref
    integer,dimension(:)  ,allocatable::isp
    integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
    real(KIND=8),dimension(1:8,1:3)::xc
    real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
    character(LEN=5)::nchar,ncharcpu
    character(LEN=80)::ordering
    character(LEN=128)::nomfich,repository,outfich,filetype='bin'
    logical::ok,ok_part,ok_cell
    real(KIND=8),dimension(:),allocatable::bound_key
    logical,dimension(:),allocatable::cpu_read
    integer,dimension(:),allocatable::cpu_list
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

    !! input from the main program
    real*8,dimension(1:6),intent(in)::arr_range

    !! return cell data to the main program
    integer,intent(out)::ntot
    type(cell_struct),dimension(:),allocatable,intent(out)::cell
    type(cell_struct),dimension(:),allocatable::cell_fake
    integer::ncount

    !! set range
    xmin=arr_range(1) ; xmax=arr_range(2)
    ymin=arr_range(3) ; ymax=arr_range(4)
    zmin=arr_range(5) ; zmax=arr_range(6)

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
    read(10)ngrid_current
    read(10)boxlen
    close(10)
    twotondim=2**ndim
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

    allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
    allocate(ngridlevel(1:ncpu,1:nlevelmax))
    if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

    if(ndim==2)then
       write(*,*)'Output file contains 2D data'
       write(*,*)'Aborting'
       stop
    endif

    nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
    open(unit=10,file=nomfich,form='formatted',status='old')
    read(10,*)
    read(10,*)
    read(10,'(A13,I11)')temp_label,levelmin
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)

    read(10,*)
    read(10,'(A13,E23.15)')temp_label,t
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)

    read(10,'(A14,A80)')temp_label,ordering
    read(10,*)
    allocate(cpu_list(1:ncpu))
    if(TRIM(ordering).eq.'hilbert')then
       allocate(bound_key(0:ncpu))
       allocate(cpu_read(1:ncpu))
       cpu_read=.false.
       do impi=1,ncpu
          read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
       end do
    endif
    close(10)

    !-----------------------
    ! Map parameters
    !-----------------------
    if(lmax==0)then
       lmax=nlevelmax
    endif
    xxmin=xmin ; xxmax=xmax
    yymin=ymin ; yymax=ymax
    zzmin=zmin ; zzmax=zmax

    if(TRIM(ordering).eq.'hilbert')then

       dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
       do ilevel=1,lmax
          dx=0.5d0**ilevel
          if(dx.lt.dmax)exit
       end do
       lmin=ilevel
       bit_length=lmin-1
       maxdom=2**bit_length
       imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
       if(bit_length>0)then
          imin=int(xxmin*dble(maxdom))
          imax=imin+1
          jmin=int(yymin*dble(maxdom))
          jmax=jmin+1
          kmin=int(zzmin*dble(maxdom))
          kmax=kmin+1
       endif

       dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
       ndom=1
       if(bit_length>0)ndom=8
       idom(1)=imin; idom(2)=imax
       idom(3)=imin; idom(4)=imax
       idom(5)=imin; idom(6)=imax
       idom(7)=imin; idom(8)=imax
       jdom(1)=jmin; jdom(2)=jmin
       jdom(3)=jmax; jdom(4)=jmax
       jdom(5)=jmin; jdom(6)=jmin
       jdom(7)=jmax; jdom(8)=jmax
       kdom(1)=kmin; kdom(2)=kmin
       kdom(3)=kmin; kdom(4)=kmin
       kdom(5)=kmax; kdom(6)=kmax
       kdom(7)=kmax; kdom(8)=kmax
     
       do i=1,ndom
          if(bit_length>0)then
             call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
          else
             order_min=0.0d0
          endif
          bounding_min(i)=(order_min)*dkey
          bounding_max(i)=(order_min+1.0D0)*dkey
       end do
     
       cpu_min=0; cpu_max=0
       do impi=1,ncpu
          do i=1,ndom
             if (   bound_key(impi-1).le.bounding_min(i).and.&
                  & bound_key(impi  ).gt.bounding_min(i))then
                cpu_min(i)=impi
             endif
             if (   bound_key(impi-1).lt.bounding_max(i).and.&
                  & bound_key(impi  ).ge.bounding_max(i))then
                cpu_max(i)=impi
             endif
          end do
       end do
     
       ncpu_read=0
       do i=1,ndom
          do j=cpu_min(i),cpu_max(i)
             if(.not. cpu_read(j))then
                ncpu_read=ncpu_read+1
                cpu_list(ncpu_read)=j
                cpu_read(j)=.true.
             endif
          enddo
       enddo
    else
       ncpu_read=ncpu
       do j=1,ncpu
          cpu_list(j)=j
       end do
    end  if

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

    ntot=0
    ! count the total number of cells at first
    ! Loop over processor files
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

       ! Loop over levels
       do ilevel=1,lmax
         ngrida=ngridfile(icpu,ilevel)
         if(ngrida>0)then
           ! Loop over cells
           do ind=1,twotondim
             ntot=ntot+ngrida
           enddo
           ! End loop over cell
         endif
       enddo

       close(10)
    end do
    ! End loop over cpus

    !! output
    allocate(cell_fake(1:ntot)) ; ncount=0

    ! Loop over processor files
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
             allocate(var(1:ngrida,1:twotondim,1:nvarh))
             allocate(x  (1:ngrida,1:ndim))
             allocate(rho(1:ngrida))
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
                      if(j.eq.icpu)then
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
                   ok_cell= (x(i,1))>=xmin .and. &
                          & (x(i,2))>=ymin .and. &
                          & (x(i,3))>=zmin .and. &
                          & (x(i,1))<=xmax .and. &
                          & (x(i,2))<=ymax .and. & 
                          & (x(i,3))<=zmax
                   if(ok_cell)then
                     ncount=ncount+1
                     cell_fake(ncount)%xx=x(i,1:3)       !! cell center
                     cell_fake(ncount)%rr=var(i,ind,1)   !! density
                     cell_fake(ncount)%tt=var(i,ind,5)   !! pressure
                     cell_fake(ncount)%zz=var(i,ind,6)   !! metallicity
                     cell_fake(ncount)%level=ilevel      !! AMR level
                     cell_fake(ncount)%refinement=ref(i) !! refinement
                     cell_fake(ncount)%dx=1d0/2**ilevel
                   endif
                end do

             end do
             ! End loop over cell

             deallocate(xg,son,var,ref,rho,x)
          endif

       end do
       ! End loop over levels

       close(10)
       close(11)

    end do
    ! End loop over cpus

    close(21)

    allocate(cell(1:ncount))
    cell(1:ncount)=cell_fake(1:ncount) ; ntot=ncount
    deallocate(cell_fake)

  end subroutine amr2cell

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine part2cube(repository,arr_range,part,ntot)
    !! This reads particle data from RAMSES snapshot.
    !! It is just a simple modification of part2cube.f90.

    implicit none
    integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,n_frw,nstar
    integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
    integer::nx=0,ny=0,nz=0,ix,iy,iz,ixp1,iyp1,izp1,idim,jdim,kdim,ncpu_read
    real(KIND=8)::mtot,ddx,ddy,ddz,dex,dey,dez,t,time,time_tot,time_simu,weight
    real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
    real(KIND=8)::aexp,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d
    integer::imin,imax,jmin,jmax,kmin,kmax,lmin
    real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,dz,deltax
    real(KIND=4),dimension(:,:,:),allocatable::toto
    real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
    real(KIND=8),dimension(:,:,:),allocatable::cube
    real(KIND=8),dimension(:,:),allocatable::x
    real(KIND=8),dimension(:)  ,allocatable::m,age
    character(LEN=1)::proj='z'
    character(LEN=5)::nchar,ncharcpu
    character(LEN=80)::ordering
    character(LEN=80)::GMGM
    character(LEN=128)::nomfich,repository,outfich
    logical::ok,ok_part,periodic=.false.,star=.false.,ageweight=.false.
    integer::impi,ndom,bit_length,maxdom
    integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
    real(KIND=8),dimension(1:8)::bounding_min,bounding_max
    real(KIND=8)::dkey,order_min,dmax
    real(kind=8),dimension(:),allocatable::bound_key
    logical,dimension(:),allocatable::cpu_read
    integer,dimension(:),allocatable::cpu_list

    !! input from the main program
    real*8,dimension(1:6),intent(in)::arr_range

    !! return particle data to the main program
    integer::ind1,ind2
    type(part_struct),dimension(:),allocatable,intent(out)::part
    type(part_struct),dimension(:),allocatable::part_fake
    integer,intent(out)::ntot

    !! set range
    xmin=arr_range(1) ; xmax=arr_range(2)
    ymin=arr_range(3) ; ymax=arr_range(4)
    zmin=arr_range(5) ; zmax=arr_range(6)

    !-----------------------------------------------
    ! Lecture du fichier particules au format RAMSES
    !-----------------------------------------------
    ipos=INDEX(repository,'output_')
    nchar=repository(ipos+7:ipos+13)
    nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
    inquire(file=nomfich, exist=ok) ! verify input file 
    if ( .not. ok ) then
       print *,TRIM(nomfich)//' not found.'
       stop
    endif

    nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
    inquire(file=nomfich, exist=ok) ! verify input file 
    if ( .not. ok ) then
       print *,TRIM(nomfich)//' not found.'
       stop
    endif

    open(unit=10,file=nomfich,form='formatted',status='old')
    read(10,'(A13,I11)')GMGM,ncpu
    read(10,'(A13,I11)')GMGM,ndim
    read(10,'(A13,I11)')GMGM,levelmin
    read(10,'(A13,I11)')GMGM,levelmax
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*) !! GMGM, boxlen
    read(10,'(A13,E23.15)')GMGM,t
    read(10,'(A13,E23.15)')GMGM,aexp
    read(10,'(A13,E23.15)')GMGM,h0
    read(10,'(A13,E23.15)')GMGM,omega_m
    read(10,'(A13,E23.15)')GMGM,omega_l
    read(10,'(A13,E23.15)')GMGM,omega_k
    read(10,'(A13,E23.15)')GMGM,omega_b
    read(10,'(A13,E23.15)')GMGM,unit_l
    read(10,'(A13,E23.15)')GMGM,unit_d
    read(10,'(A13,E23.15)')GMGM,unit_t
    read(10,*)
    read(10,'(A14,A80)')GMGM,ordering
    read(10,*)
    allocate(cpu_list(1:ncpu))
    if(TRIM(ordering).eq.'hilbert')then
       allocate(bound_key(0:ncpu))
       allocate(cpu_read(1:ncpu))
         cpu_read=.false.
       do impi=1,ncpu
          read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
       end do
    endif
    close(10)

    !-----------------------
    ! Cosmological model
    !-----------------------
    ! Allocate look-up tables
    n_frw=1000
    allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
    allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
  
    ! Compute Friedman model look up table
    call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
         & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)

    ! Find neighboring expansion factors
    i=1
    do while(aexp_frw(i)>aexp.and.i<n_frw)
       i=i+1
    end do
    ! Interploate time
    time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
         & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))

    !-----------------------
    ! Map parameters
    !-----------------------
    if(nx==0)then
       nx=2**levelmin
    endif
    if(ny==0)then
       ny=nx
    end if
    if(nz==0)then
       nz=nx
    end if
    allocate(cube(0:nx,0:ny,0:nz))
    cube=0.0d0
    idim=1
    jdim=2
    kdim=3
    xxmin=xmin ; xxmax=xmax
    yymin=ymin ; yymax=ymax
    zzmin=zmin ; zzmax=zmax
    dx=(xxmax-xxmin)/dble(nx)
    dy=(yymax-yymin)/dble(ny)
    dz=(zzmax-zzmin)/dble(nz)

    if(TRIM(ordering).eq.'hilbert')then

       dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
       do ilevel=1,levelmax
          deltax=0.5d0**ilevel
          if(deltax.lt.dmax)exit
       end do
       lmin=ilevel
       bit_length=lmin-1
       maxdom=2**bit_length
       imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
       if(bit_length>0)then
          imin=int(xmin*dble(maxdom))
          imax=imin+1
          jmin=int(ymin*dble(maxdom))
          jmax=jmin+1
          kmin=int(zmin*dble(maxdom))
          kmax=kmin+1
       endif
     
       dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
       ndom=1
       if(bit_length>0)ndom=8
       idom(1)=imin; idom(2)=imax
       idom(3)=imin; idom(4)=imax
       idom(5)=imin; idom(6)=imax
       idom(7)=imin; idom(8)=imax
       jdom(1)=jmin; jdom(2)=jmin
       jdom(3)=jmax; jdom(4)=jmax
       jdom(5)=jmin; jdom(6)=jmin
       jdom(7)=jmax; jdom(8)=jmax
       kdom(1)=kmin; kdom(2)=kmin
       kdom(3)=kmin; kdom(4)=kmin
       kdom(5)=kmax; kdom(6)=kmax
       kdom(7)=kmax; kdom(8)=kmax
     
       do i=1,ndom
          if(bit_length>0)then
             call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
          else
             order_min=0.0d0
          endif
          bounding_min(i)=(order_min)*dkey
          bounding_max(i)=(order_min+1.0D0)*dkey
       end do
       cpu_min=0; cpu_max=0
       do impi=1,ncpu
          do i=1,ndom
             if (   bound_key(impi-1).le.bounding_min(i).and.&
                  & bound_key(impi  ).gt.bounding_min(i))then
                cpu_min(i)=impi
             endif
             if (   bound_key(impi-1).lt.bounding_max(i).and.&
                  & bound_key(impi  ).ge.bounding_max(i))then
                cpu_max(i)=impi
             endif
          end do
       end do
     
       ncpu_read=0
       do i=1,ndom
          do j=cpu_min(i),cpu_max(i)
             if(.not. cpu_read(j))then
                ncpu_read=ncpu_read+1
                cpu_list(ncpu_read)=j
                cpu_read(j)=.true.
             endif
          enddo
       enddo
    else
       ncpu_read=ncpu
       do j=1,ncpu
          cpu_list(j)=j
       end do
    end  if

    !! find the total number of particles
    npart=0
    do k=1,ncpu_read
       icpu=cpu_list(k)
       call title(icpu,ncharcpu)
       nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1)ncpu2
       read(1)ndim2
       read(1)npart2
       read(1)
       read(1)nstar
       close(1)
       npart=npart+npart2
    end do

    !! store the variable
    allocate(part_fake(1:npart))
    ind1=1 ; ind2=0

    do k=1,ncpu_read
       icpu=cpu_list(k)
       call title(icpu,ncharcpu)
       nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1)ncpu2
       read(1)ndim2
       read(1)npart2
         ind2=ind2+npart2
       read(1)
       read(1)
       read(1)
       read(1)
       read(1)
       allocate(m(1:npart2))
       if(nstar>0)allocate(age(1:npart2))
       allocate(x(1:npart2,1:ndim2))
       ! Read position
       do i=1,ndim
          read(1)m
          part_fake(ind1:ind2)%xp(i)=m
       end do
       ! Skip velocity
       do i=1,ndim
          read(1)m
          part_fake(ind1:ind2)%vp(i)=m
       end do
       ! Read mass
       read(1)m
         part_fake(ind1:ind2)%mp=m
       if(nstar>0)then
          read(1) ! Skip identity
          read(1) ! Skip level
          read(1)age
          part_fake(ind1:ind2)%ap=age

          !! metallicity
          read(1)m
          part_fake(ind1:ind2)%zp=m
       endif
       ind1=ind1+npart2

       close(1)
       deallocate(x,m)
       if(nstar>0)deallocate(age)
    end do

    ntot=0
    do i=1,ind2
      if(part_fake(i)%xp(1).ge.xmin .and. part_fake(i)%xp(1).le.xmax .and. &
       & part_fake(i)%xp(2).ge.ymin .and. part_fake(i)%xp(2).le.ymax .and. &
       & part_fake(i)%xp(3).ge.zmin .and. part_fake(i)%xp(3).le.zmax .and. &
       & part_fake(i)%ap   .ne.0)ntot=ntot+1
    enddo

    allocate(part(1:ntot))
    ntot=0
    do i=1,ind2
      if(part_fake(i)%xp(1).ge.xmin .and. part_fake(i)%xp(1).le.xmax .and. &
       & part_fake(i)%xp(2).ge.ymin .and. part_fake(i)%xp(2).le.ymax .and. &
       & part_fake(i)%xp(3).ge.zmin .and. part_fake(i)%xp(3).le.zmax .and. &
       & part_fake(i)%ap   .ne.0)then
        ntot=ntot+1
        part(ntot)=part_fake(i)
      endif
    enddo

    deallocate(part_fake)

  end subroutine part2cube

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine title(n,nchar)

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

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine hilbert3d(x,y,z,order,bit_length,npoint)

    implicit none

    integer     ,INTENT(IN)                     ::bit_length,npoint
    integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
    ! real(kind=8),INTENT(OUT),dimension(1:npoint)::order
    real(kind=8),INTENT(OUT)::order

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
       ! order(ip)=0.
       ! do i=0,3*bit_length-1
       !    b0=0 ; if(i_bit_mask(i))b0=1
       !    order(ip)=order(ip)+dble(b0)*dble(2)**i
       ! end do

       ! save Hilbert key as double precision real
       order=0.
       do i=0,3*bit_length-1
          b0=0 ; if(i_bit_mask(i))b0=1
          order=order+dble(b0)*dble(2)**i
       end do
    end do

  end subroutine hilbert3d

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine rd_catalog(path,gal,ngal)

    implicit none

    integer::i
    character(len=200)::path
    integer::eof,nline,ngal

    type(catalog),dimension(:),allocatable::gal

    !! count the number of lines at first
    nline=0
    open(unit=10,file=path,action='read',status='old')
    do
      read(10,*,iostat=eof)
      if(eof<0)goto 101
      nline=nline+1
    enddo
    101 continue

    rewind 10

    !! read the file
    allocate(gal(1:nline-1)) !! rule out the header
    read(10,*) !! skip
    do i=1,nline-1
      read(10,*)gal(i)%np,gal(i)%id,gal(i)%level,gal(i)%host,gal(i)%sub,  &
              & gal(i)%nsub,gal(i)%nextsub,gal(i)%m,gal(i)%mvir,gal(i)%r, &
              & gal(i)%rvir,gal(i)%tvir,gal(i)%cvel,gal(i)%xx(1:3),       &
              & gal(i)%vv(1:3),gal(i)%ax,gal(i)%ay,gal(i)%az,gal(i)%sp
    enddo

    ngal=nline-1
    close(10)

  end subroutine rd_catalog

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
       & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

    implicit none
    integer::ntable
    real(kind=8)::O_mat_0, O_vac_0, O_k_0
    real(kind=8)::alpha,axp_min,age_tot
    real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
    ! ######################################################!
    ! This subroutine assumes that axp = 1 at z = 0 (today) !
    ! and that t and tau = 0 at z = 0 (today).              !
    ! axp is the expansion factor, hexp the Hubble constant !
    ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
    ! time, and t the look-back time, both in unit of 1/H0. !
    ! alpha is the required accuracy and axp_min is the     !
    ! starting expansion factor of the look-up table.       !
    ! ntable is the required size of the look-up table.     !
    ! ######################################################!
    real(kind=8)::axp_tau, axp_t
    real(kind=8)::axp_tau_pre, axp_t_pre
    real(kind=8)::dadtau, dadt
    real(kind=8)::dtau,dt
    real(kind=8)::tau,t
    integer::nstep,nout,nskip

    axp_tau = 1.0D0
    axp_t = 1.0D0
    tau = 0.0D0
    t = 0.0D0
    nstep = 0
  
    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau
     
       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt
     
    end do

    age_tot=-t

    nskip=nstep/ntable
  
    axp_t = 1.d0
    t = 0.d0
    axp_tau = 1.d0
    tau = 0.d0
    nstep = 0
    nout=0
    t_out(nout)=t
    tau_out(nout)=tau
    axp_out(nout)=axp_tau
    hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt
     
       if(mod(nstep,nskip)==0)then
          nout=nout+1
          t_out(nout)=t
          tau_out(nout)=tau
          axp_out(nout)=axp_tau
          hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
       end if

    end do
    t_out(ntable)=t
    tau_out(ntable)=tau
    axp_out(ntable)=axp_tau
    hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  end subroutine friedman

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine ext_cell(cell,ntotal,cellout,rr,ll,ncount)
    !! variables
    implicit none

    integer::i
    integer::ll !! refinement level under consideration
    integer::ntotal !! total number of cells
    integer::ncount !! number of cells with level 'll'
    integer,dimension(1:ntotal)::level !! array containing refinement level
    logical,dimension(1:ntotal)::ref !! record whether a cell is refined or not

    type(cell_struct),dimension(1:ntotal),intent(in)::cell

    !! output
    real(kind=8),dimension(:,:),allocatable,intent(out)::cellout
    logical,dimension(:),allocatable,intent(out)::rr

    level=cell%level
    ref=cell%refinement

    ncount=0
    do i=1,ntotal
      if(level(i).eq.ll)ncount=ncount+1
    enddo
    allocate(cellout(1:ncount,1:5),rr(1:ncount))
    ncount=0
    do i=1,ntotal
      if(level(i).eq.ll)then
        ncount=ncount+1
        cellout(ncount,1)=cell(i)%rr
        cellout(ncount,2:4)=cell(i)%xx
        cellout(ncount,5)=cell(i)%tt
        rr(ncount)=ref(i)
      endif
    enddo

  end subroutine ext_cell

  !! ====================================================================== !!
  !! ====================================================================== !!

  subroutine decompose(xp,vp,mp,npart,theta,phi)
    !! This routine decompse the galaxy into two components: bulge and disk.
    !! Then it determines two rotation angles required in the SKIRT simulation.

    !! variables
    implicit none

    integer::i,j,k
    !! input
    real(kind=8)::xp(1:npart,1:3),vp(1:npart,1:3)
    real(kind=8)::mp(1:npart)
    integer::npart

    real(kind=8)::jj(1:npart,1:3) !! angular momentum
    real(kind=8)::totjx,totjy,totjz,normj,norjx,norjy,norjz
      !! total and normalized angular momentum
    real(kind=8)::jproj(1:npart),cosa(1:npart),eps(1:npart)
      !! projected angular momentum,
      !! angle between total J and individual js, circularity parameter

    real(kind=8)::GG=6.67408d-11 !! gravitational constant
    real(kind=8)::Msun=1.98855d30 !! mass of the sun
    real(kind=8)::rr(1:npart),jcir(1:npart)
      !! distance to individual particle, circular angular momentum
    real(kind=8)::minside(1:npart) !! mass inside a particle

    integer::nbulge,ndisk,ncount

    !! disk component
    real(kind=8),dimension(:,:),allocatable::jdisk

    !! angles to the line of sight in SKIRT
    real(kind=8)::theta,phi

    !! angular momentum in the net z-direction
    xp=xp*3.08567757144d19 !! m
    vp=vp*1d3              !! m/s
    mp=mp/1d11             !! divide mass for a while
    jj(:,1) = mp * ( xp(:,2)*vp(:,3) - xp(:,3)*vp(:,2) )
    jj(:,2) = mp * ( xp(:,3)*vp(:,1) - xp(:,1)*vp(:,3) )
    jj(:,3) = mp * ( xp(:,1)*vp(:,2) - xp(:,2)*vp(:,1) )
    !! net angular momentum
    totjx=sum(jj(:,1)) ; totjy=sum(jj(:,2)) ; totjz=sum(jj(:,3))
    !! normalization
    normj=sqrt(totjx**2+totjy**2+totjz**2)
    norjx=totjx/normj ; norjy=totjy/normj ; norjz=totjz/normj
    !! projection to the net angular momentum vector
    jproj=( jj(:,1)*norjx + jj(:,2)*norjy + jj(:,3)*norjz ) / mp

    !! cosine angle between net J and individual j vectors
    cosa=jproj/sqrt(jj(:,1)**2+jj(:,2)**2+jj(:,3)**2)

    !! circular angular momentum
    mp=mp*1d11 !! return to its original value
    rr=sqrt(xp(:,1)**2+xp(:,2)**2+xp(:,3)**2)
    minside=0

    !! sort the data in the order of distance
    call quicksort3(rr,mp,cosa,jproj,jcir,1,size(rr))

    do i=2,npart
      minside(i)=minside(i-1)+mp(i-1)
      jcir(i)=rr(i)*sqrt(GG*minside(i)*Msun/rr(i))
    enddo

    !! orbital circularity parameter
    eps=jproj/jcir

    !! decompose bulge and disk
    ndisk=0
    do i=1,npart
      !! Scannapieco et al. 2009
      if( (eps(i).gt.0.5) .and. (cosa(i).gt.0.7) )ndisk=ndisk+1
    enddo
    nbulge=npart-ndisk
    !! angular momenta of disk particles
    allocate(jdisk(1:ndisk,1:3)) ; ncount=0
    do i=1,ndisk
      if( (eps(i).gt.0.5) .and. (cosa(i).gt.0.7) )then
        ncount=ncount+1
        jdisk(ncount,1:3)=jj(i,1:3)
      endif
    enddo
    totjx=sum(jdisk(:,1)) ; totjy=sum(jdisk(:,2)) ; totjz=sum(jdisk(:,3))
    normj=sqrt(totjx**2+totjy**2+totjz**2)
    norjx=totjx/normj ; norjy=totjy/normj ; norjz=totjz/normj

    theta=acos(norjz)
    if(norjy.ge.0)phi= acos(norjx/sqrt(norjx**2+norjy**2))
    if(norjy.lt.0)phi=-acos(norjx/sqrt(norjx**2+norjy**2))

    !! radian to degree
    !theta=theta*180./3.14159265359
    !phi=phi*180./3.14159265359

    !! return to the original units
    xp=xp/3.08567757144d19

  end subroutine decompose

  !! ====================================================================== !!
  !! ====================================================================== !!

  subroutine nodecompose(xp,vp,mp,npart,theta,phi)
    !! Determines the orientation of the galaxy without decomposing its component.

    !! variables
    implicit none

    integer::i,j,k
    !! input
    real(kind=8)::xp(1:npart,1:3),vp(1:npart,1:3)
    real(kind=8)::mp(1:npart)
    integer::npart

    real(kind=8)::jj(1:npart,1:3) !! angular momentum
    real(kind=8)::totjx,totjy,totjz,normj,norjx,norjy,norjz
      !! total and normalized angular momentum
    real(kind=8)::jproj(1:npart),cosa(1:npart),eps(1:npart)
      !! projected angular momentum,
      !! angle between total J and individual js, circularity parameter

    real(kind=8)::GG=6.67408d-11 !! gravitational constant
    real(kind=8)::Msun=1.98855d30 !! mass of the sun
    real(kind=8)::rr(1:npart),jcir(1:npart)
      !! distance to individual particle, circular angular momentum
    real(kind=8)::minside(1:npart) !! mass inside a particle

    integer::nbulge,ndisk,ncount

    !! disk component
    real(kind=8),dimension(:,:),allocatable::jdisk

    !! angles to the line of sight in SKIRT
    real(kind=8)::theta,phi

    !! angular momentum in the net z-direction
    xp=xp*3.08567757144d19 !! m
    vp=vp*1d3              !! m/s
    mp=mp/1d11             !! divide mass for a while
    jj(:,1) = mp * ( xp(:,2)*vp(:,3) - xp(:,3)*vp(:,2) )
    jj(:,2) = mp * ( xp(:,3)*vp(:,1) - xp(:,1)*vp(:,3) )
    jj(:,3) = mp * ( xp(:,1)*vp(:,2) - xp(:,2)*vp(:,1) )
    !! net angular momentum
    totjx=sum(jj(:,1)) ; totjy=sum(jj(:,2)) ; totjz=sum(jj(:,3))
    !! normalization
    normj=sqrt(totjx**2+totjy**2+totjz**2)
    norjx=totjx/normj ; norjy=totjy/normj ; norjz=totjz/normj

    theta=acos(norjz)
    if(norjy.ge.0)phi= acos(norjx/sqrt(norjx**2+norjy**2))
    if(norjy.lt.0)phi=-acos(norjx/sqrt(norjx**2+norjy**2))

    !! return to the original units
    xp=xp/3.08567757144d19

  end subroutine nodecompose

  !! ====================================================================== !!
  !! ====================================================================== !!

end module ramski_mods

  !! ====================================================================== !!
  !! ====================================================================== !!


!! == MODULES ============================================================== !!
!! ========================================================================= !!



!! ========================================================================= !!
!! == MAIN PROGRAM ========================================================= !!

program ramski_v4

  !! call modules
    use ramski_mods
    use tree_mods

  !! variables
    implicit none

    !! loop indices
    integer::a,i,j,k

    !! get argument
    integer::n
    character(len=8)::opt
    character(len=500)::arg
    integer::nsnap

  !! RAMSES parameters
    real*8::Lbox=100              !! simulation box size (Mpc)
    real*8::H0=70.4               !! hubble constant
    real*8::hh                    !! dimensionless hubble parameter
    real*8::aexp                  !! scale factor

  !! SKIRT instrument paramters
    real*8::ps=0.04               !! plate scale (arcsec/pixel)
    integer::npix                 !! the number of pixels
    real*8::d_ang                 !! angular diameter distance - Mpc
      integer::n_integration=5000 !! number of points for integration
    real*8::smooth=50             !! softening length - not important (pc)
    integer::extr=1               !! extract only one out of 'extr' particles

    !! SKIRT dust parameters
    real*8::Tdust=2d4 !! upper limit of T of the gas containing dust
    real*8::fdust=0.4 !! fraction of dust (0.4 comes from Saftly+ 2015)

  !! set range
    real*8::xmin,xmax,ymin,ymax,zmin,zmax !! range of interest
    integer::nx,ny,nz             !! the number of cells along each direction
    integer::lmin,lmax            !! AMR level to be considered in SKIRT
                                  !! (may differ to those in RAMSES)
    real*8,dimension(1:3)::gx     !! galaxy center
    real*8,dimension(1:3)::cx     !! cell center
    real*8::dx                    !! physical size of a coarse cell

  !! read RAMSES data
    character(len=200)::path      !! path to the RAMSES snapshot
    type(info_struct)::info       !! simulation information
    real*8::unit_d,unit_l,unit_t,unit_m,unit_T2,t
    real*8::omega_m,omega_l,omega_k,omega_b
    real*8,dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
    integer::n_frw=1000
    real*8::time_tot
    character(len=200)::dir       !! directory of RAMSES snapshot
    type(cell_struct),dimension(:),allocatable::cell !! cell data
    integer::ncell                                   !! the number of cells
    real*8,dimension(1:6)::arr_range                 !! input to amr2cell
    type(part_struct),dimension(:),allocatable::part !! particle data
    integer::npart                                   !! the number of particles

    !! lookup table
    integer*8,dimension(0:2**20-1)::bm
    common /LUT/ bm

    !! dust only in cold gas
    real*8::temper                !! gas temperature

    !! call the subroutine 'morton'
    !! follow the syntax of previous version - just because it's bothering
      real*8,dimension(:,:),allocatable::cell1
      integer::ncell1
      integer,dimension(:),allocatable::ncell_lev !! # of cells of each level

    !! store the sorted cell data in an array of structure 'morton'
    type(morton_struct),dimension(:),allocatable::arr_mort
    real*8,dimension(:),allocatable::dd1,dd2,dd3
    logical,dimension(:),allocatable::rr1,rr2,rr3,cc1,cc2,cc3

    !! calculate particle age
    real*8::lbt !! look-back time of the current snapshot
    real*8::t_universe !! age of the universe = 1/H0

    !! write the results
    character(len=500)::outdir,write_file
    character(len=5)::ngal_char,nsnap_char
    logical::mkdir
      !! ski files
      integer::nph=50000 !! the number of photons for Montecarlo simulation
      character(len=200)::skiname
      real*8::fov
      real*8,dimension(1:2)::inputfov
      real*8,dimension(1:3)::inputang
      real*8::minw=0.1,maxw=1.1 !! range of spectrum - micron (10000 Angstrom)
      integer::npoint=100 !! number of points in the spectrum
        !! rotation angles
        real*8::theta,phi
        real*8,dimension(:,:),allocatable::xp,vp
        real*8,dimension(:),allocatable::mp
          !! to use previous version without modifications
        real*8::imgx,imgy,imgz,imgx2,imgy2,imgz2
        !! assumed redshift, distance and field of view
        real*8::z_assume=0.05,d_assume
        real*8,dimension(1:2)::fov_assume
        integer::npix_assume
        character(len=10)::z_char
      !! divide the field of view into pieces
      integer::n_piece=1
      !!
      integer::nn

    !! TREE-related variables
    !!character(len=256)::fn

    !! count_tree
    !!integer::n_halos_all,nsteps,flist_index,slist_index
    !!logical::big_run
    !!integer::n_all_fathers,n_all_sons
    !! load_tree
    !!integer,dimension(:),allocatable::fatherID,fatherIDx
    !!integer,dimension(:),allocatable::sonID
    !!real*8,dimension(:),allocatable::fatherMass
    !!integer,dimension(:,:),allocatable::i_arr
    !!real*8,dimension(:,:),allocatable::f_arr
    !!real*8,dimension(:),allocatable::aexp_arr,omega_t_arr,age_univ_arr

    !! galaxies in the last snapshot
    integer::nstart,nend,tstart
    integer::ngal !! used to read the tree
    real*8::masscut=1d10

    !! history of each galaxy
    !!real*8,dimension(:,:,:),allocatable::gal !! mass and position
    !!integer::ind,ind1,ind2,ind3,ind4
    !!real*8,dimension(:),allocatable::fMass
                        !! father mass - fine the most massive one
    !!integer,dimension(:),allocatable::find
                        !! father index
    !!integer::n_fathers
    integer::bad_nouts(102)

    !! find maximum father mass
    integer,dimension(1:1)::tmp_ind

    !! unit_l
    character(len=80)::GMGM

    !! faceon and edgeon
    integer::faceon=0,edgeon=0
 
    !! NH center ill defined from the tree.
    !! Read it from a text file
    real(KIND=8),dimension(:,:), allocatable::gxs
    integer::nlines

    !! define parameter here instead of using script
    dir='/storage1/NewHorizon/OUTPUT_DIR/' !! RAMSES snapshot
    outdir='/home/hopung/Work/NH/JP' !! output directory
    !!fn='/home/hopung/Work/NH/GalaxyMaker/gal/tree.dat' !! path to the tree
    !!fn_centers='centers_599_13.txt'

    nlines = 560
    allocate(gxs(1:nlines,1:3))
    open(12, file="centers_599_13.txt", status="old")

    ! read in values
    do i=1,nlines
        read(12,*) gxs(i,:)
    enddo
    close(12)

    bad_nouts =(/106,114,119,124,130,136,141,147,159,165,172,178,184,190,197,207,213&
               &,219,224,229,235,240,245,250,254,258,262,266,270,274,280,283,287,290&
               &,294,297,301,305,312,319,433,435,438,440,444,447,449,452,454,457,460&
               &,462,465,467,470,473,476,479,482,484,487,489,492,494,497,499,501,505&
               &,508,509,511,514,517,522,525,527,530,532,535,538,540,544,547,549,552&
               &,554,557,560,562,565,568,570,573,576,578,581,584,587,589,592,595,598/)
    lmin=12     !! minimum AMR level (must be smaller than lmax)
    lmax=21     !! maximum AMR level (must be equal to lmax in RAMSES)
    Tdust=10000 !! dust survival temperature
    fdust=0.4   !! dust fraction - 0.4 from Saftly et al. 2015
    smooth=50   !! stellar particle smoothing length
    !extr=5      !! extract one particle out of a group of partices
    extr=1
    nph=50000   !! number of photons to be used in SKIRT
      nph=2d5
    ps=0.04     !! plate scale - arcsec/pixel
    minw=0.1    !! SED minimum wavelength (micron)
    maxw=1.4    !! SED maximum wavelength (micron)
    npoint=120  !! number of points in SED
    fov=60      !! field of view (kpc)
                !! This will be slightly adjusted according to place scale.
      inputfov(1)=fov ; inputfov(2)=fov
    n_piece=5   !! divide the field of view into n_piece x n_piece equivalent
                !! regions (if n_piece=4, FoV splits into 4x4=16 aread)
                !! (Use this for high quality image.)
n_piece=1
    lbox=100    !! simulation box size (Mpc)
    h0=70.4     !! Hubble constant
    masscut=3d10
    nstart=110  !! first snapshot number in the tree
    nend=599    !! last snapshot to be considered
    z_assume=1.00
    faceon=1    !! if you want this output, set any non-zero integer
    edgeon=1    !! if you want this output, set any non-zero integer

    !! READ MORTON LOOK-UP TABLE
    open(unit=10,file='LUT.txt',form='formatted',status='old')
    do i=0,2**20-1
      read(10,*)bm(i)
    enddo
    close(10)
    write(*,*)"loading LUT done"
    !! simulation time - for paticle ages
      !! path to the last RAMSES snapshot - just to read a info file
      write(nsnap_char,'(I0.5)')nend
      path=trim(dir)//'output_'//nsnap_char
      !! read info
      call rd_info(path,info)
      omega_m=info%omega_m ; omega_l=info%omega_l ; omega_k=info%omega_k
        !! keep it simple
    write(*,*)"loading info done"
    !! Allocate look-up tables
    allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
    allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
    !! Compute Friedman model look up table
    call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
         & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)

    !! READ the TREE
    !! count_tree
    !!call count_tree(fn,n_halos_all,flist_index,slist_index,nsteps,big_run,&
    !!               &n_all_fathers,n_all_sons)
    !! load_tree
    !!call load_tree(fn,fatherID,fatherIDx,sonID,fatherMass,&
    !!              &i_arr,f_arr,aexp_arr,omega_t_arr,age_univ_arr,&
    !!              &n_halos_all,n_all_fathers,n_all_sons,big_run,nsteps)
    !!  f_arr(:,1)=f_arr(:,1)*1d11 !! mass
    !!write(*,*)"loading tree done"
    !! galaxies in the last snapshot
    !!ngal=0
    !!do i=1,n_halos_all
    !!  if( i_arr(i,11).eq.(nend-tstart+1) .and. f_arr(i,1).ge.masscut)then
    !!    ngal=ngal+1
    !!    if(ngal.eq.1)ind1=i
    !!    ind2=i
    !!  endif
    !!enddo
    !!allocate(gal(ngal,nstart:nend,1:4)) ; gal=0

    !!ngal=0
    !!do i=ind1,ind2 !! do not need to run over the all haloes in the tree
    !!  if( i_arr(i,11).eq.nend-tstart+1 .and. f_arr(i,1).ge.masscut )then
    !!    ngal=ngal+1
    !!    ind=i
    !!    do j=nend,nstart,-1
    !!      gal(ngal,j,1)=f_arr(ind,1)     !! mass
    !!      gal(ngal,j,2:4)=f_arr(ind,3:5) !! position - physical Mpc
    !!
    !!      n_fathers=i_arr(ind,10)      !! number of fathers
    !!      ind3=i_arr(ind,13)           !! flist_index
    !!      ind4=ind3+n_fathers-1        !! note the subtraction
    !!
    !!      allocate(fMass(1:n_fathers),find(1:n_fathers))
    !!      fMass(1:n_fathers)=fatherMass(ind3:ind4)               
    !!      find(1:n_fathers)=fatherIDx(ind3:ind4)
    !!
    !!      tmp_ind=maxloc(fMass)
    !!      ind=find(tmp_ind(1))-1
    !!
    !!      deallocate(fMass,find)
    !!    enddo
    !!  endif
    !!enddo

    !! MAIN LOOP
    !! loop over snapshots
    do a=nend,nstart,-1
      if ( ANY( bad_nouts==a ) .or. (a .gt. 160)) then
      CYCLE
      endif
      !! store the snapshot number as a string
      write(nsnap_char,'(I0.5)')a
      gx(1:3)=gxs(nend-a+1,:)

      !! path to the RAMSES snapshot
      path=trim(dir)//'output_'//nsnap_char
      write(*,*) "main loop: path to snapshot", path
      !! read info
      call rd_info(path,info)
        aexp=info%aexp ; H0=info%H0 ; hh=H0/100d0
        unit_d=info%unit_d ; unit_l=info%unit_l ; unit_t=info%unit_t
        unit_m=info%unit_m
        unit_T2=(unit_l/unit_t)**2*1.67d-24/1.38062d-16
        omega_m=info%omega_m ; omega_l=info%omega_l ; omega_k=info%omega_k
        t_universe=1./H0 * 3.086d19 / 3.154d7 !! yr

      !! calculate angular diameter distance - assumed redshift
      call cosmo(z_assume,omega_m,omega_L,H0,n_integration,d_assume)
      npix_assume=inputfov(1)*1d3/(d_assume*1d6)/3.14159265359*180.*3600./ps+1
      do while(mod(npix_assume,n_piece).ne.0)
        npix_assume=npix_assume+1
      enddo
        !! above - the number of pixels
        !! calculate the field of view with the number of pixels
        fov_assume(1)=real(npix_assume)*ps/3600./180.*&
                     & 3.14159265359*d_assume*1d3
        fov_assume(2)=real(npix_assume)*ps/3600./180.*&
                     & 3.14159265359*d_assume*1d3

      !! look-back time of the current snapshot - for particle ages
      do i=1,n_frw-1
        if(info%t.le.tau_frw(i) .and. info%t.gt.tau_frw(i+1))lbt=  &
         & ( t_frw(i+1)-t_frw(i) ) / ( tau_frw(i+1)-tau_frw(i) ) * &
         & ( info%t - tau_frw(i) ) + t_frw(i)
      enddo

      !! loop over galaxies
!!      do i=1,ngal !!- JP
      !!do i=13,13
        !! store the galaxy number as a string
        write(ngal_char,'(I0.5)')i
        !gx(1:3)=gal(i,a,2:4)*3.086d24/unit_l+0.5
        write(*,*)gx
        !! coordinate of a galaxy on the grid
        xmin=int(gx(1)*2**lmin)/real(2**lmin,8) ; xmax=xmin+1./2**lmin
        ymin=int(gx(2)*2**lmin)/real(2**lmin,8) ; ymax=ymin+1./2**lmin
        zmin=int(gx(3)*2**lmin)/real(2**lmin,8) ; zmax=zmin+1./2**lmin
          cx=(/ (xmin+xmax)/2d0 , (ymin+ymax)/2d0 , (zmin+zmax)/2d0 /)
          !! the nubmer of cells along each direction
          nx=1 ; ny=1 ; nz=1
        !! guarantee that 3x3x3 cells are extracted
        xmin=xmin-1d0/2**lmin ; xmax=xmax+1d0/2**lmin
        ymin=ymin-1d0/2**lmin ; ymax=ymax+1d0/2**lmin
        zmin=zmin-1d0/2**lmin ; zmax=zmax+1d0/2**lmin
          !! the nubmer of cells along each direction
          nx=nx+2 ; ny=ny+2 ; nz=nz+2 !! (+1 on each side)

        dx=(xmax-xmin)*unit_l/3.086d21 !! physical size of the region
        do while (dx.le.inputfov(1))
          xmin=xmin-1d0/2**lmin ; xmax=xmax+1d0/2**lmin
          ymin=ymin-1d0/2**lmin ; ymax=ymax+1d0/2**lmin
          zmin=zmin-1d0/2**lmin ; zmax=zmax+1d0/2**lmin
          dx=(xmax-xmin)*unit_l/3.086d21
          nx=nx+2 ; ny=ny+2 ; nz=nz+2 !! (+1 on each side)
        enddo

        !! set range
        arr_range=(/ xmin , xmax , ymin , ymax , zmin , zmax /)

        !! read cells
        call amr2cell(path,arr_range,cell,ncell)

        !! read particles
        call part2cube(path,arr_range,part,npart)
          !! mass - Msun

        !! convert gas mass into dust mass and destroy dust in hot gas
        do j=1,ncell
          temper=cell(j)%tt/cell(j)%rr*unit_T2 !! temperature of gas

          !! density in Msun/pc3
          cell(j)%rr=cell(j)%rr*unit_d !! g/cm3
          cell(j)%rr=cell(j)%rr/1.989d33*3.086d18*3.086d18*3.086d18 !! Msun/pc3
          cell(j)%rr=cell(j)%rr*cell(j)%zz*fdust !! gas mass x Z x fdust

          if(temper .gt. Tdust)cell(j)%rr=0 !! dust is destroyed in hot gas
        enddo

        !! sort the cells in Morton order
        allocate(arr_mort(lmin:lmax),ncell_lev(lmin:lmax))
        do j=lmin,lmax
          call ext_cell(cell,ncell,cell1,rr1,j,ncell1)
          !! cell1 is ALL cells at level j
          call morton(cell1,rr1,ncell1,lmin,j,nx,ny,nz)
          ncell_lev(j)=ncell1
          allocate(arr_mort(j)%dd(1:ncell1),arr_mort(j)%rr(1:ncell1),&
                 & arr_mort(j)%cc(1:ncell1))
          arr_mort(j)%dd=cell1(:,1)
          arr_mort(j)%rr=rr1(:)
          arr_mort(j)%cc=rr1(:)
          deallocate(cell1,rr1)
        enddo

        !! link the cells of all levels
        j=lmin
        allocate(dd1(1:ncell_lev(j)),rr1(1:ncell_lev(j)),cc1(1:ncell_lev(j)))
        dd1=arr_mort(j)%dd
        rr1=arr_mort(j)%rr
        cc1=arr_mort(j)%cc

        do j=lmin+1,lmax
          allocate(dd2(1:ncell_lev(j)),rr2(1:ncell_lev(j)),cc2(1:ncell_lev(j)))
          dd2=arr_mort(j)%dd
          rr2=arr_mort(j)%rr
          cc2=arr_mort(j)%cc
          
          allocate(dd3(1:sum(ncell_lev(lmin:j))), &
                   rr3(1:sum(ncell_lev(lmin:j))), &
                   cc3(1:sum(ncell_lev(lmin:j)))  )

          !! link !!
          call mlink(dd1,rr1,cc1,dd2,rr2,cc2,dd3,rr3,cc3, &
             &       sum(ncell_lev(lmin:j-1)),ncell_lev(j))

          !! save the variables for the next loop
          if(j.ne.lmax)then
            deallocate(dd1,rr1,cc1)
            allocate(dd1(1:size(dd3)),rr1(1:size(dd3)),cc1(1:size(dd3)))
            dd1=dd3 ; rr1=rr3 ; cc1=cc3
            deallocate(dd2,rr2,cc2,dd3,rr3,cc3)
          endif
        enddo

       !! calculate particle age
        do j=1,npart
          do k=1,n_frw-1
            if(part(j)%ap.le.tau_frw(k) .and. part(j)%ap.gt.tau_frw(k+1))then
              part(j)%ap =                                                &
              & ( t_frw(k+1) - t_frw(k) ) / ( tau_frw(k+1)-tau_frw(k) ) * &
              & ( part(j)%ap-tau_frw(k) ) + t_frw(k)
              part(j)%ap=(lbt-part(j)%ap)*t_universe !! age in yr
              goto 104
            endif
          enddo
          104 continue
        enddo

       !! write the files
        !! mkdir
        write_file=trim(outdir)
        inquire(file=write_file,exist=mkdir)
        if(.not.mkdir)call system('mkdir '//trim(write_file))

        write_file=trim(outdir)//'/'//ngal_char
        inquire(file=write_file,exist=mkdir)
        if(.not.mkdir)call system('mkdir '//trim(write_file))

        write(ngal_char,'(I0.5)')i
        write_file=trim(outdir)//'/'//ngal_char//'/'//nsnap_char
        inquire(file=write_file,exist=mkdir)
        if(.not.mkdir)call system('mkdir '//trim(write_file))

        write_file=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//'/xy'
        inquire(file=write_file,exist=mkdir)
        if(.not.mkdir)call system('mkdir '//trim(write_file))

        if(faceon.ne.0)then
        write_file=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//'/faceon'
          inquire(file=write_file,exist=mkdir)
          if(.not.mkdir)call system('mkdir '//trim(write_file))
        endif

        if(edgeon.ne.0)then
        write_file=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//'/edgeon'
          inquire(file=write_file,exist=mkdir)
          if(.not.mkdir)call system('mkdir '//trim(write_file))
        endif

       !! cell
        write_file=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//&
                  &'/dust_ramski.txt'
        open(unit=10,file=write_file,action='write',status='replace')
          write(10,'(A2,3I4)')' ! ',nx,ny,nz
          do j=1,size(rr3)
            if(cc3(j).eqv..false.)write(10,*)dd3(j)
            if(cc3(j).eqv..true. )write(10,*)'! 2 2 2'
          enddo
        close(10)

      !! particle
        write_file=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//&
                  &'/part_ramski.txt'  
        open(unit=10,file=write_file,action='write',status='replace')
          do j=1,npart,extr
            write(10,*)(part(j)%xp(:)-gx(:))*unit_l/3.086d18,smooth,&
                     & part(j)%mp*real(extr)*unit_m/1.989d33,&
                       part(j)%zp,part(j)%ap
          enddo
        close(10)

      !! ski files
          inputfov=inputfov*1d3 !! pc
          fov_assume=fov_assume*1d3

          !! non-rotated
          inputang=(/ 0d0 , 0d0 , 90d0 /)
          arr_range(1:2)=(arr_range(1:2)-gx(1))*unit_l/3.086d18 !! pc
          arr_range(3:4)=(arr_range(3:4)-gx(2))*unit_l/3.086d18 !! pc
          arr_range(5:6)=(arr_range(5:6)-gx(3))*unit_l/3.086d18 !! pc
            !! assumed redshift
            do j=0,n_piece**2-1
              write(z_char,'(F4.2)')z_assume
              write(skiname,'(I0.5)')j
              skiname=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//&
                     &'/xy/xy_redshift_'//trim(z_char)//'_'//&
                     &trim(skiname)//'.ski'
              call writeski(skiname,nph,inputang,fov_assume,0d0,&
                          & 0d0,npix_assume,d_assume,minw,maxw,npoint,&
                          & arr_range,j,n_piece)
            enddo

          !! faceon and edgeon
          if(faceon.ne.0 .or. edgeon.ne.0)then

            allocate(xp(1:npart,1:3),vp(1:npart,1:3),mp(1:npart))

            do j=1,npart
              xp(j,1:3)=(part(j)%xp(1:3)-gx(1:3))*unit_l/3.086d21
              vp(j,1:3)=part(j)%vp(1:3)*unit_l/unit_t/1d5
              mp(j)    =part(j)%mp*unit_m/1.989d33
            enddo
            call nodecompose(xp,vp,mp,npart,theta,phi) !! find rotation angles

            imgx=gx(1)*unit_l/3.086d18
            imgy=gx(2)*unit_l/3.086d18
            imgz=gx(3)*unit_l/3.086d18
              !! rotate
              !! about y-axis at first
                imgx2= cos(theta)*imgx+sin(theta)*imgz
                imgy2=imgy
                imgz2=-sin(theta)*imgx+cos(theta)*imgz
              !! then, about z-axis
                imgx=cos(phi)*imgx2-sin(phi)*imgy2
                imgy=sin(phi)*imgx2+cos(phi)*imgy2
                imgz=imgz2
            
              imgx=0 ; imgy=0 ; imgz=0

              !! radian to degree
                theta=theta*180./3.14159265359
                phi=phi*180./3.14159265359

              !! faceon
              if(faceon.ne.0)then
                inputang=(/ theta , phi , 90d0 /)
                  !! assumed redshift
                  do j=0,n_piece**2-1
                    write(skiname,'(I0.5)')j    
                    skiname=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//&
                           &'/faceon/faceon_redshift_'//trim(z_char)//'_'//&
                           &trim(skiname)//'.ski'
                    call writeski(skiname,nph,inputang,fov_assume,imgx,imgy,&
                               & npix_assume,d_assume,minw,maxw,npoint,&
                               & arr_range,j,n_piece)
                  enddo
              endif

              !! edgeon
              if(edgeon.ne.0)then
                if(theta.ge.90.0)inputang=(/ theta-90d0 , phi , 90d0 /)
                if(theta.lt.90.0)inputang=(/ theta+90d0 , phi , 90d0 /)
                  !! assumed redshift
                  do j=0,n_piece**2-1
                    write(skiname,'(I0.5)')j
                    skiname=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//&
                           &'/edgeon/edgeon_redshift_'//trim(z_char)//'_'//&
                           &trim(skiname)//'.ski'
                    call writeski(skiname,nph,inputang,fov_assume,imgx,imgy,&
                                & npix_assume,d_assume,minw,maxw,npoint,&
                                & arr_range,j,n_piece)
                  enddo
              endif

            deallocate(xp,vp,mp)
          endif

        deallocate(cell,part,arr_mort,ncell_lev)
        deallocate(dd1,dd2,dd3,rr1,rr2,rr3,cc1,cc2,cc3)
        inputfov=inputfov/1d3 !! pc
        fov_assume=fov_assume/1d3

      !!enddo !! loop over galaxies
    enddo !! loop over snapshots

end program ramski_v4

!! == MAIN PROGRAM ========================================================= !!
!! ========================================================================= !!





!! ========================================================================= !!
!! == SUBROUTINES ========================================================== !!

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine morton(cell,ref,ncount,lmin,lmax,nx,ny,nz)
    !! Convert the real-valued coordinated to integers and  sort the data - !!
    !! in Morton order. --------------------------------------------------- !!

    !! variables
      implicit none

      !! lookup table
      integer*8,dimension(0:2**20-1)::bm
      common /LUT/ bm

      !! bit manipulation
      integer::i,j
      integer::ncount,lmin,lmax
      integer::nx,ny,nz,mindx,mindy,mindz
      real(kind=8),dimension(1:ncount,1:5)::cell
      real(kind=8),dimension(1:ncount)::xx,yy,zz,dd,pp
      integer*8,dimension(1:ncount)::intx,inty,intz !! converted from x, y and z
      integer*8,dimension(1:ncount)::indx,indy,indz !! indices x, y and z
      integer*8,dimension(1:ncount)::indc         !! indices of coarse cells
      integer*8,dimension(1:ncount)::bits
      integer::mult
      integer::nelem=1
      logical,dimension(1:ncount)::ref !! cell refinement

      integer::ncw,ind1,ind2
      integer,dimension(:),allocatable::countref !! count # of refined cells
        !! extract only cells within each coarse cells
        integer*8,dimension(:),allocatable::bitsref
        real(kind=8),dimension(:),allocatable::densref !! density
        real(kind=8),dimension(:),allocatable::presref !! pressure
        logical,dimension(:),allocatable::openref !! open or close
        real(kind=8),dimension(:),allocatable::xref,yref,zref !! cell position

    !! 
      real(kind=8)::minx,miny,minz

    !! cell center
      xx(:)=cell(:,2) ; yy(:)=cell(:,3) ; zz(:)=cell(:,4)
      dd(:)=cell(:,1) ; pp(:)=cell(:,5)

    !! convert real to integer
      mult=2**(lmax+1)

      intx=nint(xx*mult)
      inty=nint(yy*mult)
      intz=nint(zz*mult)
        !! For example, if max=4 then mult=2**5=32
        !! 0.28125, 0.34375, 0.40625 and 0.46875 become 9, 11, 13 and 15.

    !! Shift the integers.
      do i=1,ncount
        !! 9, 11, 13 and 15 become 1, 1, 1 and 1.
        indx(i)=intx(i)/2**(lmax+1-lmin)
        indy(i)=inty(i)/2**(lmax+1-lmin)
        indz(i)=intz(i)/2**(lmax+1-lmin)
        !! 9, 11, 13 and 15 become 1, 3, 5 and 7.
        !!                      (9-8, 11-8, 13-8 and 15-8)
        intx(i)=intx(i)-indx(i)*2**(lmax+1-lmin)
        inty(i)=inty(i)-indy(i)*2**(lmax+1-lmin)
        intz(i)=intz(i)-indz(i)*2**(lmax+1-lmin)
        !! bit interleaving - using lookup table
        bits(i)=bm(intx(i))+ishft(bm(inty(i)),1)+ishft(bm(intz(i)),2)
          !! x, y and z - there are no overlapping bits
        !! index of coarse cell
      enddo

    !! find coarse cells where refined ones reside
      do i=1,ncount
        indc(i)=ny*nx*indz(i)+nx*indy(i)+indx(i)
      enddo
      !! Note that 'indc' has nothing to do with Morton order.
      !! It follows Raster order.

    !! Sort the data in Raster order at first.
      call quicksort(indc,bits,xx,yy,zz,dd,pp,ref,1,size(indc))

    !! count the number of coarse cells that contains refined ones
      do i=1,ncount
        if(i.ne.1 .and. indc(i-1).ne.indc(i))nelem=nelem+1
      enddo

      allocate(countref(1:nelem)) ; countref=0
      ncw=1
      do i=1,ncount
        if(i.ne.1 .and. indc(i-1).ne.indc(i))ncw=ncw+1
        countref(ncw)=countref(ncw)+1
      enddo

      ind1=1
      ind2=countref(1)

      do i=1,nelem
        allocate(bitsref(1:countref(i)),densref(1:countref(i)))
        allocate(openref(1:countref(i)),presref(1:countref(i)))
        allocate(xref(1:countref(i)),yref(1:countref(i)),zref(1:countref(i)))
        bitsref(:)=bits(ind1:ind2)
        xref(:)=xx(ind1:ind2) ; yref(:)=yy(ind1:ind2) ; zref(:)=zz(ind1:ind2)
        openref(:)=ref(ind1:ind2)
        densref(:)=dd(ind1:ind2) ; presref(:)=pp(ind1:ind2)
        call quicksort2(bitsref,xref,yref,zref,densref,presref, &
           &            openref,1,size(bitsref))
          !! do not need to record bitsref

        xx(ind1:ind2)=xref(:) ; yy(ind1:ind2)=yref(:) ; zz(ind1:ind2)=zref(:)
        ref(ind1:ind2)=openref(:)
        dd(ind1:ind2)=densref(:)
        pp(ind1:ind2)=presref(:)

        if(i.ne.nelem)then
          ind1=ind1+countref(i)
          ind2=ind2+countref(i+1)
        endif
        deallocate(bitsref,densref,openref,xref,yref,zref,presref)
      enddo

      !! return
      cell(:,2)=xx(:) ; cell(:,3)=yy(:) ; cell(:,4)=zz(:)
      cell(:,1)=dd(:) ; cell(:,5)=pp(:)

  end subroutine morton

  !! ===================================================================== !!
  !! ===================================================================== !!

  recursive subroutine quicksort(a,b,c,d,e,f,g,h,first,last)

    implicit none
    integer*8::a(*)
    integer*8::b(*)
    real(kind=8)::c(*),d(*),e(*),f(*),g(*)
    logical::h(*)
    integer::x,t
    real(kind=8)::tf
    logical::tl
    integer*8::tb
    integer::first,last
    integer::i,j

    x=a((first+last)/2)
    i=first
    j=last
    do
      do while(a(i)<x)
        i=i+1
      enddo
      do while(x<a(j))
        j=j-1
      enddo
      if(i>=j)exit
      t=a(i) ; a(i)=a(j) ; a(j)=t
      tb=b(i) ; b(i)=b(j) ; b(j)=tb
      tf=c(i) ; c(i)=c(j) ; c(j)=tf
      tf=d(i) ; d(i)=d(j) ; d(j)=tf
      tf=e(i) ; e(i)=e(j) ; e(j)=tf
      tf=f(i) ; f(i)=f(j) ; f(j)=tf
      tf=g(i) ; g(i)=g(j) ; g(j)=tf
      tl=h(i) ; h(i)=h(j) ; h(j)=tl
      i=i+1
      j=j-1
    enddo
    if(first<i-1)call quicksort(a,b,c,d,e,f,g,h,first,i-1)
    if(j+1<last)call quicksort(a,b,c,d,e,f,g,h,j+1,last)

  end subroutine quicksort

  !! ===================================================================== !!
  !! ===================================================================== !!

  recursive subroutine quicksort2(a,b,c,d,e,f,g,first,last)

    implicit none
    integer*8::a(*)
    real(kind=8)::b(*),c(*),d(*),e(*),f(*)
    logical::g(*)
    integer*8::x,t
    real(kind=8)::tf
    logical::tl
    integer::first,last
    integer::i,j

    x=a((first+last)/2)
    i=first
    j=last
    do
      do while(a(i)<x)
        i=i+1
      enddo
      do while(x<a(j))
        j=j-1
      enddo
      if(i>=j)exit
      t=a(i) ; a(i)=a(j) ; a(j)=t
      tf=b(i) ; b(i)=b(j) ; b(j)=tf
      tf=c(i) ; c(i)=c(j) ; c(j)=tf
      tf=d(i) ; d(i)=d(j) ; d(j)=tf
      tf=e(i) ; e(i)=e(j) ; e(j)=tf
      tf=f(i) ; f(i)=f(j) ; f(j)=tf
      tl=g(i) ; g(i)=g(j) ; g(j)=tl
      i=i+1
      j=j-1
    enddo
    if(first<i-1)call quicksort2(a,b,c,d,e,f,g,first,i-1)
    if(j+1<last)call quicksort2(a,b,c,d,e,f,g,j+1,last)

  end subroutine quicksort2

  !! ===================================================================== !!
  !! ===================================================================== !!

  recursive subroutine quicksort3(a,b,c,d,e,first,last)

    implicit none
    real(kind=8)::a(*),b(*),c(*),d(*),e(*)
    real(kind=8)::x
    real(kind=8)::tf
    integer::first,last
    integer::i,j

    x=a((first+last)/2)
    i=first
    j=last
    do
      do while(a(i)<x)
        i=i+1
      enddo
      do while(x<a(j))
        j=j-1
      enddo
      if(i>=j)exit
      tf=a(i) ; a(i)=a(j) ; a(j)=tf
      tf=b(i) ; b(i)=b(j) ; b(j)=tf
      tf=c(i) ; c(i)=c(j) ; c(j)=tf
      tf=d(i) ; d(i)=d(j) ; d(j)=tf
      tf=e(i) ; e(i)=e(j) ; e(j)=tf
      i=i+1
      j=j-1
    enddo
    if(first<i-1)call quicksort3(a,b,c,d,e,first,i-1)
    if(j+1<last)call quicksort3(a,b,c,d,e,j+1,last)

  end subroutine quicksort3

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine mlink(tr1,to1,tc1,tr2,to2,tc2,tr3,to3,tc3,levgrid,levgrid2)

    !! -------------------------------------------------------------------- !!
    !!     1     2     3     4     5     6     7    ...     level   i       !! 
    !!           ^                       ^                                  !!
    !!       1  2  3  4              9 10 11 12     ...     level i+1       !!
    !!       5  6  7  8             13 14 15 16                             !!
    !! -------------------------------------------------------------------  !!
    !!     1     2    11    12    13    14    23                            !!
    !!           ^                       ^                                  !!
    !!       3  4  5  6             15 16 17 18                             !!
    !!       7  8  9 10             19 20 21 22                             !!
    !! -------------------------------------------------------------------- !!

    implicit none

    integer::i,j,k
    integer::levgrid,levgrid2
    real(kind=8),dimension(1:levgrid)::tr1
    real(kind=8),dimension(1:levgrid2)::tr2
    real(kind=8),dimension(1:levgrid+levgrid2)::tr3
    logical,dimension(1:levgrid)::to1
    logical,dimension(1:levgrid)::tc1
    logical,dimension(1:levgrid+levgrid2)::to3
    logical,dimension(1:levgrid2)::to2
    logical,dimension(1:levgrid2)::tc2
    logical,dimension(1:levgrid+levgrid2)::tc3

    !! indicies
    integer::ind1,ind2

    k=0
    do i=1,levgrid
      to3(i+8*k)=to1(i)
      tr3(i+8*k)=tr1(i)
      tc3(i+8*k)=tc1(i)
      if(to1(i).eqv..true.)then
        to3(i+8*k)=.false.
        k=k+1
        to3(i+8*k-7:i+8*k)=to2(8*k-7:8*k)
        tr3(i+8*k-7:i+8*k)=tr2(8*k-7:8*k)
        tc3(i+8*k-7:i+8*k)=tc2(8*k-7:8*k)
      endif
    enddo

    return
  
  end subroutine mlink

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine writeski(skiname,nphoton,inputang,inputfov,cenx,ceny,&
                     &npix,d_ang,minw,maxw,npoint,amrrange,ind,n_piece)

    !! variables
      implicit none

      character(len=128)::skiname

      integer::nphoton
      character(len=10)::npstr

      real(kind=8),dimension(1:3)::inputang    
      real(kind=8),dimension(1:2)::inputfov
      real(kind=8)::theta,phi,omega
      character(len=10)::thestr,phistr,omestr
      real(kind=8)::fovx,fovy
      character(len=15)::fovxstr,fovystr
      real(kind=8)::cenx,ceny
      character(len=15)::cenxstr,cenystr
      integer::npix
      character(len=10)::npixstr
      real(kind=8)::d_ang !! Mpc
      character(len=10)::dangstr
      character(len=50)::istname

      real(kind=8)::minw,maxw
      character(len=10)::minwstr,maxwstr
      integer::npoint
      character(len=10)::npointstr

      real(kind=8),dimension(1:6)::amrrange
      real(kind=8)::minx,maxx
      real(kind=8)::miny,maxy
      real(kind=8)::minz,maxz
      character(len=15)::minxstr,maxxstr,minystr,maxystr,minzstr,maxzstr

      integer::n_piece,i,ind1,ind2,ind
      character(len=200)::filename

      theta=inputang(1) ; phi=inputang(2) ; omega=inputang(3)
      fovx=inputfov(1) ; fovy=inputfov(2)
      minx=amrrange(1) ; maxx=amrrange(2)
      miny=amrrange(3) ; maxy=amrrange(4)
      minz=amrrange(5) ; maxz=amrrange(6)

    open(unit=10,file=skiname,action='write',status='replace')
      write(10,'(A)')'<?xml version="1.0" encoding="UTF-8"?>'
      write(10,'(A)')'<!--SKIRT radiative transfer simulations - © 2012-2014 Astronomical Observatory, Ghent University-->'
      write(10,'(A)')'<skirt-simulation-hierarchy type="MonteCarloSimulation" format="6.1" producer="SKIRT v7.4 (git 828-dd92d35'//&
                    &' built on Oct 17 2017 at 15:34:26)" time="2017-10-19T14:54:14">'
      !! number of photons
      write(npstr,'(I10)')nphoton
      !! =================
      write(10,'(A)')'    <PanMonteCarloSimulation packages="'//adjustl(npstr)//'" minWeightReduction="1e4" minScattEvents="0" s'//&
                    &'cattBias="0.5" continuousScattering="false">'
      write(10,'(A)')'        <random type="Random">'
      write(10,'(A)')'            <Random seed="12341234"/>'
      write(10,'(A)')'        </random>'
      write(10,'(A)')'        <units type="Units">'
      write(10,'(A)')'            <ExtragalacticUnits fluxOutputStyle="Wavelength"/>'
      write(10,'(A)')'        </units>'
      write(10,'(A)')'        <instrumentSystem type="InstrumentSystem">'
      write(10,'(A)')'            <InstrumentSystem>'
      write(10,'(A)')'                <instruments type="Instrument">'
      !! angles and position
      write(thestr,'(F10.3)')theta ; write(phistr,'(F10.3)')phi
      write(omestr,'(F10.3)')omega
      write(fovxstr,'(F15.3)')fovx/real(n_piece)
      write(fovystr,'(F15.3)')fovy/real(n_piece)
      write(dangstr,'(F10.3)')d_ang
      !! ===================
        ind1=mod(ind,n_piece)
        ind2=ind/n_piece
        write(cenxstr,'(F15.3)')(ind1+0.5)/real(n_piece)*fovx-fovx/2d0
        write(cenystr,'(F15.3)')(ind2+0.5)/real(n_piece)*fovy-fovy/2d0
        write(npixstr,'(I10)')npix/n_piece
        write(istname,'(I0.5)')ind
        write(10,'(A)')'                    <SimpleInstrument instrumentName="'//adjustl(trim(istname))//'" distance="'//&
                      &adjustl(trim(dangstr))//&
                      &' Mpc" inclination="'//adjustl(trim(thestr))//' deg" azimuth="'//adjustl(trim(phistr))//' deg"'//&
                      &' positionAngle="'//adjustl(trim(omestr))//' deg" fieldOfViewX="'//adjustl(trim(fovxstr))//&
                      &' pc" pixelsX="'//adjustl(trim(npixstr))//'" centerX="'//adjustl(trim(cenxstr))//&
                      &' pc" '//'fieldOfViewY="'//adjustl(trim(fovystr))//' pc" pixelsY="'//adjustl(trim(npixstr))//&
                      &'" centerY="'//adjustl(trim(cenystr))//' pc"/>'
      write(10,'(A)')'                </instruments>'
      write(10,'(A)')'            </InstrumentSystem>'
      write(10,'(A)')'        </instrumentSystem>'
      write(10,'(A)')'        <wavelengthGrid type="PanWavelengthGrid">'
      !! spectral range and resolution
      write(minwstr,'(F10.3)')minw
      write(maxwstr,'(F10.3)')maxw
      write(npointstr,'(I10)')npoint
      !! =============================
      write(10,'(A)')'            <LogWavelengthGrid writeWavelengths="true" minWavelength="'//adjustl(trim(minwstr))//&
                    &' micron" maxWavelength="'//adjustl(trim(maxwstr))//' micron" points="'//adjustl(trim(npointstr))//&
                    &'"/>'
      write(10,'(A)')'        </wavelengthGrid>'
      write(10,'(A)')'        <stellarSystem type="StellarSystem">'
      write(10,'(A)')'            <StellarSystem emissionBias="0.5">'
      write(10,'(A)')'                <components type="StellarComp">'
      write(10,'(A)')'                    <SPHStellarComp filename="../part_ramski.txt" velocity="false" writeLuminosities="false">'
      write(10,'(A)')'                        <sedFamily type="SEDFamily">'
      write(10,'(A)')'                            <BruzualCharlotSEDFamily/>'
      write(10,'(A)')'                        </sedFamily>'
      write(10,'(A)')'                    </SPHStellarComp>'
      write(10,'(A)')'                </components>'
      write(10,'(A)')'            </StellarSystem>'
      write(10,'(A)')'        </stellarSystem>'
      write(10,'(A)')'        <dustSystem type="PanDustSystem">'
      write(10,'(A)')'            <PanDustSystem sampleCount="100" writeConvergence="false" writeDensity="false" writeDepthMap="'//&
                    &'false" writeQuality="false" writeCellProperties="false" writeCellsCrossed="false" emissionBias="0.5" emiss'//&
                    &'ionBoost="1" selfAbsorption="false" cycles="0" writeEmissivity="false" writeTemperature="false" writeISRF='//&
                    &'"false">'
      write(10,'(A)')'                <dustDistribution type="DustDistribution">'
      !! AMR cells range
      write(minxstr,'(F15.3)')minx ; write(maxxstr,'(F15.3)')maxx
      write(minystr,'(F15.3)')miny ; write(maxystr,'(F15.3)')maxy
      write(minzstr,'(F15.3)')minz ; write(maxzstr,'(F15.3)')maxz
      !! ===============
      write(10,'(A)')'                    <AdaptiveMeshDustDistribution minX="'//adjustl(trim(minxstr))//' pc" maxX="'//&
                    &adjustl(trim(maxxstr))//' pc" minY="'//adjustl(trim(minystr))//' pc" maxY="'//adjustl(trim(maxystr))//&
                    &' pc" minZ="'//adjustl(trim(minzstr))//' pc" maxZ="'//adjustl(trim(maxzstr))//' pc" densityUnits="1 Msun/pc3">'
      write(10,'(A)')'                        <adaptiveMeshFile type="AdaptiveMeshFile">'
      write(10,'(A)')'                            <AdaptiveMeshAsciiFile filename="../dust_ramski.txt"/>'
      write(10,'(A)')'                        </adaptiveMeshFile>'
      write(10,'(A)')'                        <components type="MeshDustComponent">'
      write(10,'(A)')'                            <MeshDustComponent densityIndex="0" multiplierIndex="-1" densityFraction="1">'
      write(10,'(A)')'                                <mix type="DustMix">'
      write(10,'(A)')'                                    <MeanZubkoDustMix writeMix="false" writeMeanMix="false"/>'
      write(10,'(A)')'                                </mix>'
      write(10,'(A)')'                            </MeshDustComponent>'
      write(10,'(A)')'                        </components>'
      write(10,'(A)')'                    </AdaptiveMeshDustDistribution>'
      write(10,'(A)')'                </dustDistribution>'
      write(10,'(A)')'                <dustGrid type="DustGrid">'
      write(10,'(A)')'                    <AdaptiveMeshDustGrid writeGrid="false"/>'
      write(10,'(A)')'                </dustGrid>'
      write(10,'(A)')'                <dustEmissivity type="DustEmissivity">'
      write(10,'(A)')'                    <GreyBodyDustEmissivity/>'
      write(10,'(A)')'                </dustEmissivity>'
      write(10,'(A)')'                <dustLib type="DustLib">'
      write(10,'(A)')'                    <AllCellsDustLib/>'
      write(10,'(A)')'                </dustLib>'
      write(10,'(A)')'            </PanDustSystem>'
      write(10,'(A)')'        </dustSystem>'
      write(10,'(A)')'    </PanMonteCarloSimulation>'
      write(10,'(A)')'</skirt-simulation-hierarchy>'
    close(10)

  end subroutine writeski

  !! ===================================================================== !!
  !! ===================================================================== !!

  subroutine cosmo(zz,omega_m,omega_L,H0,nn,d_ang)

    implicit none

    integer::i

    real*8::zz,omega_m,omega_L,d_ang,dz,H0
    integer::nn
    real*8,dimension(:),allocatable::xx
    real*8::cc=3d5

    dz=zz/real(nn) ; d_ang=0
    allocate(xx(0:nn))
    do i=0,nn
      d_ang=d_ang+dz/sqrt(omega_m*(1+(real(i+0.5)*dz))**3+omega_L)
    enddo
    d_ang=d_ang*cc/H0/(1.+zz)

  end subroutine cosmo

  !! ===================================================================== !!
  !! ===================================================================== !!

!! == SUBROUTINES ========================================================== !!
!! ========================================================================= !!





!! ========================================================================= !!
!! == FUNCTION ============================================================= !!

  !! ===================================================================== !!
  !! ===================================================================== !!

  function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
    real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
    dadtau = axp_tau*axp_tau*axp_tau *  &
         &   ( O_mat_0 + &
         &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
         &     O_k_0   * axp_tau )
    dadtau = sqrt(dadtau)
    return
  end function dadtau

  !! ===================================================================== !!
  !! ===================================================================== !!

  function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
    real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
    dadt   = (1.0D0/axp_t)* &
         &   ( O_mat_0 + &
         &     O_vac_0 * axp_t*axp_t*axp_t + &
         &     O_k_0   * axp_t )
    dadt = sqrt(dadt)
    return
  end function dadt

  !! ===================================================================== !!
  !! ===================================================================== !!

!! == FUNCTION ============================================================= !!
!! ========================================================================= !!
