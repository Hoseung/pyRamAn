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

end module ramski_mods

  !! ====================================================================== !!
  !! ====================================================================== !!

module tree_mods

  contains

  subroutine count_tree(fn, n_halos_all, flist_index, slist_index, &
                      & nsteps, big_run, n_all_fathers, n_all_sons)
    implicit none
    character(LEN=256), INTENT(IN)::fn

    integer::nfathers, nsons
    integer, dimension(:), allocatable::nb_of_halos, nb_of_subhalos
    integer, INTENT(OUT)::flist_index, slist_index, n_halos_all, nsteps

    integer(KIND=4)::i,j, nhals_now
    logical, INTENT(IN)::big_run

    integer,INTENT(OUT)::n_all_fathers,n_all_sons

    n_all_fathers=0 ; n_all_sons=0

    open(unit=1,file=fn,status='old',form='unformatted')

    read(1)nsteps
    allocate(nb_of_halos(1:nsteps))
    allocate(nb_of_subhalos(1:nsteps))
    read(1)nb_of_halos(1:nsteps),nb_of_subhalos(1:nsteps)

    n_halos_all = sum(nb_of_halos) + sum(nb_of_subhalos)
    read(1)!aexp_arr(1:nsteps)
    read(1)!omega_t_arr(1:nsteps)
    read(1)!age_univ_arr(1:nsteps)
    flist_index=0
    slist_index=0
    do i=1,nsteps
        nhals_now=nb_of_halos(i)+nb_of_subhalos(i)
        do j=1,nhals_now

            read(1)!id_tmp!i_arr(idx,2) !id!id
            read(1)!i_arr(idx,3) !bushID
            read(1)!i_arr(idx,4) !st
            read(1)!i_arr(idx,5:9) ! hosts
            read(1)!f_arr(idx,1) ! m
            read(1)!macc ! macc alone is double
            read(1)!f_arr(idx,3:5)!xp
            read(1)!f_arr(idx,6:8)!vp
            read(1)!f_arr(idx,9:11)!lp
            read(1)!f_arr(idx,12:15)!abc
            read(1)!f_arr(idx,16:18)!energy
            read(1)!f_arr(idx,19)!spin
            read(1)nfathers!nfathers
            n_all_fathers=n_all_fathers+nfathers
            flist_index = flist_index + nfathers

            read(1)!Father ID
            read(1)!Father Mass
            read(1)nsons
            n_all_sons=n_all_sons+nsons
            if (nsons .gt. 0) then
                read(1)!id_tmp!son ID
                !write(*,*)id_tmp
            endif
            slist_index = slist_index + nsons
            read(1)! (r,m,t)_vir, c_vel
            read(1)! profile
            if (.not. big_run) then
                read(1)!id_tmp! particle ID
            endif

        enddo
    enddo
    deallocate(nb_of_halos, nb_of_subhalos)
    close(1)

  end subroutine

  subroutine load_tree(fn, fatherID, fatherIDx, sonID, fatherMass, &
                   & i_arr, f_arr, aexp_arr, omega_t_arr, age_univ_arr, &
                   & n_halos_all, n_all_fathers, n_all_sons, big_run, nsteps)

    implicit none
    character(LEN=256), INTENT(IN)::fn
    logical, INTENT(IN)::big_run
    integer, INTENT(IN)::n_all_sons, n_all_fathers, n_halos_all, nsteps

    integer(KIND=4):: nsons, flist_index, slist_index
    integer(KIND=4):: i,j, k,nhals_now, n_fathers, idx_old, nhals_old, idx
    real(KIND=8)::macc

    integer, dimension(:), allocatable::fid_tmp, nb_of_halos, &
                                      & nb_of_subhalos, sonID_tmp

    ! +1 so that in python the array content start from index 1.
    integer, dimension(:), allocatable, INTENT(OUT)::fatherID, fatherIDx
    integer, dimension(:), allocatable, INTENT(OUT)::sonID
    real(KIND=8), dimension(:), allocatable, INTENT(OUT) ::fatherMass
    integer(KIND=4), dimension(:,:), allocatable, INTENT(OUT) ::i_arr
    real(KIND=8), dimension(:,:), allocatable, INTENT(OUT) ::f_arr
    real(KIND=8), dimension(1:nsteps), INTENT(OUT) ::aexp_arr, &
                                    & omega_t_arr, age_univ_arr

    integer::n_fathers_max

    n_fathers_max = 5000
    allocate(fid_tmp(1:n_fathers_max)) ! At most 100 fathers for a halo.
    allocate(sonID_tmp(1:n_fathers_max)) ! At most 100 fathers for a halo.

    open(unit=1,file=fn,status='old',form='unformatted')

    read(1)
    allocate(nb_of_halos(1:nsteps))
    allocate(nb_of_subhalos(1:nsteps))
    read(1)nb_of_halos(1:nsteps),nb_of_subhalos(1:nsteps)

    read(1)aexp_arr(1:nsteps)
    read(1)omega_t_arr(1:nsteps)
    read(1)age_univ_arr(1:nsteps)

    flist_index=1
    slist_index=1
    nhals_old = 0
    idx=1

    allocate(i_arr(1:n_halos_all,1:15))
    allocate(f_arr(1:n_halos_all,1:25))
    allocate(fatherID(1:n_all_fathers+1),fatherIDx(1:n_all_fathers+1))
    allocate(sonID(1:n_all_sons))
    allocate(fatherMass(1:n_all_fathers))

    do i=1,nsteps
        nhals_now=nb_of_halos(i)+nb_of_subhalos(i)
        idx_old = idx - nhals_old
        i_arr(idx:idx+nhals_now, 11) = i ! current step
        do j=1,nhals_now
            i_arr(idx,1) = idx
            read(1)i_arr(idx,2) !id!id
            read(1)i_arr(idx,3) !bushID
            read(1)i_arr(idx,4) !st
            read(1)i_arr(idx,5:9) ! hosts
            read(1)f_arr(idx,1) ! m
            read(1)macc ! macc alone is double
            f_arr(idx,2) = real(macc)
            read(1)f_arr(idx,3:5)!xp
            read(1)f_arr(idx,6:8)!vp
            read(1)f_arr(idx,9:11)!lp
            read(1)f_arr(idx,12:15)!abc
            read(1)f_arr(idx,16:18)!energy
            read(1)f_arr(idx,19)!spin
            read(1)n_fathers!nfathers
            if (n_fathers .gt. n_fathers_max) then 
                write(*,*) "Increase n_fathers_max above ",n_fathers
            endif
            read(1)fid_tmp(1:n_fathers)
            fatherID(flist_index:flist_index+n_fathers-1)=fid_tmp(1:n_fathers)
            read(1)fatherMass(flist_index:flist_index+n_fathers-1)
            if (i > 0) then
                ! from the previous step.
                do k=1,n_fathers
                    if (fid_tmp(k) > 0) then
                        fid_tmp(k) = i_arr(idx_old+fid_tmp(k),1)
                    endif
                enddo
            endif
            i_arr(idx,10) = n_fathers
            fatherIDx(flist_index:flist_index+n_fathers-1)=fid_tmp(1:n_fathers)
            i_arr(idx,13) = flist_index
            flist_index = flist_index + n_fathers

            read(1)nsons
            i_arr(idx,14) = nsons
            if (nsons .gt. 0) then
            !read(1)sonID(slist_index:slist_index+nsons-1)!sonID_tmp(1:nsons-1)
                read(1)sonID_tmp(1:nsons-1)
                ! id into idx.
                if (i .lt. nsteps) then
                    sonID(slist_index:slist_index+nsons-1)=&
                      &sonID_tmp(1:nsons-1)+idx_old+nhals_now+nhals_old-1
                    i_arr(idx,15) = slist_index
                    !slist_index = slist_index+nsons
                endif
            endif
            slist_index = slist_index + nsons
            read(1)f_arr(idx,20:23) ! virial
            read(1)f_arr(idx,24:25) ! rho
            if (.not. big_run) then
                read(1)i_arr(idx,12) ! np
            endif
            idx = idx +1
        enddo
        nhals_old = nhals_now


    enddo
    deallocate(nb_of_halos, nb_of_subhalos)
    deallocate(fid_tmp)
    close(1)

  end subroutine

end module tree_mods


!! == MODULES ============================================================== !!
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
      integer::nx,ny,nz!,mindx,mindy,mindz
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
      !real(kind=8)::minx,miny,minz

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
