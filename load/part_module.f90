module part_load
  real,allocatable,dimension(:,:)::b
  real(KIND=8),allocatable,dimension(:,:)::star_arr
  real(KIND=8),allocatable,dimension(:,:)::dm_arr
contains

  subroutine count_part(ndm_actual, nstar_actual, nsink_actual, repository, xmin, xmax, ymin, ymax, zmin, zmax)
    !--------------------------------------------------------------------------
    ! Ce programme calcule la carte de densite surfacique projetee
    ! des particules de matiere noire d'une simulation RAMSES. 
    ! Version F90 par R. Teyssier le 01/04/01.
    !--------------------------------------------------------------------------
    implicit none
  
    integer,INTENT(OUT)::ndm_actual, nstar_actual, nsink_actual
  
    integer::ncpu,ndim,npart,i,j,k,icpu,ipos,nstar
    integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel
    integer::ncpu_read
    real(KIND=8), INTENT(IN)::xmin,xmax,ymin,ymax,zmin,zmax
    integer::imin,imax,jmin,jmax,kmin,kmax,lmin!,npart_actual
    real(KIND=8)::xxmin,xxmax,yymin,yymax,dx,dy,deltax,boxlen
    real(KIND=8),dimension(:,:),allocatable::x
    real(KIND=8),dimension(:)  ,allocatable::age
    integer,dimension(:)  ,allocatable::id
    character(LEN=5)::nchar,ncharcpu
    character(LEN=80)::ordering
    character(LEN=128)::nomfich,repository
    logical::ok, ok_part, ok_sink, ok_star, ok_dm
    
    ! CPU list
    integer::impi,ndom,bit_length,maxdom
    integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
    real(KIND=8),dimension(1:8)::bounding,bounding_min,bounding_max
    real(KIND=8)::dkey,order_min,dmax
    real(kind=8)::xx,yy,zz,aexp
    real(kind=8),dimension(:),allocatable::bound_key
    logical,dimension(:),allocatable::cpu_read
    !  integer,dimension(:),allocatable, INTENT(OUT)::cpu_list 
    !  -> No!, you can not return an assumed-size array. 
      integer,dimension(:),allocatable::cpu_list 
    !-----------------------------------------------
    ! Lecture du fichier particules au format RAMSES
    !-----------------------------------------------
    ipos=INDEX(repository,'output_')
    nchar=repository(ipos+7:ipos+13)
  
    write(*,*)'Reading info'
  
    nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
    write(*,*)nomfich
    open(unit=10,file=nomfich,form='formatted',status='old')
    read(10,'(13x,I11)')ncpu
    read(10,'(13x,I11)')ndim
    read(10,'(13x,I11)')levelmin
    read(10,'(13x,I11)')levelmax
    read(10,*)
    read(10,*)
    read(10,*)
  
    read(10,'(13x,E23.15)')boxlen
    read(10,'(13x,E23.15)')
    read(10,'(13x,E23.15)')aexp
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
  
    read(10,'(14x,A80)'),ordering
    write(*,'(" ordering type=",A20)'),TRIM(ordering)
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
    ! Get cpu list
    !-----------------------
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
             call hilbert3d(idom(i),jdom(i),kdom(i),bounding(1),bit_length,1)
             order_min=bounding(1)
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
  
  ! read part
  !  npart_actual=0
    ndm_actual=0
    nstar_actual=0
    nsink_actual=0
  
    do k=1,ncpu_read
       icpu=cpu_list(k)
       call title(icpu,ncharcpu)
       nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
       open(unit=1,file=nomfich,status='old',form='unformatted')
  !     write(*,*)'Processing file '//TRIM(nomfich)
       read(1)ncpu2
       read(1)ndim2
       read(1)npart2
       read(1)
       read(1)nstar
       read(1)
       read(1)
       read(1)
  
  !     allocate(m(1:npart2))
       allocate(id(1:npart2))
       if(nstar>0)then ! age and id to distinguish star, DM, and sink.
          allocate(age(1:npart2))
       endif
       allocate(x(1:npart2,1:ndim2)) 
       ! Position to select particles insdie the region of interest
       ! Read position
  !     do i=1,ndim
       read(1)x(1:npart2,1)
       read(1)x(1:npart2,2)
       read(1)x(1:npart2,3)
  !     end do
  
       ! Skip velocity
       do i=1,ndim
          read(1) !age
       end do
       ! Skip mass
       read(1) !age
       if(nstar>0)then
          read(1)id
          read(1) ! Skip level
          read(1)age
  	! ignore metal
       endif
  
       close(1)
  
       do i=1,npart2
            ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                &   x(i,3)>=zmin.and.x(i,3)<=zmax)
  
            if(ok_part.and.(age(i).eq.0.0d0).and.(id(i)>0)) then
  	      ndm_actual = ndm_actual + 1
            elseif(ok_part.and.(age(i).ne.0.0d0).and.(id(i)>0)) then
                nstar_actual = nstar_actual + 1
  	  elseif(ok_part.and.(age(i).eq.0.0d0).and.(id(i)<0)) then
                nsink_actual = nsink_actual + 1
            endif
       enddo
  !     npart_actual = npart_actual + npart2
       deallocate(x,id)
       if(nstar>0)deallocate(age)
    enddo
  
    end subroutine
  
  
    subroutine load_part(nstar_actual, ndm_actual, nsink_actual, repository, &
               & xmin, xmax, ymin, ymax, zmin, zmax)
  
    real(KIND=8),dimension(:,:),allocatable::x,v
    real(KIND=8),dimension(:)  ,allocatable::m,age,metal
    integer,dimension(:)  ,allocatable::id
    character(LEN=5)::nchar,ncharcpu
    character(LEN=80)::ordering
    character(LEN=128)::nomfich,repository
  
    real(KIND=8), INTENT(IN)::xmin,xmax,ymin,ymax,zmin,zmax
    integer,INTENT(IN)::ndm_actual, nstar_actual, nsink_actual
  
  
    ! CPU list
    integer::impi,ndom,bit_length,maxdom
    integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
    real(KIND=8),dimension(1:8)::bounding,bounding_min,bounding_max
    real(KIND=8)::dkey,order_min,dmax
    real(kind=8)::xx,yy,zz,aexp
    real(kind=8),dimension(:),allocatable::bound_key
    logical,dimension(:),allocatable::cpu_read
    integer,dimension(:),allocatable::cpu_list
  
  !  real(KIND=8),INTENT(INOUT),dimension(nstar_actual,10)::star_arr
  !  real(KIND=8),INTENT(INOUT),dimension(ndm_actual,8)::dm_arr
  !  real(KIND=8),INTENT(INOUT),dimension(nsink_actual,8)::sink_arr
  
  !  real(KIND=8),INTENT(OUT),dimension(:,:),allocatable::sink_arr
  
    logical::ok, ok_part, ok_sink, ok_star, ok_dm
    integer::ncpu,ndim,npart,i,j,k,icpu,ipos,nstar
    integer::i_dm, i_star
    integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel
    
  !  allocate(star_arr(1:nstar_actual,10)) ! id,x,y,z,vx,vy,vz,mass,age,metal
  !  allocate(dm_arr(1:ndm_actual,8)) ! id,x,y,z,vx,vy,vz,mass
  
  !  if(nsink_actual > 0) allocate(sink_arr(1:nsink_actual,8)) ! id,x,y,z,vx,vy,vz,mass
  
    ! CPU list
  
    ipos=INDEX(repository,'output_')
    nchar=repository(ipos+7:ipos+13)
  
    write(*,*)'Reading info'
  
    nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
    write(*,*)nomfich
    open(unit=10,file=nomfich,form='formatted',status='old')
    read(10,'(13x,I11)')ncpu
    read(10,'(13x,I11)')ndim
    read(10,'(13x,I11)')levelmin
    read(10,'(13x,I11)')levelmax
    read(10,*)
    read(10,*)
    read(10,*)
  
    read(10,'(13x,E23.15)')boxlen
    read(10,'(13x,E23.15)')
    read(10,'(13x,E23.15)')aexp
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
  
    read(10,'(14x,A80)'),ordering
    write(*,'(" ordering type=",A20)'),TRIM(ordering)
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
    ! Get cpu list
    !-----------------------
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
             call hilbert3d(idom(i),jdom(i),kdom(i),bounding(1),bit_length,1)
             order_min=bounding(1)
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
  
    ! Store particles
    i_dm = 0
    i_star = 0
    i_sink = 0
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
       read(1)
       read(1)
       read(1)
  
       allocate(x(1:npart2,1:ndim2)) 
       allocate(m(1:npart2))
       allocate(v(1:npart2,1:ndim2)) 
       allocate(id(1:npart2))
  !    if(nstar>0)then ! age and id to distinguish star, DM, and sink.
          allocate(age(1:npart2))
          allocate(metal(1:npart2))
  !    endif
       ! Position to select particles insdie the region of interest
       ! Read position
       do i=1,ndim
         read(1)x(1:npart2,i)
  !	read(1)x(1:npart2,2)
  !	read(1)x(1:npart2,3)
  !       read(1)m
  !       x(1:npart2,0)=m/boxlen
       end do
  
       ! Read velocity
       do i=1,ndim
         read(1)v(1:npart2,i)
  !	read(1)v(1:npart2,2)
  !	read(1)v(1:npart2,3)
       end do

     ! Skip mass
       read(1)m
       read(1)id
       if(nstar>0)then
          read(1) ! Skip level
          read(1)age
          read(1)metal
       endif
       close(1)
   
       do i=1,npart2
          ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
              &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
              &   x(i,3)>=zmin.and.x(i,3)<=zmax)
   
          if(ok_part.and.(age(i).eq.0.0d0))then ! Sink particles are also considered as DM.
             dm_arr(i_dm,1) = x(i,1)
             dm_arr(i_dm,2) = x(i,2)
             dm_arr(i_dm,3) = x(i,3)
             dm_arr(i_dm,4) = v(i,1)
             dm_arr(i_dm,5) = v(i,2)
             dm_arr(i_dm,6) = v(i,3)
             dm_arr(i_dm,7) = m(i)
             dm_arr(i_dm,8) = id(i)
             i_dm = i_dm + 1
          elseif(ok_part.and.(age(i).ne.0.0d0).and.(id(i)>0))then
             star_arr(i_star,1) = x(i,1)
             star_arr(i_star,2) = x(i,2)
             star_arr(i_star,3) = x(i,3)
             star_arr(i_star,4) = v(i,1)
             star_arr(i_star,5) = v(i,2)
             star_arr(i_star,6) = v(i,3)
             star_arr(i_star,7) = m(i)
             star_arr(i_star,8) = id(i)
             star_arr(i_star,9) = age(i)
             star_arr(i_star,10) = metal(i)
             i_star = i_star + 1
!  	elseif(ok_part.and.(age(i).eq.0.0d0).and.(id(i)<0))then
!  	   sink_arr(i_sink,1) = x(i,1)
!  	   sink_arr(i_sink,2) = x(i,2)
!  	   sink_arr(i_sink,3) = x(i,3)
!  	   sink_arr(i_sink,4) = v(i,1)
!  	   sink_arr(i_sink,5) = v(i,2)
!  	   sink_arr(i_sink,6) = v(i,3)
!  	   sink_arr(i_sink,7) = m(i)
!	   sink_arr(i_sink,8) = id(i)
!           i_sink = i_sink + 1
          endif
       end do
       deallocate(x,v,m,id)
       if(nstar>0)deallocate(age,metal)
    end do

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
  
  
end module part_load
  
  
  
  
  
  
  
  
  
