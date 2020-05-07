module amr_commons
    integer,parameter::qdp=kind(1.0_16) ! real*16, QUADHILBERT
    !integer,parameter::qdp=kind(1.0_8) ! real*8
    integer,parameter::dp=kind(1.0D0) ! real*8
    real(qdp),allocatable,dimension(:)::bound_key
    logical::verbose=.true.
end module

! Load amr head with QUADHILBERT

program load_amr
    use amr_commons
    integer::nboundary2, ndomain
    integer::ncpu2,ndim2,nx2,ny2,nz2,ngridmax2,nlevelmax2
    integer::noutput2,iout2,ifout2,ilun!,info
    integer::MAXOUT=1000
    real(dp),dimension(1:1000)::aout2=1.1d0
    real(dp),dimension(1:1000)::tout2=0.0d0
    logical::ok,debug,simple_boundary
    character(LEN=128)::ordering2
    character(LEN=80)::fileloc
    integer,parameter::tag=1100
    
    call get_command_argument(1,fileloc)
    write(*,*)TRIM(fileloc)
    
    ! USE
    ! ./load_amr /home/hoseung/Work/data/FORNAX/output_00113/amr_00113.out00001
    !

    debug=.true.
    simple_boundary=.false.
    ndomain=2400
    allocate(bound_key (0:ndomain))

    !!!!! Taken from init_amr
           
     inquire(file=fileloc, exist=ok)
     if(.not. ok)then
        write(*,*)'Restart failed:'
        write(*,*)'File '//TRIM(fileloc)//' not found'
        !call clean_stop
     end if
     if(debug)write(*,*)'amr.tmp opened for processor ',myid
     open(unit=ilun,file=fileloc,form='unformatted')
     ! Read grid variables
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)nx2,ny2,nz2
     read(ilun)nlevelmax2
     read(ilun)ngridmax2
     read(ilun)nboundary2
     read(ilun)ngrid_current
     read(ilun)boxlen
     if(ncpu2.ne.ncpu)then
        write(*,*)'Number of processes not compatible'
        write(*,*)'ncpu should be set equal to',ncpu2
        !call clean_stop
     end if
     ! Read time variables
     read(ilun)noutput2,iout2,ifout2
     if(noutput2>MAXOUT)then
       write(*,*) 'Error: noutput>MAXOUT'
       !call clean_stop
     end if
     read(ilun)tout2(1:noutput2)
     read(ilun)aout2(1:noutput2)
     ! Check compatibility with current parameters
     if((ndim2.ne.ndim).or.(nx2.ne.nx).or.(ny2.ne.ny).or.(nz2.ne.nz).or.&
          & (nboundary2.ne.nboundary).or.(nlevelmax2>nlevelmax).or.&
          & (ngrid_current>ngridmax).or.(noutput2>noutput) )then
        write(*,*)'File amr.tmp is not compatible with namelist'
        write(*,*)'         ndim   nx   ny   nz nlevelmax noutput   ngridmax nboundary'
        write(*,'("amr.tmp  =",4(I4,1x),5x,I4,4x,I4,3x,I8)')&
             & ndim2,nx2,ny2,nz2,nlevelmax2,noutput2,ngrid_current,nboundary2
        write(*,'("namelist =",4(I4,1x),5x,I4,4x,I4,3x,I8)')&
             & ndim ,nx ,ny ,nz ,nlevelmax ,noutput, ngridmax     ,nboundary
        if(myid==1)write(*,*)'Restart failed'
        !call clean_stop
     end if
     ! Old output times
     !tout(1:noutput2)=tout2(1:noutput2)
     !aout(1:noutput2)=aout2(1:noutput2)
     iout=iout2
     ifout=ifout2
     read(ilun)t
     read(ilun)!dtold(1:nlevelmax2)
     read(ilun)!dtnew(1:nlevelmax2)
     read(ilun)nstep,nstep_coarse
     nstep_coarse_old=nstep_coarse
     read(ilun)const,mass_tot_0,rho_tot
     read(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
     read(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
     read(ilun)mass_sph

     if(myid==1)write(*,*)'Restarting at t=',t,' nstep_coarse=',nstep_coarse

          ! Read levels variables
     read(ilun)!headl(1:ncpu,1:nlevelmax2)
     read(ilun)!taill(1:ncpu,1:nlevelmax2)
     read(ilun)!numbl(1:ncpu,1:nlevelmax2)
     read(ilun)!numbtot(1:10,1:nlevelmax2)
     ! Read boundary linked list
     if(simple_boundary)then
        read(ilun)!headb(1:nboundary,1:nlevelmax2)
        read(ilun)!tailb(1:nboundary,1:nlevelmax2)
        read(ilun)!numbb(1:nboundary,1:nlevelmax2)
     end if
     ! Read free memory
     read(ilun)headf,tailf,numbf,used_mem,used_mem_tot
     headf=ngrid_current+1
     tailf=ngridmax
     numbf=ngridmax-ngrid_current
     ! Read cpu boundaries
     read(ilun)ordering2
     !if(ordering2.ne.ordering)then
     !   if(myid==1)write(*,*)'Ordering is uncompatible'
     !   !call clean_stop
     !endif
    ! If ordering == 'hilbert'
     read(ilun)bound_key(0:ndomain)
    

     !!!!! Taken from init_amr
     ! print in ASCII
     write(*,'("   DOMAIN   ind_min                 ind_max")')
     do idom=1,ndomain
        write(*,'(I8,1X,E23.15,1X,E23.15)')idom,bound_key(idom-1),bound_key(idom)
     end do
    
     !!!!! Taken from dump_info
end program

