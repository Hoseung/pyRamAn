!! ========================================================================= !!
!! === RAMSKI - RAMSES wearing SKIRT ======================================= !!
!! === This code makes the files required for SKIRT simulation. ============ !!
!! === Written by Jongwon Park (Yonsei University) ========================= !!
!! ========================================================================= !!

!! ========================================================================= !!
!! == MAIN PROGRAM ========================================================= !!

subroutine ramski_v4(dir, outdir, fn, gx, nstart)

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
    real*8,dimension(1:3), intent(in)::gx     !! galaxy center in code unit
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
    character(len=200), INTENT(IN)::dir, outdir       !! directory of RAMSES snapshot
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
    character(len=500)::write_file
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
    character(len=256), INTENT(in)::fn

    !! count_tree
    integer::n_halos_all,nsteps,flist_index,slist_index
    logical::big_run
    integer::n_all_fathers,n_all_sons
    !! load_tree
    integer,dimension(:),allocatable::fatherID,fatherIDx
    integer,dimension(:),allocatable::sonID
    real*8,dimension(:),allocatable::fatherMass
    integer,dimension(:,:),allocatable::i_arr
    real*8,dimension(:,:),allocatable::f_arr
    real*8,dimension(:),allocatable::aexp_arr,omega_t_arr,age_univ_arr

    !! galaxies in the last snapshot
    integer, INTENT(in)::nstart
    integer::nend, tstart
    integer::ngal !! used to read the tree
    real*8::masscut=1d10

    !! history of each galaxy
    real*8,dimension(:,:,:),allocatable::gal !! mass and position
    integer::ind,ind1,ind2,ind3,ind4
    real*8,dimension(:),allocatable::fMass
                        !! father mass - fine the most massive one
    integer,dimension(:),allocatable::find
                        !! father index
    integer::n_fathers

    !! find maximum father mass
    integer,dimension(1:1)::tmp_ind

    !! unit_l
    character(len=80)::GMGM

    !! faceon and edgeon
    integer::faceon=0,edgeon=0

    !! define parameter here instead of using script
    !dir='/storage1/NewHorizon/OUTPUT_DIR/' !! RAMSES snapshot
    !outdir = dir//'output/'
    !outdir='/home/jongwon/analysis/mock/for_nature' !! output directory
    !fn='/storage1/NewHorizon/GalaxyMaker/gal/tree.dat' !! path to the tree

    lmin=12     !! minimum AMR level (must be smaller than lmax)
    lmax=21     !! maximum AMR level (must be equal to lmax in RAMSES)
    Tdust=10000 !! dust survival temperature
    fdust=0.4   !! dust fraction - 0.4 from Saftly et al. 2015
    !smooth=50   !! stellar particle smoothing length
    !extr=5      !! extract one particle out of a group of partices
    !extr=2
    !nph=50000   !! number of photons to be used in SKIRT
    !nph=2d5
    ps=0.04     !! plate scale - arcsec/pixel
    minw=0.1    !! SED minimum wavelength (micron)
    maxw=1.4    !! SED maximum wavelength (micron)
    npoint=120  !! number of points in SED
    fov=60      !! field of view (kpc)
                !! This will be slightly adjusted according to place scale.
    inputfov(1)=fov ; inputfov(2)=fov
    !n_piece=5   !! divide the field of view into n_piece x n_piece equivalent
                !! regions (if n_piece=4, FoV splits into 4x4=16 aread)
                !! (Use this for high quality image.)
    n_piece=1
    lbox=100    !! simulation box size (Mpc)
    h0=70.4     !! Hubble constant
    masscut=1d10
    !nstart=509  !! first snapshot number in the tree
    !nend=509    !! last snapshot to be considered
    nend = nstart
    tstart = nstart
    z_assume=1.00
    faceon=1    !! if you want this output, set any non-zero integer
    edgeon=1    !! if you want this output, set any non-zero integer

    !gal(ngal,j,2:4)=f_arr(ind,3:5) !! position - physical Mpc


    !! READ MORTON LOOK-UP TABLE
    open(unit=10,file='LUT.txt',form='formatted',status='old')
    do i=0,2**20-1
      read(10,*)bm(i)
    enddo
    close(10)

    !! simulation time - for paticle ages
      !! path to the last RAMSES snapshot - just to read a info file
      write(nsnap_char,'(I0.5)')nend
      path=trim(dir)//'output_'//nsnap_char
      !! read info
      call rd_info(path,info)
      omega_m=info%omega_m ; omega_l=info%omega_l ; omega_k=info%omega_k
        !! keep it simple

    !! Allocate look-up tables
    allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
    allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
    !! Compute Friedman model look up table
    call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
         & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)

    !! READ the TREE
    !! count_tree
    !call count_tree(fn,n_halos_all,flist_index,slist_index,nsteps,big_run,&
    !               &n_all_fathers,n_all_sons)
    !! load_tree
    !call load_tree(fn,fatherID,fatherIDx,sonID,fatherMass,&
    !              &i_arr,f_arr,aexp_arr,omega_t_arr,age_univ_arr,&
    !              &n_halos_all,n_all_fathers,n_all_sons,big_run,nsteps)
    !  f_arr(:,1)=f_arr(:,1)*1d11 !! mass

    !! galaxies in the last snapshot
    !ngal=0
    !do i=1,n_halos_all
    !  if( i_arr(i,11).eq.(nend-tstart+1) .and. f_arr(i,1).ge.masscut)then
    !    ngal=ngal+1
    !    if(ngal.eq.1)ind1=i
    !    ind2=i
    !  endif
    !enddo
    !allocate(gal(ngal,nstart:nend,1:4)) ; gal=0

    
    !ngal=0
    !do i=ind1,ind2 !! do not need to run over the all haloes in the tree
    !  if( i_arr(i,11).eq.nend-tstart+1 .and. f_arr(i,1).ge.masscut )then
    !    ngal=ngal+1
    !    ind=i
    !    do j=nend,nstart,-1
    !      gal(ngal,j,1)=f_arr(ind,1)     !! mass
    !      gal(ngal,j,2:4)=f_arr(ind,3:5) !! position - physical Mpc
!    
!          n_fathers=i_arr(ind,10)      !! number of fathers
!          ind3=i_arr(ind,13)           !! flist_index
!          ind4=ind3+n_fathers-1        !! note the subtraction
!
!          allocate(fMass(1:n_fathers),find(1:n_fathers))
!          fMass(1:n_fathers)=fatherMass(ind3:ind4)               
!          find(1:n_fathers)=fatherIDx(ind3:ind4)
!
!          tmp_ind=maxloc(fMass)
!          ind=find(tmp_ind(1))-1
!
!          deallocate(fMass,find)
!        enddo
!      endif
!    enddo

    !! MAIN LOOP
    a = nstart

      !! store the snapshot number as a string
      write(nsnap_char,'(I0.5)')a

      !! path to the RAMSES snapshot
      path=trim(dir)//'output_'//nsnap_char

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

        !! store the galaxy number as a string
        write(ngal_char,'(I0.5)')i

        !gx(1:3)=gal(i,a,2:4)*3.086d24/unit_l+0.5
        print *, "Galaxy center", gx(1), gx(2), gx(3)
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
        print *, "amr2cell done"
        !! read particles
        !call part2cube(path,arr_range,part,npart)
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
        print *, "morton done"
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
        print *, "mlink done"
       !! calculate particle age
       ! do j=1,npart
       !   do k=1,n_frw-1
       !     if(part(j)%ap.le.tau_frw(k) .and. part(j)%ap.gt.tau_frw(k+1))then
       !       part(j)%ap =                                                &
       !       & ( t_frw(k+1) - t_frw(k) ) / ( tau_frw(k+1)-tau_frw(k) ) * &
       !       & ( part(j)%ap-tau_frw(k) ) + t_frw(k)
       !       part(j)%ap=(lbt-part(j)%ap)*t_universe !! age in yr
       !       goto 104
       !     endif
       !   enddo
       !   104 continue
       ! enddo

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
        print *, "writing..."
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
      !  write_file=trim(outdir)//'/'//ngal_char//'/'//nsnap_char//&
      !            &'/part_ramski.txt'  
      !  open(unit=10,file=write_file,action='write',status='replace')
      !    do j=1,npart,extr
      !      write(10,*)(part(j)%xp(:)-gx(:))*unit_l/3.086d18,smooth,&
      !               & part(j)%mp*real(extr)*unit_m/1.989d33,&
      !                 part(j)%zp,part(j)%ap
      !    enddo
      !  close(10)

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

        deallocate(cell,arr_mort,ncell_lev)
        deallocate(dd1,dd2,dd3,rr1,rr2,rr3,cc1,cc2,cc3)
        inputfov=inputfov/1d3 !! pc
        fov_assume=fov_assume/1d3

end subroutine ramski_v4

!! == MAIN PROGRAM ========================================================= !!
!! ========================================================================= !!
! ========================================================================= !!
