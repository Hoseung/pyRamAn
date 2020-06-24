
module readhtm
   implicit none
   logical         :: galaxy, dp=.false.
   integer         :: ihalo, nhalo, npart_tot, ipart

   integer(kind=4), dimension(:,:), allocatable  :: integer_table
   real(kind=4),    dimension(:,:), allocatable  :: real_table
   real(kind=8),    dimension(:,:), allocatable  :: real_table_dp

   integer(kind=4), dimension(:), allocatable :: part_ids

contains
!#########################################################
   subroutine read_halo_brick(unitfile)!,aexp)
!#########################################################
      implicit none
      integer(kind=4) :: unitfile
      !real(kind=4)    :: aexp

      read(unitfile) integer_table(ihalo,1) ! nparts
      npart_tot = npart_tot + integer_table(ihalo,1)
      read(unitfile)                ! ID of members

      read(unitfile) integer_table(ihalo,2) ! my_number
      read(unitfile) !integer_table(ihalo,3) ! my_timestep

      ! change the time step number to match the number of branching (time resolution)
      ! you decided your tree is going to have
      !integer_table(ihalo,3) = 0 ! meaningless in this program (TK)
      read(unitfile) integer_table(ihalo,3),integer_table(ihalo,4),integer_table(ihalo,5),&
            integer_table(ihalo,6),integer_table(ihalo,7)
      ! level, hosthalo, hostsub, nbsub, nextsub

      ! Old bug if structure tree does not exist, make, shouldn't occur with
      if(integer_table(ihalo,4).le.0) integer_table(ihalo,4) = integer_table(ihalo,1)
      if(integer_table(ihalo,3).le.0) integer_table(ihalo,3) = 1
      if(.not. dp) then
         !real_table(ihalo,1) = aexp ! aexp
         read(unitfile) real_table(ihalo,1) ! m
         read(unitfile) real_table(ihalo,2:4) ! x
         read(unitfile) real_table(ihalo,5:7) ! v
         read(unitfile) real_table(ihalo,8:10) ! L
         read(unitfile) real_table(ihalo,11:14) ! r, a, b, c
         read(unitfile) real_table(ihalo,15:17) ! ek, ep, et
         read(unitfile) real_table(ihalo,18) ! spin
         if(galaxy) then
            read(unitfile) real_table(ihalo,19:21)
            read(unitfile) real_table(ihalo,22:25)
            read(unitfile) real_table(ihalo,26:27)
            read(unitfile) ! some additional surface profiles data ...
            read(unitfile)
            read(unitfile)
         else
            read(unitfile) real_table(ihalo,20:23) ! rvir, mvir, tvir, cvel
            read(unitfile) real_table(ihalo,24:25) ! rho0, rc
         end if
      else
         !real_table_dp(ihalo,1) = aexp ! aexp
         read(unitfile) real_table_dp(ihalo,1) ! m
         read(unitfile) real_table_dp(ihalo,2:4) ! x
         read(unitfile) real_table_dp(ihalo,5:7) ! v
         read(unitfile) real_table_dp(ihalo,8:10) ! L
         read(unitfile) real_table_dp(ihalo,11:14) ! r, a, b, c
         read(unitfile) real_table_dp(ihalo,15:17) ! ek, ep, et
         read(unitfile) real_table_dp(ihalo,18) ! spin
         if(galaxy) then
            read(unitfile) real_table_dp(ihalo,19:21)
            read(unitfile) real_table_dp(ihalo,22:25)
            read(unitfile) real_table_dp(ihalo,26:27)
            read(unitfile) ! some additional surface profiles data ...
            read(unitfile)
            read(unitfile)
         else
            read(unitfile) real_table_dp(ihalo,19:22) ! rvir, mvir, tvir, cvel
            read(unitfile) real_table_dp(ihalo,23:24) ! rho0, rc
         end if
      end if

      return
   end subroutine read_halo_brick

!#########################################################
   subroutine read_halo_member(unitfile)
!#########################################################
      implicit none
      integer(kind=4) :: unitfile, npart_halo
      !real(kind=4)    :: aexp

      read(unitfile) npart_halo ! nparts
      read(unitfile) part_ids(ipart:ipart+npart_halo-1) ! ID of members
      ipart = ipart + npart_halo

      call skip_read(unitfile, 10)

      if(galaxy) then
         call skip_read(unitfile, 6)
      else
         call skip_read(unitfile, 2)
      end if

      return
   end subroutine read_halo_member

!#########################################################
   subroutine allocate_table(nhalo, nint, nreal, dp)
!#########################################################
      implicit none
      integer(kind=4) :: nhalo, nint, nreal
      logical :: dp

      call close()
      allocate(integer_table(1:nhalo, 1:nint))
      if(dp) then
         allocate(real_table_dp(1:nhalo, 1:nreal))
      else
         allocate(real_table(1:nhalo, 1:nreal))
      end if

   end subroutine allocate_table

!#########################################################
   subroutine close()
!#########################################################
      implicit none
      if(allocated(integer_table))deallocate(integer_table)
      if(allocated(real_table))deallocate(real_table)
      if(allocated(real_table_dp))deallocate(real_table_dp)
      if(allocated(part_ids))deallocate(part_ids)

   end subroutine close


!#########################################################
   subroutine read_tree_brick(halofile)
!#########################################################
      implicit none
      integer(kind=4)  :: unitfile,nhalo_snap,nmem=0,i
      real(kind=4)    :: aexp, massp, omega_t, age_univ
      integer(kind=4) :: nbodies, nb_of_halos, nb_of_subhalos!, nhalo
      character(len=128):: halofile

      unitfile = 55
      !write(*,*) '---> processing...  ', trim(halofile)
      open(unit=unitfile, file=halofile, form='unformatted', status='old')
      read(unitfile) nbodies
      read(unitfile) massp
      read(unitfile) aexp
      read(unitfile) omega_t
      read(unitfile) age_univ
      read(unitfile) nb_of_halos, nb_of_subhalos
      nhalo_snap = nb_of_subhalos + nb_of_halos
      !print *, '---> No. of halos (subhalos):', nb_of_halos, nb_of_subhalos

      do i=1, nhalo_snap
         call read_halo_brick(unitfile)!, aexp)
         ihalo = ihalo + 1
      enddo
      close(unitfile)

      !print *, '---> read tree_bricks: ok'
      return
   end subroutine read_tree_brick

!#########################################################
   subroutine read_member_brick(halofile)
!#########################################################
      implicit none
      integer(kind=4)  :: unitfile,ierr,nhalo_snap,nmem=0,i
      !real(kind=4)    :: aexp, massp, omega_t, age_univ
      integer(kind=4) :: nbodies, nb_of_halos, nb_of_subhalos!, nhalo
      character(len=128):: halofile

      unitfile = 55
      !write(*,*) '---> processing...  ', trim(halofile)
      open(unit=unitfile, file=halofile, form='unformatted', status='old')
      call skip_read(unitfile, 5)
      read(unitfile) nb_of_halos, nb_of_subhalos
      nhalo_snap = nb_of_subhalos + nb_of_halos

      do i=1, nhalo_snap
         call read_halo_member(unitfile)
      enddo
      close(unitfile)

      return
   end subroutine read_member_brick


!#########################################################
   subroutine read_bricks(repository,galaxy_ini,start,end,read_members,dp_ini)
!#########################################################
      implicit none
      integer(kind=4)::iout,nout, nb_of_halos, nb_of_subhalos
      character(len=5)::snout
      character(len=128)::halofile
      logical::ok_exist
      character(len=10)::iout_format

      character(len=128),intent(in)::repository
      logical,intent(in)::galaxy_ini,read_members
      integer(kind=4),intent(in):: start, end
      logical,intent(in):: dp_ini

      nhalo = 0
      nout = 0
      npart_tot = 0

      ihalo = 1
      ipart = 1

      galaxy = galaxy_ini
      dp = dp_ini

      iout_format='(I0.3)'

      do iout=start,end-1
         write(snout,TRIM(iout_format)) iout
         halofile = TRIM(repository)//snout
         inquire(file=halofile,exist=ok_exist)
         if(ok_exist)then
            open(unit=55, file=halofile, form='unformatted', status='old')
            call skip_read(55, 5)
            read(55) nb_of_halos, nb_of_subhalos
            nhalo = nhalo + nb_of_subhalos + nb_of_halos
            nout = nout + 1
            close(55)
         end if
      end do

      !write(*,*)'Number of bricks found:', nout
      !write(*,*)'Total number of halos:', nhalo

      if(galaxy)then
         call allocate_table(nhalo, 7, 27, dp)
      else
         call allocate_table(nhalo, 7, 23, dp)
      end if

      do iout=start,end-1
         write(snout,TRIM(iout_format)) iout
         halofile = TRIM(repository)//snout
         inquire(file=halofile,exist=ok_exist)
         if(ok_exist)then
            call read_tree_brick(halofile)
         end if
      end do


      if(read_members)then
         allocate(part_ids(1:npart_tot))
         do iout=start,end-1
            write(snout,TRIM(iout_format)) iout
            halofile = TRIM(repository)//snout
            inquire(file=halofile,exist=ok_exist)
            if(ok_exist)then
               call read_member_brick(halofile)
            end if
         end do
      end if

   end subroutine read_bricks

!#########################################################
   subroutine read_single_tree(treefile,galaxy_ini, dp_ini)
!#########################################################
      implicit none
      integer(kind=4) :: unitfile, nsteps, st, j
      character(len=10) :: iout_format
      logical :: ok_exist

      integer(kind=4), dimension(:), allocatable :: nb_of_halos, nb_of_subhalos
      real(kind=4),    dimension(:), allocatable :: aexp, omega_t, age_univ

      character(len=128),intent(in) :: treefile
      logical,intent(in) :: galaxy_ini, dp_ini

      ihalo = 1
      galaxy = galaxy_ini
      dp = dp_ini
      unitfile = 55

      inquire(file=treefile, exist=ok_exist)
      if(ok_exist)then
         open(unit=unitfile,file=treefile, form='unformatted', status='old')
         read(unitfile) nsteps

         write(*,*)'Number of steps found:', nsteps

         allocate(nb_of_halos(1:nsteps))
         allocate(nb_of_subhalos(1:nsteps))

         allocate(aexp(1:nsteps))
         allocate(omega_t(1:nsteps))
         allocate(age_univ(1:nsteps))

         read(unitfile) nb_of_halos, nb_of_subhalos

         nhalo = 0
         do st = 1, nsteps
            nhalo = nhalo + nb_of_halos(st) + nb_of_subhalos(st)
         end do

         write(*,*)'Total number of halos:', nhalo
         call allocate_table(nhalo, 18, 30, dp)

         read(unitfile) aexp
         read(unitfile) omega_t
         read(unitfile) age_univ

         write(6, '(a)', advance='no') 'Progress: '

         do st = 1, nsteps
            call progress_bar(st, nsteps)
            do j = 1, nb_of_halos(st)+nb_of_subhalos(st)
               call read_halo(unitfile, st, aexp(st), age_univ(st))
               ihalo = ihalo + 1
            end do
         end do
         close(unitfile)
      else
         write(*,*) "ERROR: Tree file not found."
      end if

   end subroutine read_single_tree


!#########################################################
   subroutine read_halo(unitfile, st, aexp, age_univ)
!#########################################################
      implicit none
      integer(kind=4) :: st, nb_fathers, nb_sons, temp_space, i, max_ind
      integer(kind=4) :: unitfile
      real(kind=4)    :: aexp, age_univ
      real(kind=8)    :: macc
      integer(kind=4), dimension(:), allocatable :: integer_temp
      real(kind=4),    dimension(:), allocatable :: real_temp

      read(unitfile) integer_table(ihalo, 1) ! my_number
      read(unitfile) !integer_table(ihalo, 2)                       ! BushID
      read(unitfile) !integer_table(ihalo, 2) ! my_timestep

      ! change the time step number to match the number of branching (time resolution)
      ! you decided your tree is going to have
      !integer_table(ihalo,3) = 0 ! meaningless in this program (TK)
      read(unitfile) integer_table(ihalo, 2),integer_table(ihalo, 3),integer_table(ihalo, 4),&
            integer_table(ihalo, 5),integer_table(ihalo, 6)
      ! level, hosthalo, hostsub, nbsub, nextsub

      !real_table(ihalo,1) = aexp ! aexp
      !real_table(ihalo,2) = age_univ ! aexp
      read(unitfile) real_table(ihalo, 1) ! m
      read(unitfile) macc ! macc (why f8?!)
      real_table(ihalo, 4) = real(macc, 4)
      read(unitfile) real_table(ihalo, 3:5) ! x
      read(unitfile) real_table(ihalo, 6:8) ! v
      read(unitfile) real_table(ihalo, 9:11) ! L
      read(unitfile) real_table(ihalo, 12:15) ! r, a, b, c
      read(unitfile) real_table(ihalo, 16:18) ! ek, ep, et
      read(unitfile) real_table(ihalo, 19) ! spin

      read(unitfile) nb_fathers
      integer_table(ihalo, 7) = nb_fathers
      integer_table(ihalo, 8:12) = -1
      real_table(ihalo, 20:24) = 0

      if(nb_fathers > 0) then
         temp_space = MAX(5, nb_fathers) ! maximum up to 5

         allocate(integer_temp(1:temp_space))
         allocate(real_temp(1:temp_space))
         integer_temp = -1
         real_temp = 0

         read(unitfile) integer_temp(1:nb_fathers) ! list_fathers
         read(unitfile) real_temp(1:nb_fathers) ! mass_fathers

         do i = 1, MIN(5, nb_fathers)
            ! father with higher contribution first
            max_ind = MAXLOC(real_temp, 1)
            integer_table(ihalo, 7+i) = integer_temp(max_ind)
            real_table(ihalo, 19+i) = real_temp(max_ind)
            real_temp(max_ind) = 0
         end do

  !         integer_table(ihalo, 9:13) = integer_temp(1:5)
  !         real_table(ihalo, 22:26) = real_temp(1:5)

         deallocate(integer_temp)
         deallocate(real_temp)
      end if

      read(unitfile) nb_sons
      integer_table(ihalo, 13) = nb_sons
      if(nb_sons > 0) then
         temp_space = MAX(5, nb_sons) ! maximum up to 5

         allocate(integer_temp(1:temp_space))
         integer_temp = -1

         read(unitfile) integer_temp(1:nb_sons) ! list_sons

         integer_table(ihalo, 14:18) = integer_temp(1:5)
         deallocate(integer_temp)
      else
         integer_table(ihalo, 14:18) = -1
      end if
      read(unitfile) real_table(ihalo,25:28) ! rvir, mvir, tvir, cvel
      read(unitfile) real_table(ihalo,29:30) ! rho0, rc
   end subroutine read_halo


!#####################################################################
   subroutine skip_read(unit,nskip)
!#####################################################################
      ! skip the given number of reads

      implicit none
      integer,intent(in) :: unit,nskip
      integer :: i
      do i=1,nskip
         read(unit)
      end do
   end subroutine skip_read

!#####################################################################
    subroutine progress_bar(iteration, maximum)
!#####################################################################
        implicit none
        integer :: iteration, maximum
        if(iteration == maximum) then
            write(6, '(a)') '> - Done!'
        elseif(MOD(iteration, MAX(1, maximum/50)) == 0) then
            write(6, '(a)', advance='no') '>'
        end if
    end subroutine progress_bar

end module readhtm
