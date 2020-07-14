module readr
    implicit none
    ! Tables, used for communication with Pyhton
    ! cell: x(3), varh(5~7)
    ! part: x(3), v(3), m, (age, metal)
    ! sink: m, x(3), v(3), j(3), dMBHoverdt, dMEdoverdt, dMsmbh, d_avgptr, c_avgptr, v_avgptr, Esave,
    ! bhspin, spinmag, eps_sink, rho_star, rho_dm, low_star, low_dm, fast_star, fast_dm
    real(kind=8),    dimension(:,:),   allocatable::real_table

    ! cell: level, cpu
    ! part: id, level, cpu
    ! sink: id
    integer(kind=4), dimension(:,:),   allocatable::integer_table

    ! part: long id
    integer(kind=8), dimension(:,:),   allocatable::long_table

    ! part: family, tag
    integer(kind=1), dimension(:,:),   allocatable::byte_table

    integer :: nhvar

    ! Some useful header informations...
    real(kind=8) :: aexp

contains

!#####################################################################
    subroutine read_part(repo, iout, cpu_list, mode, xlims, verbose, longint)
!#####################################################################
        implicit none

        integer :: i, icpu, idim
        integer :: ncpu, ndim
        integer :: npart, nstar_int, nsink
        integer(kind=8) :: npart_tot, nstar, npart_c, ilong

        real(kind=8), dimension(:), INTENT(IN):: xlims

        integer :: part_n, nreal, nint, nbyte, nlong
        integer :: pint
        logical :: ok

        character(len=128) :: file_path

        character(len=128),    intent(in) :: repo
        integer,               intent(in) :: iout
        integer, dimension(:), intent(in) :: cpu_list
        character(len=10),     intent(in) :: mode
        logical,               intent(in) :: verbose
        logical,               intent(in) :: longint

        part_n = 30

        file_path = part_filename(repo, iout, cpu_list(1), mode)

        ! Step 1: Verify there is file
        inquire(file=file_path, exist=ok)
        if ( .not. ok ) then
            print *, file_path, ' not found.'
            stop
        endif

        ! Step 2: Count the total number of particles.
        open(unit=part_n, file=file_path, status='old', form='unformatted')
        read(part_n) ncpu
        read(part_n) ndim
        call skip_read(part_n, 2)
        if(longint) then
            read(part_n) nstar
        else
            read(part_n) nstar_int
            nstar = nstar_int
        end if
        call skip_read(part_n, 2)
        read(part_n) nsink
        close(part_n)

        npart_tot = 0
        do i = 1, SIZE(cpu_list)
            icpu = cpu_list(i)

            open(unit=part_n, file=part_filename(repo, iout, icpu, mode), status='old', form='unformatted')
            call skip_read(part_n, 2)
            read(part_n) npart
            close(part_n)
            npart_tot = npart_tot + npart
        end do

        ! Set coulum spaces of each datatype for different versions of RAMSES
        if(mode == 'nh' .or. mode == 'yzics') then ! New Horizon / YZiCS / Horizon-AGN
            if(nstar > 0 .or. nsink > 0) then
                nreal = 2*ndim + 3
            else
                nreal = 2*ndim + 1
            end if
            nint = 3
            nbyte = 0
        elseif(mode == 'iap' .or. mode == 'gem' .or. mode == 'fornax' .or. mode == 'none') then ! New RAMSES version that includes family, tag
            nreal = 2*ndim + 3
            nint = 3
            nbyte = 2
        end if

        if(longint) then
            nint = nint-1
            nlong = 1
        else
            nlong = 0
        end if

        call clear_all()

        ! Allocate space for particle data
        allocate(real_table(1:npart_tot, 1:nreal))
        allocate(integer_table(1:npart_tot, 1:nint))
        allocate(byte_table(1:npart_tot, 1:nbyte))
        if(longint)allocate(long_table(1:npart_tot, 1:nlong))

        ! Step 3: Read the actual particle data
        ! Current position for particle
        npart_c = 1

        if(verbose)write(6, '(a)', advance='no') 'Progress: '
        do i = 1, SIZE(cpu_list)
            if(verbose)call progress_bar(i, SIZE(cpu_list))
            icpu = cpu_list(i)

            open(unit=part_n, file=part_filename(repo, iout, icpu, mode), status='old', form='unformatted')
            ! Skip headers
            call skip_read(part_n, 2)
            read(part_n) npart
            call skip_read(part_n, 5)

            ! Read position(3), velocity(3), mass
            do idim = 1, 2*ndim+1
                read(part_n) real_table(npart_c:npart_c+npart-1, idim)
            end do

            ! Read id
            pint=1
            if(longint) then
                read(part_n) long_table(npart_c:npart_c+npart-1, 1)
            else
                read(part_n) integer_table(npart_c:npart_c+npart-1, pint)
                pint = pint+1
            end if
            ! Read level
            read(part_n) integer_table(npart_c:npart_c+npart-1, pint)
            pint = pint+1
            if(mode == 'nh' .or. mode == 'yzics') then
                ! If star or sink particles are activated, RAMSES adds epoch, metallicity information for particles.
                if(nstar > 0 .or. nsink > 0) then
                    read(part_n) real_table(npart_c:npart_c+npart-1, 2*ndim+2)
                    read(part_n) real_table(npart_c:npart_c+npart-1, 2*ndim+3)
                end if

                ! Add CPU information
                integer_table(npart_c:npart_c+npart-1, pint) = icpu

            elseif(mode == 'iap' .or. mode == 'gem' .or. mode == 'none' .or. mode == 'fornax') then
                ! family, tag
                read(part_n) byte_table(npart_c:npart_c+npart-1, 1)
                read(part_n) byte_table(npart_c:npart_c+npart-1, 2)


                ! If star or sink particles are activated, RAMSES adds epoch, metallicity information for particles.
                if(nstar > 0 .or. nsink > 0) then
                    read(part_n) real_table(npart_c:npart_c+npart-1, 2*ndim+2)
                    read(part_n) real_table(npart_c:npart_c+npart-1, 2*ndim+3)
                else
                    real_table(npart_c:npart_c+npart-1, 2*ndim+2:2*ndim+3) = 0d0
                end if

                ! Add CPU information
                integer_table(npart_c:npart_c+npart-1, pint) = icpu
            end if
            npart_c = npart_c + npart
            close(part_n)
        end do

        !!  region check. 
        do ilong=1,npart_c
            if (real_table(i,1)>=xlims(1).and.real_table(i,1)<=xlims(2).and. &
              & real_table(i,2)>=xlims(3).and.real_table(i,2)<=xlims(4).and. &
              & real_table(i,3)>=xlims(5).and.real_table(i,3)<=xlims(6)) then
                byte_table(i,1)=-128
              end if
        end do
    end subroutine read_part

!#####################################################################
    subroutine read_sinkprop(repo, iprop, drag_part, mode)
!#####################################################################
        implicit none

        integer :: nsink, ndim, i, ireal, iint
        integer :: sink_n, nreal, nint

        character(len=128),    intent(in) :: repo
        integer,               intent(in) :: iprop
        logical,               intent(in) :: drag_part
        character(len=10),     intent(in) :: mode

        sink_n = 40
        call clear_all()

        open(unit=sink_n, file=sinkprop_filename(repo, iprop), status='old', form='unformatted')
        read(sink_n) nsink
        read(sink_n) ndim
        read(sink_n) aexp
        call skip_read(sink_n, 3)

        if(drag_part) then
            nreal = 34
            nint = 3
        else
            nreal = 22
            nint = 1
        end if
        if(mode == 'fornax') nreal = nreal + 1
        allocate(real_table(1:nsink, 1:nreal))
        allocate(integer_table(1:nsink, 1:nint))

        iint = 1
        ireal = 1

        read(sink_n) integer_table(:,iint)
        iint = iint + 1
        do i=1,22
            read(sink_n) real_table(:, ireal)
            ireal = ireal + 1
        end do
        if(mode == 'fornax') then
            read(sink_n) real_table(:, ireal)
            ireal = ireal + 1
        end if
        if(drag_part) then
            do i=1,8
                read(sink_n) real_table(:, ireal)
                ireal = ireal + 1
            end do
            do i=1,2
                read(sink_n) integer_table(:, iint)
                iint = iint + 1
            end do
            do i=1,4
                read(sink_n) real_table(:, ireal)
                ireal = ireal + 1
            end do
        end if
        close(sink_n)

    end subroutine read_sinkprop

!#####################################################################
    subroutine read_sink(repo, iout, icpu, levelmin, nlevelmax)
!#####################################################################
        implicit none

        integer :: nsink, nindsink, ndim, i
        integer :: sink_n, nreal, nint, nstat

        character(len=128),    intent(in) :: repo
        integer,               intent(in) :: iout, icpu, levelmin, nlevelmax

        ndim = 3
        sink_n = 50
        call clear_all()

        open(unit=sink_n, file=sink_filename(repo, iout, icpu), form='unformatted')
        rewind(sink_n)
                
        read(sink_n) nsink
        read(sink_n) nindsink
        
        nstat = (ndim*2+1)*(nlevelmax-levelmin+1)
        nreal = 20+nstat
        nint = 1

        allocate(real_table(1:nsink, 1:nreal))
        allocate(integer_table(1:nsink, 1:nint))

        if(nsink>0) then
            read(sink_n) integer_table(1:nsink, 1) ! id
            do i=1,20
                read(sink_n) real_table(:, i) ! mass, pos*3, vel*3, t, dM*3, Esave, j*3, spin*3, spinmag, eps
            end do
            do i=21,21+nstat-1
                read(sink_n) real_table(:, i) ! stats
            end do
        endif
        close(sink_n)

    end subroutine read_sink

!#####################################################################
    subroutine clear_all()
!#####################################################################
        ! Clean old data table is used more then once
        implicit none
        if(allocated(real_table)) deallocate(real_table)
        if(allocated(integer_table)) deallocate(integer_table)
        if(allocated(byte_table)) deallocate(byte_table)
        if(allocated(long_table)) deallocate(long_table)
    end subroutine clear_all

!#####################################################################
    subroutine skip_read(unit,nskip)
!#####################################################################
        ! skip the given number of reads

        implicit none
        integer,intent(in) :: unit, nskip
        integer :: i
        do i=1,nskip
            read(unit)
        end do
    end subroutine skip_read

!#####################################################################
    character(len=5) function charind(iout)
!#####################################################################
        implicit none
        integer, intent(in) :: iout

        write(charind, '(I0.5)') iout

    end function charind

!#####################################################################
    character(len=128) function part_filename(repo, iout, icpu, mode)
!#####################################################################
        implicit none
        character(len=128), intent(in)  :: repo
        integer,            intent(in)  :: iout, icpu
        character(len=10),  intent(in)  :: mode

        if(mode == 'ng') then
            part_filename = TRIM(repo)//'/output_'//charind(iout)//'/part.out'//charind(icpu)
        else
            part_filename = TRIM(repo)//'/output_'//charind(iout)//'/part_'//charind(iout)//'.out'//charind(icpu)
        end if
    end function part_filename

!#####################################################################
    character(len=128) function sinkprop_filename(repo, iout)
!#####################################################################
        implicit none
        character(len=128), intent(in)  :: repo
        integer,            intent(in)  :: iout

        sinkprop_filename = TRIM(repo)//'/sink_'//charind(iout)//'.dat'
    end function sinkprop_filename

!#####################################################################
    character(len=128) function sink_filename(repo, iout, icpu)
!#####################################################################
        implicit none
        character(len=128), intent(in)  :: repo
        integer,            intent(in)  :: iout, icpu

        sink_filename = TRIM(repo)//'/output_'//charind(iout)//'/sink_'//charind(iout)//'.out'//charind(icpu)
    end function sink_filename

!#####################################################################
    subroutine progress_bar(iteration, maximum)
!#####################################################################
        implicit none
        integer :: iteration, maximum
        if(iteration == maximum) then
            write(6, '(a)') '# - Done!'
        elseif(MOD(iteration, MAX(1, maximum/50)) == 0) then
            write(6, '(a)', advance='no') '#'
        end if
    end subroutine progress_bar
end module readr

