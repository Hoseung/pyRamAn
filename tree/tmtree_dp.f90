subroutine count_tree(fn, n_halos_all,  flist_index, slist_index, nsteps, big_run)
    implicit none
    character(LEN=256), INTENT(IN)::fn

    integer::nfathers, nsons
    integer, dimension(:), allocatable::nb_of_halos, nb_of_subhalos
    integer, INTENT(OUT)::nsteps, n_halos_all
    integer(KIND=8), INTENT(OUT)::flist_index, slist_index

    integer(KIND=4)::i,j, nhals_now
    logical, INTENT(IN)::big_run

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
            flist_index = flist_index + nfathers

            read(1)!Father ID
            read(1)!Father Mass
            read(1)nsons
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

        write(*,*)flist_index
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
    integer, INTENT(IN)::n_all_sons, n_all_fathers, nsteps
    integer(KIND=8), INTENT(IN)::n_halos_all

    integer(KIND=4):: nsons
    integer(KIND=8):: flist_index, slist_index
    integer(KIND=4):: i,j, k,nhals_now, n_fathers, idx_old, nhals_old, idx
    real(KIND=8)::macc

    integer, dimension(:), allocatable::fid_tmp, nb_of_halos, nb_of_subhalos, sonID_tmp

    ! +1 so that in python the array content start from index 1.
    integer, dimension(1:n_all_fathers+1), INTENT(OUT)::fatherID, fatherIDx
    integer, dimension(1:n_all_sons), INTENT(OUT)::sonID
    real(KIND=8), dimension(1:n_all_fathers), INTENT(OUT) ::fatherMass
    integer(KIND=4), dimension(1:n_halos_all,1:15), INTENT(OUT) ::i_arr
    real(KIND=8), dimension(1:n_halos_all,1:25), INTENT(OUT) ::f_arr
    real(KIND=8), dimension(1:nsteps), INTENT(OUT) ::aexp_arr, omega_t_arr, age_univ_arr

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
                    sonID(slist_index:slist_index+nsons-1)=sonID_tmp(1:nsons-1)+idx_old+nhals_now+nhals_old-1
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
