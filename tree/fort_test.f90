program test_fort

!  use tmtree
  implicit none


integer::n_all_sons=30008625
character(52)::fn="/home/hoseung/Work/data/Horizon-AGN/halo/DM/tree.dat"
integer::n_halos_all=19185146
integer::nsteps=62
integer::n_all_fathers=47755013


integer, allocatable::fatherID(:), fatherIDx(:)
integer, allocatable::sonID(:)
real(KIND=4), allocatable ::fatherMass(:)
integer(KIND=4), allocatable ::i_arr(:,:)
real(KIND=4), allocatable ::f_arr(:,:)
real(KIND=4), allocatable ::aexp_arr(:), omega_t_arr(:), age_univ_arr(:)

allocate(fatherID(n_all_fathers+1))
allocate(fatherIDx(n_all_fathers+1))
allocate(sonID(n_all_sons))
allocate(fatherMass(n_all_fathers))
allocate(i_arr(n_halos_all,14))
allocate(f_arr(n_halos_all,25))
allocate(aexp_arr(nsteps))
allocate(omega_t_arr(nsteps))
allocate(age_univ_arr(nsteps))

call load_tree(fn, fatherID, fatherIDx, sonID, fatherMass, &
               & i_arr, f_arr, aexp_arr, omega_t_arr, age_univ_arr, &
               & n_halos_all, n_all_fathers, n_all_sons, .true., nsteps)

deallocate(fatherID, fatherIDx, sonID, i_arr, f_arr, aexp_arr, omega_t_arr, age_univ_arr)

end program 
