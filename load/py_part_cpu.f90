module part_load_interface
  use iso_c_binding, only: c_double, c_int, c_char
  use part_load_module, only:count_part, load_part
  implicit none
contains

subroutine c_count_part(ndm_actual, nstar_actual, nsink_actual, repository, xmin, xmax, ymin, ymax, zmin, zmax, cpu_list)
  integer(c_int),INTENT(OUT)::ndm_actual, nstar_actual, nsink_actual
  character(128, kind=c_char), INTENT(IN)::repository
  real(c_double), INTENT(IN)::xmin,xmax,ymin,ymax,zmin,zmax
  integer(c_int),dimension(:),INTENT(IN)::cpu_list
  call count_part(ndm_actual, nstar_actual, nsink_actual, repository, xmin, xmax, ymin, ymax, zmin, zmax, cpu_list)
end subroutine


  subroutine c_load_part(star_float, star_int, dm_float, dm_int, nstar_actual, &
             & ndm_actual, nsink_actual, repository, &
             & xmin, xmax, ymin, ymax, zmin, zmax, read_metal, cpu_list)

  real(c_double), INTENT(IN)::xmin,xmax,ymin,ymax,zmin,zmax
  integer(c_int),INTENT(IN)::ndm_actual, nstar_actual, nsink_actual
  integer(c_int),dimension(:),INTENT(IN)::cpu_list
  integer(c_int), INTENT(IN):: read_metal !# was logical
  character(128, kind=c_char), INTENT(IN)::repository

  real(c_double),INTENT(OUT),dimension(nstar_actual,9)::star_float
  real(c_double),INTENT(OUT),dimension(ndm_actual + nsink_actual,7)::dm_float
  integer(c_int),INTENT(OUT),dimension(nstar_actual)::star_int
  integer(c_int),INTENT(OUT),dimension(ndm_actual + nsink_actual)::dm_int

  call load_part(star_float, star_int, dm_float, dm_int, nstar_actual, &
                & ndm_actual, nsink_actual, repository, &
                & xmin, xmax,ymin, ymax, zmin, zmax, read_metal, cpu_list)

end subroutine

end module
