program test_hydro
    !se amr2cell_fun
    implicit none
    character*100   :: work_dir
    real(kind=8)    :: xmi, xma, ymi, yma, zmi, zma
    integer(kind=4) :: lmax, ngridtot, nvarh
    integer(kind=4), dimension(3) :: cpus


    work_dir="/home/hoseung/Work/data/NewHorizon/snapshots/output_00250"
    xmi=0.49
    xma=0.5
    ymi=0.49
    yma=0.5
    zmi=0.49
    zma=0.5
    lmax = 22
    cpus(1:3) =(/100,1100,2100 /)
    

    call a2c_count(ngridtot, nvarh, work_dir, xmi, xma, ymi, yma, zmi, zma, lmax, cpus)

end program
