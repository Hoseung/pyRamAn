import load
#work_dir = "/home/hoseung/Work/data/036370/snapshots/output_00080/"
work_dir = './output_00038'
# subroutine a2c_count(ngridtot, nvarh, repository, xmin, xmax, ymin, ymax, zmin, zmax, lmax)
#ng = 10
out= load.a2c.a2c_count(work_dir, 0.5, 0.6, 0.5, 0.6, 0.5, 0.6, 8)
print("done")
ngridtot = out[0]
nvarh=out[1]
print(ngridtot)
# subroutine a2c_load(xarr, dxarr, varr, repository, xmin, xmax, ymin, ymax, zmin, zmax, lmax, ngridtot, nvarh)
cell = load.a2c.a2c_load(work_dir, 0.5, 0.6, 0.5, 0.6, 0.5, 0.6, 8, ngridtot, nvarh)



