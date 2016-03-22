

dim_filter = file_lines('/home/inchani/Analysis/mock_img/sed/filter_JHK.dat')

dim_filter = dim_filter

filter_fn = fltarr(dim_filter,3) ; J H K
dim = fix((2480.0-1040.0)/10.) +1
lambda_arr = fltarr(dim)
filter_arr = dblarr(dim,3)


for i=0ul, dim-1 do begin
lambda_arr(i) = 1040.0 + i*10.
endfor


a=0.  & b=0. & j = -1 & start = -1

openr,lun,'/home/inchani/Analysis/mock_img/sed/filter_JHK.dat',/get_lun

for i=0ul, dim_filter-1 do begin
readf,lun,a,b

if b eq 0.000 then begin
start = start * (-1)
if start eq 1 then begin
j = j+1
print, j
ind = where( lambda_arr eq a)
k = ind(0)

endif
endif

bb = b * 1d0
bb = (bb * 1d-3)*1d3

filter_arr(k,j) = bb
k = k +1

endfor
free_lun,lun



openw,lun,'/home/inchani/Analysis/mock_img/sed/new_Filter_JHK.dat',/get_lun

for i=0ul, dim-1 do begin
printf,lun,lambda_arr(i)*1d1,filter_arr(i,0),filter_arr(i,1),filter_arr(i,2)
endfor

free_lun,lun




end
