
fn_sed = 'bc2003_0500_sed.txt'


dim_age = 0ul & dim_sed =0ul & dim_lambda = 0ul

dim_age = file_lines('/home/inchani/Analysis/mock_img/sed/bc2003/original/ages_bc2003.dat')
dim_lambda = file_lines('/home/inchani/Analysis/mock_img/sed/bc2003/original/lambda_bc2003.dat')
dim_ref_lambda = file_lines('/home/inchani/Analysis/mock_img/sed/bc2003/original/ref_lambda_arr.dat')
dim_sed = dim_age * dim_lambda

check_sed = file_lines('/home/inchani/Analysis/mock_img/sed/bc2003/original/'+fn_sed)
if check_sed ne dim_sed then begin
print, 'error rrrrr'
exit
endif

sed_arr = dblarr(dim_lambda,dim_age)
lambda_arr = dblarr(dim_lambda)
ref_lambda_arr = dblarr(dim_ref_lambda)
age_arr = dblarr(dim_age)

a=0d0
openr,lun,'/home/inchani/Analysis/mock_img/sed/bc2003/original/ref_lambda_arr.dat',/get_lun
for i=0ul,dim_ref_lambda-1 do begin
	readf,lun,a
	ref_lambda_arr(i) =a
endfor
free_lun,lun

a=0d0
openr,lun,'/home/inchani/Analysis/mock_img/sed/bc2003/original/lambda_bc2003.dat',/get_lun
for i=0ul,dim_lambda-1 do begin
	readf,lun,a
	lambda_arr(i) =a
endfor
free_lun,lun


a=0d0
openr,lun,'/home/inchani/Analysis/mock_img/sed/bc2003/original/'+fn_sed,/get_lun
for i=0ul,dim_age-1 do begin
for j=0ul,dim_lambda-1 do begin
	readf,lun,a
	sed_arr(j,i) = a
endfor
endfor

openw,lun,'/home/inchani/Analysis/mock_img/sed/bc2003/'+fn_sed,/get_lun
for i=0ul,dim_age-1 do begin
ref = 1ul
for j=0ul,dim_lambda-1 do begin
	if (j mod 5) eq 0 then begin
	if ref_lambda_arr(ref-1) eq lambda_arr(j) then begin
		a=sed_arr(j,i)
		printf,lun,a
	endif else begin
		print,'error'
		exit
	endelse
	ref = ref + 1
	endif
endfor
if ref-1 ne dim_ref_lambda then begin
		print, 'warning!!!!'
		exit
endif

endfor
free_lun,lun









end
