import subprocess

nout = 368
snout=str(nout).zfill(3)
"""
fn_part = 'md5_'+snout+'part'+'.txt'
fn_hydro = 'md5_'+snout+'hydro'+'.txt'
fn_amr = 'md5_'+snout+'arm'+'.txt'
fn_grav = 'md5_'+snout+'grav'+'.txt'

fn_result_part = "result" + fn_part
fn_result_hydro = "result" + fn_hydro
fn_result_amr = "result" + fn_amr
fn_result_grav = "result" + fn_grav

subprocess.call(["md5sum -c", fn_part > fn_result_part) 
subprocess.call(["md5sum -c", fn_hydro > fn_result_hydro) 
subprocess.call(["md5sum -c", fn_amr > fn_result_amr) 
subprocess.call(["md5sum -c", fn_grav > fn_result_grav) 
"""

species=["part", "amr", "hydro", "grav"]
for i in range(4):
	fn = 'md5_'+snout+species[i]+'.txt'
	fn_result = "result" + fn
	subprocess.call(["md5sum -c", fn > fn_result) 



# Look for non-OK result 
"""
stderr.txt: OK
stdout.txt: FAILED
md5sum: WARNING: 1 computed checksum did NOT match
"""
