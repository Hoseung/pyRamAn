import os
from sys import stdin

def edit_file(filename = 'twobody.pyf'):
    
    os.path.isfile(filename)
    if not os.path.isfile(filename):
        print('The file %s does not exist. Check the f2py compilation step. Returning...'%filename)
        return
    
    cmd = 'cp twobody.pyf twobody.pyf.bak'
    os.system(cmd)
    
    f1 = open(filename+'.bak', 'r')
    f2 = open(filename, 'w')
    for line in f1:
        if 'check' in line:
            arr = line.split(',')
            pos_i = 3
            for i in range(len(arr)):
                if 'check' in arr[i]: pos_i = i
            newline = arr[0]
            for i in range(1,len(arr)):
                if i!=pos_i:
                    if i==pos_i-1: newline = newline+',intent(hide)'
                    else: newline = newline+','+arr[i]
            f2.write(newline)
        else:
            f2.write(line)
    f1.close()
    f2.close()
    
    cmd = 'rm -rf twobody.pyf.bak'
    os.system(cmd)
    
    return

#########################################################################
edit_file()
