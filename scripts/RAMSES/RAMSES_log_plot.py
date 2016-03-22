
# coding: utf-8

# In[3]:

#import matplotlib
#matplotlib.use("Agg") # This must be called before pyplot.


import numpy as np
import matplotlib.pyplot as plt

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


class Log():
    """
        Let's assume that there are at least some lines.
    """
    def __init__(self):
        self.M_st=[]
        self.M_econs=[]#None] * 1000
        self.M_epot=[]#None] * 1000
        self.M_ekin=[]#None] * 1000
        self.M_eint=[]#None] * 1000
        
        self.F_st=[]
        self.F_Error1=[]#None] * 10000
        self.F_Error2=[]#None] * 10000
        self.t=[]
        self.dt=[]
        self.a=[]
        self.mem_g=[]
        self.mem_p=[]
        self.L=[]
        self.ST=[]


# In[4]:

if __name__ == "__main__":
    wdir = input("work directory (including the tailing '/'): \n")
    filename = input("log file name: \n")#"C29172_h.o2516439"
    #nlines = file_len(filename)
    nlevelmax = 20
    ll = Log()
    with open(wdir + filename, 'r') as f:
        for line in f.readlines():
            #line = f.readline()        
            if "Main step" in line:
                ll.M_st.append(int(line[11:17]))
                ll.M_econs.append(float(line[24:34]))
                ll.M_epot.append(float(line[39:49]))
                ll.M_ekin.append(float(line[54:64])) 
                ll.M_eint.append(float(line[69:79]))               
                #print(int(line[11:17]), float(line[24:34]), float(line[39:49]), float(line[54:64]), float(line[69:79]))
            elif "Fine step" in line:
                ll.F_st.append(int(line[11:17]))    
                ll.t.append(float(line[20:33])) 
                ll.dt.append(float(line[36:47])) 
                ll.a.append(float(line[49:60])) 
                ll.mem_g.append(float(line[64:68])) 
                ll.mem_p.append(float(line[70:74])) 
                #print(int(line[11:17]), float(line[20:33]), float(line[36:47]), float(line[49:60]),
                #      float(line[64:68]), float(line[70:74]))
                #line = f.readline() # read next line
            elif "  ==> Level=" in line:
                ll.L.append(int(line[15:18])) 
                ll.ST.append(int(line[24:29])) 
                ll.F_Error1.append(float(line[37:47]))
                try:
                    ll.F_Error2.append(float(line[47:57])) 
                except:
                    ll.F_Error2.append(0.0)
                #print(int(line[15:18]), int(line[24:29]), float(line[37:47]))#, float(line[47:57]))

    from matplotlib.ticker import MaxNLocator

    fig, axs = plt.subplots(4,4)
    fig.set_size_inches(16,12)
    axs = axs.ravel()
    #plt.locator_params(axis = 'x', nbins = 4)
    #plt.locator_params(axis = 'y', nbins = 5)

    for i, name in enumerate(ll.__dict__.keys()):
        if name in ["F_Error1", "F_Error2", "dt"]:
            axs[i].plot(np.log10(getattr(ll, name)))
        else:
            axs[i].plot(getattr(ll, name))
        axs[i].set_title(name)
        axs[i].xaxis.set_major_locator(MaxNLocator(4))
        axs[i].yaxis.set_major_locator(MaxNLocator(5))

    plt.tight_layout()
    plt.suptitle(filename)
    plt.show()
    #plt.savefig(wdir + "log_plots.png", dpi=200)

