def trajectory_multi(td,thisgal,save=None, out_dir='./'):
    import matplotlib.colors as colors
    import matplotlib as mpl
    import numpy as np
    import matplotlib.pyplot as plt

    color = td['aexp']/max(td['aexp'])

    marker_area = np.sqrt(td['rvir'])

    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(2,2,1)
    ax1.scatter(td['x'], td['y'], s = marker_area,facecolors='none', c=color)
    ax1.set_title('x-y')

    ax2 = fig.add_subplot(2,2,2)
    ax2.scatter(td['y'], td['z'], s = marker_area,facecolors='none', c=color)
    ax2.set_title('y-z')

    ax3=fig.add_subplot(2,2,3)
    cbar = ax3.scatter(td['x'], td['z'], s = marker_area,facecolors='none', c=color)
    ax3.set_title('x-z')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(cbar, cax = cbar_ax)

    ax4 = fig.add_subplot(2,2,4)
    ax4.axis('off')

    plt.text(0.5, 0.5,"Galaxy ID: {:5d}".format(thisgal) + " \n Mass {:.3e}".format(td['mvir'][0])  \
         , color='r',fontsize=20, fontname='Courier', \
         horizontalalignment='center', verticalalignment='center',\
         transform = ax4.transAxes)
	# if not explicitly deleted, figures generated by pyplot are retained!
	#

# add a colorbar to fig not axies.
    if save:
        fn_out = out_dir+str(thisgal)+'3D_multi.png'
        import os
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        plt.savefig(fn_out)
    else :
        plt.show()

    plt.close()
	# cla = clear current axis


def trajectory3D(td,thisgal,save=None):
    '''
    thisgal = name of halo. integer
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    
    x = td['x']
    y = td['y']
    z = td['z']

# Single 3D plot with projections
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = '3d')
# line plot
    N = len(x)
    for i in range(N-1):
        ax.plot(x[i:i+2], y[i:i+2], z[i:i+2], color=plt.cm.jet(255*i/N))

    ax.scatter(x,y,z, c = td['aexp']/max(td['aexp']))


# cx, cy, cz makes the data 2D-like.
# ... HOW???
    cx = np.ones_like(x) * ax.get_xlim3d()[0]
    cy = np.ones_like(x) * ax.get_ylim3d()[1]
    cz = np.ones_like(z) * ax.get_zlim3d()[0]

    for i in range(N-1):
        ax.plot(x[i:i+2], y[i:i+2], cz[i:i+2], color=plt.cm.jet(255*i/N))
        ax.plot(x[i:i+2], cy[i:i+2], z[i:i+2], color=plt.cm.jet(255*i/N))
        ax.plot(cx[i:i+2], y[i:i+2], z[i:i+2], color=plt.cm.jet(255*i/N))

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    if save is None:
        plt.show()
    else:
        plt.ioff()
        fn_out = save + str(thisgal)+'3d.png'
        plt.savefig(fn_out)
        plt.close()


def plot_all(tree, thisgal, save=False, out_dir='./', fn_save=None,
             nrows=None, ncols=None,
             quantities=None, normalizer=None ):
    """
    thisgal
    """
    import math
    import matplotlib.pyplot as plt
    import numpy as np
	#color = tree['aexp']/max(tree['aexp'])
	#marker_area = np.sqrt(tree['rvir'])

	# which quantities to plot?
	# 1,10,11,12,-16,18-

    if normalizer is None:
        normalizer = np.array([1e-11,1e-11,1,1,1,1,1,1,1e-11
        ,1,1e-11,1e-11,1e-11,1e-11	,1,1,1,1,1,1,1])

    if quantities is None:
        quantities=["sam_mvir","mvir","rvir","rs","vrms","vmax"
        ,"jx","jy","jz","spin","m200b","m200c","m500c","m2500c"
        ,"xoff","voff","btoc","ctoa","ax","ay","az"]

    units = ["$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$kpc/h$","$kpc/h$","$km/s$","$km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$M * L * V$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$kpc/h$","$km/s$","","","$kpc/h$","$kpc/h$","$kpc/h$"]

    '''
    units = ["$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$kpc/h$","$kpc/h$","$km/s$","$km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$kpc/h$","$km/s$","","","$kpc/h$","$kpc/h$","$kpc/h$"]
    '''

    # where returns ([i,n,d,i,c,e,s],).
    # whith the tailing ',', this is 2D array.

    if ncols is None:
        ncols = 6
    if nrows is None:
        nrows = math.ceil( len(quantities) / ncols)
    fig = plt.figure(figsize=(8,6))
    fig.suptitle("Galaxy ID: {:5d}".format(thisgal), fontsize=12)
    for i, q in enumerate(quantities):
        ax = fig.add_subplot(nrows, ncols, i)
        ax.plot(tree[q] * normalizer[i])
        ax.set_title(q)
        ax.locator_params(nbins = 2)
        ax.yaxis.set_label_coords(-0.3, 0.5)
        ax.set_ylabel(units[i])

	#plt.tight_layout()
	#plt.subplots_adjust(left=0.15, top=0.85) #global margin
    plt.subplots_adjust(wspace=0.7,hspace=0.5)

# if not explicitly deleted, figures generated by pyplot are retained!
# add a colorbar to fig not axies.
    if save:
        import os
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        if fn_save is None:
            fn_save = str(thisgal)+'halo_all.png'
            print("Default otput name is {0}".format(out_dir + fn_save))
        plt.savefig(out_dir + fn_save, bbox_inches='tight')
    else :
        plt.show()
    plt.close()


def plot_all_multiPDF(tree, thisgal, out_dir='./', fn_save=None,
             nrows=None, ncols=None,
             quantities=None, normalizer=None):
    import datetime
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    # Create the PdfPages object to which we will save the pages:

    if ncols is None:
        ncols = 3
    if nrows is None:
        nrows = 3

    if normalizer is None:
        normalizer = np.array([1e-11,1e-11,1,1,1,1,1,1,1e-11
        ,1,1e-11,1e-11,1e-11,1e-11	,1,1,1,1,1,1,1])

    if quantities is None:
        quantities=["sam_mvir","mvir","rvir","rs","vrms","vmax"
        ,"jx","jy","jz","spin","m200b","m200c","m500c","m2500c"
        ,"xoff","voff","btoc","ctoa","ax","ay","az"]

    units = ["$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$kpc/h$","$kpc/h$","$km/s$","$km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$M * L * V$"
    ,"$10^{11}M_{\odot}/h * Mpc/h * km/s$"
    ,"$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$10^{11}M_{\odot}/h$","$10^{11}M_{\odot}/h$"
    ,"$kpc/h$","$km/s$","","","$kpc/h$","$kpc/h$","$kpc/h$"]

#    with PdfPages(out_dir + fn_save) as pdf:
    pp = PdfPages(out_dir + fn_save)
    plt.ioff()
    for i, q in enumerate(quantities):
        if (i % (ncols*nrows)) == 0:
            fig = plt.figure(figsize=(11.69, 8.27), dpi=100)
            fig.suptitle("Galaxy ID: {:5d}".format(thisgal), fontsize=12)
        ax = fig.add_subplot(nrows, ncols, i % (ncols*nrows))
        ax.plot(tree[q] * normalizer[i])
        ax.set_title(q)
        ax.locator_params(nbins = 2)
        ax.yaxis.set_label_coords(-0.3, 0.5)
        ax.set_ylabel(units[i])
        if ((i + 1) % (ncols*nrows)) == 0 or (i + 1) == len(quantities):
            plt.subplots_adjust(wspace=0.7,hspace=0.5)
            pp.savefig(fig)
            plt.close()
#            plt.savefig(pp)
            #fig.suptitle("Galaxy ID: {:5d}".format(thisgal), fontsize=12)
            # We can also set the file's metadata via the PdfPages object:
    d = pp.infodict()
    d['Title'] = 'Multipage PDF Example'
    d['Author'] = u'Jouni K. Sepp\xe4nen'
    d['Subject'] = 'How to create a multipage pdf file and set its metadata'
    d['Keywords'] = 'PdfPages multipage keywords author title subject'
    d['CreationDate'] = datetime.datetime(2009, 11, 13)
    d['ModDate'] = datetime.datetime.today()
    pp.close()