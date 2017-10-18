
# coding: utf-8
from fig1_funcs import *
# ## View details of odd galaxies.


# In[52]:
from analysis.cal_lambda import *

wdir = "../"

odd_list = np.array([160500811,  160500745, 3566300002, 2917600065, 2917600001,
                     1000200001, 3999000032, 3999000015, 3999000001, 3641300444,
                     3641500003])

clusterlist = np.array([np.int(odd/1e5) for odd in odd_list])

gals = []
for cluster in np.unique(clusterlist):
    i_this_cluster = clusterlist == cluster
    idgal_list = odd_list[i_this_cluster] - cluster*1e5
    print(cluster, "  ", idgal_list)
    gals.extend(do_main(wdir=wdir, w_wdir=wdir, w_galaxy_plot=True,
               w_galaxy_plot_dir= wdir +'galaxy_plot/',
               cluster = cluster,
               idgal_list=idgal_list))

#gals.extend(do_main(wdir=wdir, w_wdir=wdir, w_galaxy_plot=True,
#               w_galaxy_plot_dir= wdir +'galaxy_plot/',
#               idgal_list=np.array([1])))

#gals.extend(do_main(wdir=wdir, w_wdir=wdir, w_galaxy_plot=True,
#               w_galaxy_plot_dir= wdir +'galaxy_plot/',
#               idgal_list=np.array([444])))


plot_gals_all(gals, npix=60,
         fn_save = "./Odd_gal_plot_3")

