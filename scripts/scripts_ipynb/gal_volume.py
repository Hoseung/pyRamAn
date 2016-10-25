
# coding: utf-8

# # Vispy First try
# 
# It works! 
# But the volume center is weird. 

# In[ ]:

# prepare a 3D volume data


# In[20]:

import load
import numpy as np
import tree.halomodule as hmo
import galaxy

wdir = './29172/'
nout = 180

info = load.info.Info(nout=nout, base=wdir)

gcat = hmo.Halo(base=wdir, is_gal=True, verbose=False, nout=nout)

gg = gcat.data[10]
galid = gg['id']

gm = load.rd_GM.rd_gal(nout, galid, wdir=wdir)
gm.cell = load.rd_GM.rd_cell(nout, galid, wdir=wdir)

gal = galaxy.galaxy.Galaxy(halo = gg, info=info)
good_gal = gal.mk_gal(gm.star, None, gm.cell, unit_conversion="GM", verbose=False)

celldata=gal.cell

points = np.stack((celldata['x'], celldata['y'], celldata['z'])).T


# In[29]:

xmin = -45
xmax = 45
npoints = 50


# In[30]:

# interpolate irregular data into uniform grid.
from scipy.interpolate import griddata
grid_x, grid_y, grid_z = np.mgrid[xmin:xmax:1j*npoints,
                                  xmin:xmax:1j*npoints,
                                  xmin:xmax:1j*npoints]
vol1 = griddata(points, celldata['metal'], (grid_x, grid_y, grid_z), method='nearest')


# In[32]:

from itertools import cycle

import numpy as np

from vispy import app, scene, io
from vispy.color import get_colormaps, BaseColormap

# Read volume
#vol1 = np.load(io.load_data_file('volume/stent.npz'))['arr_0']
#vol2 = np.load(io.load_data_file('brain/mri.npz'))['data']
#vol2 = np.flipud(np.rollaxis(vol2, 1))

# Prepare canvas
canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
canvas.measure_fps()

# Set up a viewbox to display the image with interactive pan/zoom
view = canvas.central_widget.add_view()

# Set whether we are emulating a 3D texture
emulate_texture = False

# Create the volume visuals, only one is visible
volume1 = scene.visuals.Volume(vol1, parent=view.scene, threshold=0.225,
                               emulate_texture=emulate_texture)

# move the center of camera rotation by(or to?) (0,0,0)
volume1.transform = scene.STTransform(translate=(0, 0, 0))
#volume2 = scene.visuals.Volume(vol2, parent=view.scene, threshold=0.2,
#                               emulate_texture=emulate_texture)
#volume2.visible = False

# Create two cameras (1 for firstperson, 3 for 3d person)
fov = 60.
cam1 = scene.cameras.FlyCamera(parent=view.scene, fov=fov)
cam2 = scene.cameras.TurntableCamera(parent=view.scene, fov=fov)
cam3 = scene.cameras.ArcballCamera(parent=view.scene, fov=fov)
view.camera = cam2  # Select turntable at first

# create colormaps that work well for translucent and additive volume rendering
class TransFire(BaseColormap):
    glsl_map = """
    vec4 translucent_fire(float t) {
        return vec4(pow(t, 0.5), t, t*t, max(0, t*1.05 - 0.05));
    }
    """


class TransGrays(BaseColormap):
    glsl_map = """
    vec4 translucent_grays(float t) {
        return vec4(t, t, t, t*0.05);
    }
    """

# Setup colormap iterators
opaque_cmaps = cycle(get_colormaps())
translucent_cmaps = cycle([TransFire(), TransGrays()])
opaque_cmap = next(opaque_cmaps)
translucent_cmap = next(translucent_cmaps)


# Implement key presses
@canvas.events.key_press.connect
def on_key_press(event):
    global opaque_cmap, translucent_cmap
    if event.text == '1':
        cam_toggle = {cam1: cam2, cam2: cam3, cam3: cam1}
        view.camera = cam_toggle.get(view.camera, 'fly')
    elif event.text == '2':
        methods = ['mip', 'translucent', 'iso', 'additive']
        method = methods[(methods.index(volume1.method) + 1) % 4]
        print("Volume render method: %s" % method)
        cmap = opaque_cmap if method in ['mip', 'iso'] else translucent_cmap
        volume1.method = method
        volume1.cmap = cmap
        #volume2.method = method
        #volume2.cmap = cmap
#    elif event.text == '3':
#        volume1.visible = not volume1.visible
#        volume2.visible = not volume1.visible
    elif event.text == '3':
        if volume1.method in ['mip', 'iso']:
            cmap = opaque_cmap = next(opaque_cmaps)
        else:
            cmap = translucent_cmap = next(translucent_cmaps)
        volume1.cmap = cmap
        #volume2.cmap = cmap
    elif event.text == '0':
        cam1.set_range()
        cam3.set_range()
    elif event.text != '' and event.text in '[]':
        s = -0.025 if event.text == '[' else 0.025
        volume1.threshold += s
        #volume2.threshold += s
        th = volume1.threshold# if volume1.visible else volume2.threshold
        print("Isosurface threshold: %0.3f" % th)


# for testing performance
#@canvas.connect
#def on_draw(ev):
    #canvas.update()

if __name__ == '__main__':
    print(__doc__)
    app.run()


# In[ ]:



