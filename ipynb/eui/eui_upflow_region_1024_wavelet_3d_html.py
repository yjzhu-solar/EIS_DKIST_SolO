import sunpy
import sunpy.map
from matplotlib import colormaps as cm
import astropy.units as u
from astropy.visualization import ImageNormalize, AsinhStretch
from map_coalign import MapSequenceCoalign
import numpy as np
import importlib
import image_wavelet
importlib.reload(image_wavelet)
from image_wavelet import WPSImage
from glob import glob
from copy import deepcopy
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import os

eui_files = sorted(glob("../../src/EUI/HRI/euv174/20221024/coalign_step_boxcar/*.fits"))
eui_map_seq_coalign = MapSequenceCoalign(sunpy.map.Map(eui_files[:]))


time = np.arange(360)*5.
scales = 5*np.logspace(1, 7,num=30, base=2)


def global_ws_test():

    x = np.arange(30,dtype=np.float64)
    y = np.arange(30,dtype=np.float64)

    xmesh, ymesh = np.meshgrid(x, y)

    frequency_mesh = np.sqrt(xmesh**2 + ymesh**2)/10.
    time = np.linspace(0,10,2000)
    data_sequence = np.sin(2*np.pi*frequency_mesh[:,:,np.newaxis]*time[np.newaxis,np.newaxis,:])

    scales = 1e-2*np.logspace(1, 7,num=50, base=2)

    wps = WPSImage(scales, data_sequence, time)
    wps._wps_image(ncpu="max")

    xmesh_3d, ymesh_3d, period_mesh_3d = np.meshgrid(x, y, wps.period)

    plt.plot(wps.period, wps.global_ws_unbias_coi[15,15,:])
    plt.show()

    fig = go.Figure()

    fig.add_trace(go.Volume(
        x=xmesh_3d.flatten(), y=ymesh_3d.flatten(), z=period_mesh_3d.flatten(),
        value=wps.global_ws_unbias_coi_sig.flatten(),
        isomin=100,
        opacity=0.1,
        surface_count=100,
        ))
    

    
    fig.write_html('../../figs/test_figs/wavelet_test.html',full_html=True)


def generate_3d_global_ws(bottom_left, top_right, filename,
                          eui_map_seq = eui_map_seq_coalign,
                          time=time, scales=scales, projection_type='orthographic'):
    
    eui_map_seq_cut = eui_map_seq.submap(bottom_left,top_right=top_right)
    eui_map_array = eui_map_seq_cut.as_array()

    wps = WPSImage(scales, eui_map_array, time,lag1=0.72)
    wps._wps_image(ncpu="max")

    ycoord, xcoord = np.arange(eui_map_array.shape[0]),  np.arange(eui_map_array.shape[1])
    xmesh, ymesh, period_mesh = np.meshgrid(xcoord, ycoord, wps.period)

    xmesh_2d, ymesh_2d = np.meshgrid(xcoord, ycoord)

    plt.plot(wps.period, wps.global_signif)
    plt.show()


    fig = go.Figure()

    fig.add_trace(go.Volume(
        x=xmesh.flatten(), y=ymesh.flatten(), z=period_mesh.flatten()/60.,
        value=wps.global_ws_unbias_coi.flatten(),
        opacity=0.1,
        surface_count=10,
        ))


    fig.add_trace(go.Surface(z=list(-4*np.ones_like(xmesh_2d)),
                            surfacecolor=eui_map_array[:,:,180],
                            #    cmin=100,cmax=2500,
                            colorscale='gray', showscale=False))

    fig.update_layout(scene_xaxis_showticklabels=True,
                    scene_yaxis_showticklabels=True,
                    scene=dict(
                        xaxis=dict(title="Solar-X [Pixels]"),
                        yaxis=dict(title="Solar-Y [Pixels]"),
                        zaxis=dict(title="Period [Minutes]"),
                        aspectmode='manual',
                        aspectratio=dict(x=1, y=1),
                        camera=dict(projection=dict(type=projection_type))
                    )
                    )

    # fig.update_layout(scene=dict(zaxis=dict(dtick=1, type='log')))

    if os.path.exists(os.path.dirname(filename)) == False:
        os.makedirs(os.path.dirname(filename))

    fig.write_html(filename,full_html=True)



generate_3d_global_ws([500,600]*u.pix, [670,760]*u.pix, 
                      "../../figs/wavelet_test/upflow_east_1.html")

generate_3d_global_ws([850,700]*u.pix,[1050,1000]*u.pix,
                        "../../figs/wavelet_test/noupflow_3.html")



