import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import pandas as pd
import sunpy
import sunpy.map
from sunpy.coordinates import (get_earth, get_horizons_coord,
                                Helioprojective, propagate_with_solar_surface)
import sunkit_image
import sunkit_image.coalignment as coalignment
import astropy
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import astropy.constants as const
from astropy.io import fits

import cmcrameri.cm as cmcm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (AutoLocator, AutoMinorLocator, 
    FixedLocator, FixedFormatter, LogLocator, StrMethodFormatter)
from ipywidgets import interactive, widgets
from IPython.display import display, clear_output
from astropy.visualization import (AsinhStretch, LinearStretch,
        LogStretch, ImageNormalize)
import os
from sun_blinker import SunBlinker
from copy import deepcopy   
from glob import glob
import h5py

eui_files = sorted(glob("../../src/EUI/HRI/euv174/20221024/solo_L2_eui-hri*.fits"))
eui_map_seq = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)
# eui_map_seq = coalignment.mapsequence_coalign_by_rotation(eui_map_seq,layer_index=181)
# eui_template = sunpy.map.Map(eui_files[181]).submap([500,500]*u.pix,top_right=[1500,1500]*u.pix)

eis_195_velmap_derot_repro_shifted_hrifov = sunpy.map.Map("../../src/EIS/DHB_007_v2/20221025T0023/sunpymaps/eis_195_velmap_derot_repro_hrifov.fits")
eis_hhflare_195_velmap_derot_repro_hrifov = sunpy.map.Map("../../src/coalign_map/20221024/eis_hhflare_195_velmap_derot_repro_hrifov.fits")

# if os.path.exists ("../../src/EUI/HRI/euv174/20221024/coalign_shifts.h5"):
#     with h5py.File("../../src/EUI/HRI/euv174/20221024/coalign_shifts.h5","r") as f:
#         eui_map_seq_coalign_shifts_x = f["x"][()]
#         eui_map_seq_coalign_shifts_y = f["y"][()]
#     eui_map_seq_coalign_shifts = {"x":eui_map_seq_coalign_shifts_x*u.arcsec,"y":eui_map_seq_coalign_shifts_y*u.arcsec}
# else:
#     eui_map_seq_coalign_shifts = coalignment.calculate_match_template_shift(eui_map_seq,layer_index=181)
#     eui_map_seq_coalign_shifts_x = eui_map_seq_coalign_shifts["x"].value
#     eui_map_seq_coalign_shifts_y = eui_map_seq_coalign_shifts["y"].value

#     with h5py.File("../../src/EUI/HRI/euv174/20221024/coalign_shifts.h5","w") as f:
#         f.create_dataset("x",data=eui_map_seq_coalign_shifts_x)
#         f.create_dataset("y",data=eui_map_seq_coalign_shifts_y)


eui_map_seq_coalign = eui_map_seq

Txshift_hri, Tyshift_hri = (1.67083 + 1.40322)*u.arcsec,(7.60192 - 2.32321 )*u.arcsec

map_181 = eui_map_seq_coalign[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
eui_map_region_east_181 = map_181.submap([400,200]*u.pix,top_right=[800,800]*u.pix)
eui_map_region_west_181 = map_181.submap([1500,0]*u.pix,top_right=[2048,800]*u.pix)
eui_map_region_center_181 = map_181.submap([1100,700]*u.pix,top_right=[1300,1000]*u.pix)

for ii, map in enumerate(eui_map_seq_coalign[:]):
    map = map.shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_date = map.date
    with propagate_with_solar_surface():
        map = map.reproject_to(map_181.wcs)
    eui_map_region_east = map.submap([400,200]*u.pix,top_right=[800,800]*u.pix)
    eui_map_region_west = map.submap([1500,0]*u.pix,top_right=[2048,800]*u.pix)
    eui_map_region_center = map.submap([1100,700]*u.pix,top_right=[1300,1000]*u.pix)

    fig = plt.figure(figsize=(10,6),constrained_layout=True)

    ax1 = fig.add_subplot(1,3,1,projection=eui_map_region_east_181)
    ax2 = fig.add_subplot(1,3,2,projection=eui_map_region_west_181)
    ax3 = fig.add_subplot(1,3,3,projection=eui_map_region_center_181)

    ax1.imshow(eui_map_region_east.data,origin="lower",cmap="solar orbiterhri_euv174",norm=ImageNormalize(vmin=10,vmax=4000,stretch=AsinhStretch(0.1)))
    ax2.imshow(eui_map_region_west.data,origin="lower",cmap="solar orbiterhri_euv174",norm=ImageNormalize(vmin=10,vmax=4000,stretch=AsinhStretch(0.1)))
    ax3.imshow(eui_map_region_center.data,origin="lower",cmap="solar orbiterhri_euv174",norm=ImageNormalize(vmin=10,vmax=1e4,stretch=AsinhStretch(0.3)))

    # eui_map_region_east.plot(axes=ax1,norm=ImageNormalize(vmin=10,vmax=4000,stretch=AsinhStretch(0.1)))
    # eui_map_region_west.plot(axes=ax2,norm=ImageNormalize(vmin=10,vmax=4000,stretch=AsinhStretch(0.1)))
    # eui_map_region_center.plot(axes=ax3,norm=ImageNormalize(vmin=10,vmax=4000,stretch=AsinhStretch(0.3)))

    ax2.set_ylabel(' ')
    ax3.set_ylabel(' ')
    ax1.set_ylabel('Solar-Y [arcsec]')

    ax1.set_title(f"East {map_date}")
    ax2.set_title(f"West {map_date}")
    ax3.set_title(f"Center {map_date}")

    for ax_ in (ax1,ax2):
        bounds = ax_.axis()
        eis_195_velmap_derot_repro_shifted_hrifov.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
                                                                   axes=ax_)
        ax_.axis(bounds)
    
    bounds = ax3.axis()
    eis_hhflare_195_velmap_derot_repro_hrifov.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
                                                             axes=ax3)
    ax3.axis(bounds)

    for ax_ in (ax1,ax2,ax3):
        ax_.grid("on",alpha=0.5,ls=":")
        ax_.set_xlabel('Solar-X [arcsec]')

    plt.savefig(fname=os.path.join("../../figs/EUI/20221024/upflow_video_nocoalign_wcs/",f"eui_hri_20221024_{ii:03d}.png"),
                dpi=300)
    fig.clf()
    plt.close(fig)


