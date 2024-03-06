# This program uses sunkit_image 0.51, when the function calculate_match_template_shift 
# and mapsequence_coalign_by_match_template assumes the rotation matrix PC_ij is the identity matrix
# to convert the shift in arcsecs to pixels. 
# This is not the case for the EUI HRI images, but the function mapsequence_coalign_by_match_template
# converts the wrong shift in arcsecs to pixels in the wrong way as well, so the result is correct.
# The bug has been reported to the sunkit_image developers and might be fixed in the future.

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
import flow_vis
import cv2

eui_files = sorted(glob("../../src/EUI/HRI/euv174/20221024/solo_L2_eui-hri*.fits"))
eui_map_seq = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)
# eui_map_seq = coalignment.mapsequence_coalign_by_rotation(eui_map_seq,layer_index=181)
# eui_template = sunpy.map.Map(eui_files[181]).submap([500,500]*u.pix,top_right=[1500,1500]*u.pix)

eis_195_velmap_derot_repro_shifted_hrifov = sunpy.map.Map("../../src/EIS/DHB_007_v2/20221025T0023/sunpymaps/eis_195_velmap_derot_repro_hrifov.fits")
eis_hhflare_195_velmap_derot_repro_hrifov = sunpy.map.Map("../../src/coalign_map/20221024/eis_hhflare_195_velmap_derot_repro_hrifov.fits")

if os.path.exists ("../../src/EUI/HRI/euv174/20221024/coalign_shifts_step.h5"):
    with h5py.File("../../src/EUI/HRI/euv174/20221024/coalign_shifts_step.h5","r") as f:
        eui_map_seq_coalign_shifts_x = f["x"][()]
        eui_map_seq_coalign_shifts_y = f["y"][()]
    eui_map_seq_coalign_shifts = {"x":eui_map_seq_coalign_shifts_x*u.arcsec,"y":eui_map_seq_coalign_shifts_y*u.arcsec}
else:
    eui_map_exampe_shift_x = np.zeros(9)
    eui_map_exampe_shift_y = np.zeros(9)

    eui_map_example_seq = sunpy.map.Map(eui_files[21::40],sequence=True,memmap=True)

    for ii in range(3,-1,-1):
        eui_map_example_segment_ = eui_map_example_seq[ii:ii+2]
        eui_map_example_shift_ = coalignment.calculate_match_template_shift(eui_map_example_segment_,layer_index=-1)
        eui_map_exampe_shift_x[ii] = eui_map_example_shift_["x"].value[0] + eui_map_exampe_shift_x[ii+1]
        eui_map_exampe_shift_y[ii] = eui_map_example_shift_["y"].value[0] + eui_map_exampe_shift_y[ii+1]
    
    for ii in range(5,9):
        eui_map_example_segment_ = eui_map_example_seq[ii-1:ii+1]
        eui_map_example_shift_ = coalignment.calculate_match_template_shift(eui_map_example_segment_,layer_index=0)
        eui_map_exampe_shift_x[ii] = eui_map_example_shift_["x"].value[-1] + eui_map_exampe_shift_x[ii-1]
        eui_map_exampe_shift_y[ii] = eui_map_example_shift_["y"].value[-1] + eui_map_exampe_shift_y[ii-1]
    # eui_map_example_coalign = coalignment.mapsequence_coalign_by_match_template(eui_map_example_seq,shift=eui_map_example_shift)

    eui_map_seq_coalign_shifts_x = np.zeros(len(eui_files))
    eui_map_seq_coalign_shifts_y = np.zeros(len(eui_files))
    
    for ii in range(9):
        eui_map_seq_segment_ = eui_map_seq[ii*40:(ii+1)*40]
        eui_map_seq_coalign_shifts_ = coalignment.calculate_match_template_shift(eui_map_seq_segment_,layer_index=21)

        eui_map_seq_coalign_shifts_x[ii*40:(ii+1)*40] = eui_map_seq_coalign_shifts_["x"].value + eui_map_exampe_shift_x[ii]
        eui_map_seq_coalign_shifts_y[ii*40:(ii+1)*40] = eui_map_seq_coalign_shifts_["y"].value + eui_map_exampe_shift_y[ii]

    with h5py.File("../../src/EUI/HRI/euv174/20221024/coalign_shifts_step.h5","w") as f:
        f.create_dataset("x",data=eui_map_seq_coalign_shifts_x)
        f.create_dataset("y",data=eui_map_seq_coalign_shifts_y)

    eui_map_seq_coalign_shifts = {"x":eui_map_seq_coalign_shifts_x*u.arcsec,"y":eui_map_seq_coalign_shifts_y*u.arcsec}


eui_map_seq_coalign = coalignment.mapsequence_coalign_by_match_template(eui_map_seq,shift=eui_map_seq_coalign_shifts)

Txshift_hri, Tyshift_hri = (1.66986 + 2.49223)*u.arcsec,(7.60204 - 2.76366 )*u.arcsec

map_181 = eui_map_seq_coalign[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
map_181.meta["rsun_ref"] = 696000000.0
eui_map_region_east_181 = map_181.submap([400,200]*u.pix,top_right=[800,800]*u.pix)
eui_map_region_west_181 = map_181.submap([1500,0]*u.pix,top_right=[2048,800]*u.pix)
eui_map_region_center_181 = map_181.submap([1100,700]*u.pix,top_right=[1300,1000]*u.pix)

def draw_hsv_color_wheel(ax,vmin=0,vmax=1.5,samples=1024,clip_circle=True):
    xx, yy = np.meshgrid(
    np.linspace(-vmax, vmax, samples), np.linspace(-vmax, vmax, samples))

    rr = np.sqrt(xx**2 + yy**2) 

    color_wheel_rgb = flow_vis.flow_to_color(np.dstack((xx,yy)))

    if clip_circle == True:
        color_wheel_rgb[rr > vmax] = 255
        color_wheel_a = np.where(rr > vmax, 255, 255)
    else:
        color_wheel_a = np.ones(rr.shape)*255

    color_wheel_rgba = np.dstack((color_wheel_rgb,color_wheel_a))

    ax.imshow(color_wheel_rgba,extent=(-vmax,vmax,-vmax,vmax),origin="lower")

    ax.set_axis_off()


def plot_of_map(map_prev, map_next, map_date, eis_vel_map, save_dir, ii, norm, 
                bottom_left, top_right, params=[0.5, 3, 15, 3, 5, 1.2, 0], rebin=None):
    
    map_prev_cut = map_prev.submap(bottom_left, top_right=top_right)
    map_next_cut = map_next.submap(bottom_left, top_right=top_right)

    if rebin is not None:
        map_prev_cut = map_prev_cut.resample(u.Quantity(np.array(map_prev.data.shape)//2*u.pix))
        map_next_cut = map_next_cut.resample(u.Quantity(np.array(map_next.data.shape)//2*u.pix))

    prev = cv2.normalize(map_prev_cut.data, None, 0, 255, cv2.NORM_MINMAX).astype('uint8')
    next = cv2.normalize(map_next_cut.data, None, 0, 255, cv2.NORM_MINMAX).astype('uint8')
    flow = cv2.calcOpticalFlowFarneback(prev, next, None, *params)
    flow = flow - np.nanmedian(flow,axis=(0,1))

    flow_hsv_vis = flow_vis.flow_to_color(flow)

    fig = plt.figure(figsize=(10,10),layout="constrained")

    ax1 = fig.add_subplot(221,projection=map_prev_cut,)

    map_prev_cut.plot(axes=ax1,norm=norm,
                            title=f"EUI/HRI 17.4 nm {map_date}")

    ax2 = fig.add_subplot(222,projection=map_prev_cut,)
    map_prev_cut.plot(axes=ax2,norm=norm,
                            title=f"EUI/HRI 17.4 nm {map_date}")

    ax2.quiver(np.arange(0,flow.shape[1],5),np.arange(0,flow.shape[0],5),
                flow[::5,::5,0],flow[::5,::5,1],color='r')

    ax2.set_ylabel(' ')

    ax3 = fig.add_subplot(223,projection=map_prev_cut,)

    ax3.imshow(flow_hsv_vis,origin="lower")

    ax3_wheel_ax = inset_axes(ax3,width="100%",height="100%",loc="upper left",
                    bbox_to_anchor=[0,0.847,0.15,0.15],
                    bbox_transform=ax3.transAxes,
                    borderpad=0)

    draw_hsv_color_wheel(ax3_wheel_ax,clip_circle=True)

    ax3.grid("on",lw=0.5,ls=":",color="yellow",alpha=0.5)

    ax3.set_ylabel(ax1.get_ylabel())
    ax3.set_xlabel(ax1.get_xlabel())



    ax4 = fig.add_subplot(224,projection=map_prev_cut,)
    ax4.imshow((map_next_cut.data - map_prev_cut.data)/map_prev_cut.data,
               origin="lower",cmap="binary_r",vmin=-0.5,vmax=0.5)
    
    ax4.set_ylabel(' ')
    ax4.set_xlabel(ax1.get_xlabel())
    ax1.set_xlabel(' ')
    
    for ax_ in (ax1,ax2,ax3,ax4):
        bounds = ax_.axis()
        eis_vel_map.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
                                axes=ax_)
        ax_.axis(bounds)

    plt.savefig(fname=os.path.join(save_dir,f"eui_hri_20221024_{ii:03d}.png"),
                dpi=300)
    fig.clf()
    plt.close(fig)


for ii, map in enumerate(eui_map_seq_coalign[:-1]):
    map_prev = map.shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_next = eui_map_seq_coalign[ii+1].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_date = map.date
    # map = map.reproject_to(eui_map_seq_coalign[0].wcs)
    map_prev = sunpy.map.Map(map_prev.data,map_181.meta)
    map_next = sunpy.map.Map(map_next.data,map_181.meta)

    plot_of_map(map_prev, map_next, map_date, eis_195_velmap_derot_repro_shifted_hrifov,
                "../../figs/EUI/20221024/zoomin_optical_flow/east_1/",ii,
                ImageNormalize(vmin=150,vmax=1.5e3,stretch=AsinhStretch(0.4)),
                bottom_left=[500,600]*u.pix,top_right=[670,760]*u.pix)
    
    plot_of_map(map_prev, map_next, map_date, eis_195_velmap_derot_repro_shifted_hrifov,
                "../../figs/EUI/20221024/zoomin_optical_flow/west_2/",ii,
                ImageNormalize(vmin=300,vmax=2.7e3,stretch=AsinhStretch(0.4)),
                bottom_left=[1800,200]*u.pix,top_right=[2000,500]*u.pix)
    


    # eui_map_region_east = map.submap([400,200]*u.pix,top_right=[800,800]*u.pix)
    # eui_map_region_west = map.submap([1500,0]*u.pix,top_right=[2048,800]*u.pix)
    # eui_map_region_center = map.submap([1100,700]*u.pix,top_right=[1300,1000]*u.pix)

    # eui_map_region_east_zoomin_1 = eui_map_region_east.submap([100,400]*u.pix,
    #                                                         top_right=[270,560]*u.pix)
    # eui_map_region_east_zoomin_1.plot_settings["norm"] = ImageNormalize(vmin=150,vmax=1.5e3,stretch=AsinhStretch(0.4))

    # fig = plt.figure(figsize=(5,5),layout="constrained")
    # ax = fig.add_subplot(111,projection=eui_map_region_east_zoomin_1)

    # eui_map_region_east_zoomin_1.plot(axes=ax)

    # bounds = ax.axis()
    # eis_195_velmap_derot_repro_shifted_hrifov.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
    #                                                         axes=ax)
    # ax.axis(bounds)
    # ax.set_title(f"EUI/HRI {map_date}")

    # plt.savefig(fname=os.path.join("../../figs/EUI/20221024/zoomin/east_1/",f"eui_hri_20221024_{ii:03d}.png"),
    #             dpi=300)
    # fig.clf()
    # plt.close(fig)

    # eui_map_region_east_zoomin_2 = eui_map_region_east.submap([100,200]*u.pix,
    #                                                     top_right=[310,400]*u.pix)
    # eui_map_region_east_zoomin_2.plot_settings["norm"] = ImageNormalize(vmin=100,vmax=1.5e3,stretch=AsinhStretch(0.4))
    # # eui_map_region_east_zoomin_1.plot_settings["cmap"] = plt.get_cmap("sdoaia171").reversed()
    # fig = plt.figure(figsize=(5,5),layout="constrained")

    # ax = fig.add_subplot(111,projection=eui_map_region_east_zoomin_2)
    # eui_map_region_east_zoomin_2.plot(axes=ax)

    # bounds = ax.axis()
    # eis_195_velmap_derot_repro_shifted_hrifov.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
    #                                                         axes=ax)
    # ax.axis(bounds)
    # ax.set_title(f"EUI/HRI {map_date}")

    # plt.savefig(fname=os.path.join("../../figs/EUI/20221024/zoomin/east_2/",f"eui_hri_20221024_{ii:03d}.png"),
    #             dpi=300)
    # fig.clf()
    # plt.close(fig)

    # eui_map_region_center_zoomin_1 = eui_map_region_center.submap([20,40]*u.pix,
    #                                                     top_right=[160,200]*u.pix)
    # eui_map_region_center_zoomin_1.plot_settings["norm"] = ImageNormalize(vmin=5e2,vmax=1.2e4,stretch=AsinhStretch(0.4))

    # fig = plt.figure(figsize=(5,5),layout="constrained")
    # ax = fig.add_subplot(111,projection=eui_map_region_center_zoomin_1)
    # eui_map_region_center_zoomin_1.plot(axes=ax)

    # bounds = ax.axis()
    # eis_hhflare_195_velmap_derot_repro_hrifov.draw_contours(levels=[-10,-5,5,10]
    #             ,colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,axes=ax)
    # ax.axis(bounds)
    # ax.set_title(f"EUI/HRI {map_date}")

    # plt.savefig(fname=os.path.join("../../figs/EUI/20221024/zoomin/center_1/",f"eui_hri_20221024_{ii:03d}.png"),
    #             dpi=300)
    # fig.clf()
    # plt.close(fig)

    # eui_map_region_west_zoomin_1 = eui_map_region_west.submap([200,500]*u.pix,
    #                                                     top_right=[470,800]*u.pix)
    # eui_map_region_west_zoomin_1.plot_settings["norm"] = ImageNormalize(vmin=500,vmax=6e3,stretch=AsinhStretch(0.4))
    # fig = plt.figure(figsize=(5,5),layout="constrained")
    # ax = fig.add_subplot(111,projection=eui_map_region_west_zoomin_1)
    # eui_map_region_west_zoomin_1.plot(axes=ax)

    # bounds = ax.axis()
    # eis_195_velmap_derot_repro_shifted_hrifov.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
    #                                                         axes=ax)
    # ax.axis(bounds)
    # ax.set_title(f"EUI/HRI {map_date}")

    # plt.savefig(fname=os.path.join("../../figs/EUI/20221024/zoomin/west_1/",f"eui_hri_20221024_{ii:03d}.png"),
    #             dpi=300)
    # fig.clf()
    # plt.close(fig)

    # eui_map_region_west_zoomin_2 = eui_map_region_west.submap([300,200]*u.pix,
    #                                                     top_right=[500,500]*u.pix)
    # eui_map_region_west_zoomin_2.plot_settings["norm"] = ImageNormalize(vmin=300,vmax=2.7e3,stretch=AsinhStretch(0.4))
    # fig = plt.figure(figsize=(5,5),layout="constrained")
    # ax = fig.add_subplot(111,projection=eui_map_region_west_zoomin_2)
    # eui_map_region_west_zoomin_2.plot(axes=ax)

    # bounds = ax.axis()
    # eis_195_velmap_derot_repro_shifted_hrifov.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
    #                                                         axes=ax)
    # ax.axis(bounds)
    # ax.set_title(f"EUI/HRI {map_date}")

    # plt.savefig(fname=os.path.join("../../figs/EUI/20221024/zoomin/west_2/",f"eui_hri_20221024_{ii:03d}.png"),
    #             dpi=300)
    # fig.clf()
    # plt.close(fig)

    print(f"Done {ii}",flush=True)


