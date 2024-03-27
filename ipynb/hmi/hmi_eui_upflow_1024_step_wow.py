# This program uses sunkit_image 0.51, when the function calculate_match_template_shift 
# and mapsequence_coalign_by_match_template assumes the rotation matrix PC_ij is the identity matrix
# to convert the shift in arcsecs to pixels. The function only shift the image array, but don't update
# the WCS information. # At this point, this function only works for SDO cutout level 1.5 data.
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
import datetime
from watroo import wow
from tqdm import tqdm

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

def plot_colorbar(im, ax, width="3%", height="100%",loc="lower left",fontsize=10,
                  bbox_to_anchor=(1.02, 0., 1, 1),orientation="vertical"):
    clb_ax = inset_axes(ax,width=width,height=height,loc=loc,
                bbox_to_anchor=bbox_to_anchor,
                 bbox_transform=ax.transAxes,
                 borderpad=0)
    clb = plt.colorbar(im,pad = 0.05,orientation=orientation,ax=ax,cax=clb_ax)
    clb_ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    clb_ax.yaxis.get_offset_text().set_fontsize(fontsize)
    clb_ax.tick_params(labelsize=fontsize)
    return clb, clb_ax


def get_nearest_hmi_index(eui_map_date,hmi_date_obs, time_shift=datetime.timedelta(seconds=300)):
    return (np.abs(eui_map_date + time_shift - hmi_date_obs )).argmin()

def plot_eui_hmi_cutout(eui_map,hmi_map_hrifov,phi_los_map,eis_vel_map,bottom_left,top_right,
                        eui_norm,hmi_norm,eui_map_date, hmi_map_date,save_dir,figsize=(8,8),wow_filter=True):
    
    eui_map_cutout = eui_map.submap(bottom_left,top_right=top_right)
    hmi_map_cutout = hmi_map_hrifov.submap(bottom_left,top_right=top_right)
    phi_los_map_cutout = phi_los_map.submap(bottom_left,top_right=top_right)

    if wow_filter is True:
        eui_norm = ImageNormalize()
        eui_map_cutout = sunpy.map.Map(wow(eui_map_cutout.data,bilateral=1,weights=[],denoise_coefficients=[5, 2])[0],
                                       eui_map_cutout.meta)

    fig = plt.figure(figsize=figsize,layout='constrained')

    ax1 = fig.add_subplot(221,projection=eui_map_cutout)
    ax2 = fig.add_subplot(222,projection=hmi_map_cutout)
    ax3 = fig.add_subplot(223,projection=eui_map_cutout)
    ax4 = fig.add_subplot(224,projection=phi_los_map_cutout)

    im1 = eui_map_cutout.plot(axes=ax1,norm=eui_norm,cmap="solar orbiterhri_euv174",
                        title=f"EUI/HRI 17.4 nm {eui_map_date.strftime('%Y-%m-%d %H:%M:%S')}")
    im2 = hmi_map_cutout.plot(axes=ax2,norm=hmi_norm,cmap="hmimag",
                        title=f"HMI Blos {hmi_map_date.strftime('%Y-%m-%d %H:%M:%S')}")
    im3 = eui_map_cutout.plot(axes=ax3,norm=eui_norm,cmap="solar orbiterhri_euv174",
                        title=f"EUI/HRI 17.4 nm {eui_map_date.strftime('%Y-%m-%d %H:%M:%S')}")
    im4 = phi_los_map_cutout.plot(axes=ax4,norm=hmi_norm,cmap="hmimag",
                            title="PHI Blos 2022-10-24 19:15:03")
    
    clb2, clb_ax2 = plot_colorbar(im2,ax2,width="5%")
    clb4, clb_ax4 = plot_colorbar(im4,ax4,width="5%")

    cs = ax3.contourf(hmi_map_cutout.data,levels=[-1000,-10,10,1000],colors=["b","w","r"],alpha=[0.4,0,0.4])


    ax1.set_xlabel(" ")
    ax2.set_xlabel(" ")
    ax2.set_ylabel(" ")
    ax4.set_ylabel(" ")

    
    for ax_ in (ax1,ax2,ax3,ax4):
        bounds = ax_.axis()
        bounds = ax_.axis()
        eis_vel_map.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
                                axes=ax_)
        ax_.axis(bounds)

    # plt.show()
    if os.path.exists(os.path.dirname(save_dir)) is False:
        os.makedirs(os.path.dirname(save_dir))
    plt.savefig(save_dir,dpi=300,bbox_inches="tight")
    fig.clf()
    plt.close(fig)




if __name__ == '__main__':
    eui_files = sorted(glob("../../src/EUI/HRI/euv174/20221024/solo_L2_eui-hri*.fits"))
    eui_map_seq = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)

    eis_195_velmap_derot_repro_shifted_hrifov = sunpy.map.Map("../../src/EIS/DHB_007_v2/20221025T0023/sunpymaps/eis_195_velmap_derot_repro_hrifov.fits")
    eis_hhflare_195_velmap_derot_repro_hrifov = sunpy.map.Map("../../src/coalign_map/20221024/eis_hhflare_195_velmap_derot_repro_hrifov.fits")

    if os.path.exists ("../../src/EUI/HRI/euv174/20221024/coalign_shifts_step_mandal.h5"):
        with h5py.File("../../src/EUI/HRI/euv174/20221024/coalign_shifts_step_mandal.h5","r") as f:
            eui_map_seq_coalign_shifts_x = f["x"][()]
            eui_map_seq_coalign_shifts_y = f["y"][()]
        eui_map_seq_coalign_shifts = {"x":eui_map_seq_coalign_shifts_x*u.arcsec,"y":eui_map_seq_coalign_shifts_y*u.arcsec}
    # else:
    #     eui_map_exampe_shift_x = np.zeros(9)
    #     eui_map_exampe_shift_y = np.zeros(9)

    #     eui_map_example_seq = sunpy.map.Map(eui_files[21::40],sequence=True,memmap=True)

    #     for ii in range(3,-1,-1):
    #         eui_map_example_segment_ = eui_map_example_seq[ii:ii+2]
    #         eui_map_example_shift_ = coalignment.calculate_match_template_shift(eui_map_example_segment_,layer_index=-1)
    #         eui_map_exampe_shift_x[ii] = eui_map_example_shift_["x"].value[0] + eui_map_exampe_shift_x[ii+1]
    #         eui_map_exampe_shift_y[ii] = eui_map_example_shift_["y"].value[0] + eui_map_exampe_shift_y[ii+1]
        
    #     for ii in range(5,9):
    #         eui_map_example_segment_ = eui_map_example_seq[ii-1:ii+1]
    #         eui_map_example_shift_ = coalignment.calculate_match_template_shift(eui_map_example_segment_,layer_index=0)
    #         eui_map_exampe_shift_x[ii] = eui_map_example_shift_["x"].value[-1] + eui_map_exampe_shift_x[ii-1]
    #         eui_map_exampe_shift_y[ii] = eui_map_example_shift_["y"].value[-1] + eui_map_exampe_shift_y[ii-1]
    #     # eui_map_example_coalign = coalignment.mapsequence_coalign_by_match_template(eui_map_example_seq,shift=eui_map_example_shift)

    #     eui_map_seq_coalign_shifts_x = np.zeros(len(eui_files))
    #     eui_map_seq_coalign_shifts_y = np.zeros(len(eui_files))
        
    #     for ii in range(9):
    #         eui_map_seq_segment_ = eui_map_seq[ii*40:(ii+1)*40]
    #         eui_map_seq_coalign_shifts_ = coalignment.calculate_match_template_shift(eui_map_seq_segment_,layer_index=21)

    #         eui_map_seq_coalign_shifts_x[ii*40:(ii+1)*40] = eui_map_seq_coalign_shifts_["x"].value + eui_map_exampe_shift_x[ii]
    #         eui_map_seq_coalign_shifts_y[ii*40:(ii+1)*40] = eui_map_seq_coalign_shifts_["y"].value + eui_map_exampe_shift_y[ii]

    #     with h5py.File("../../src/EUI/HRI/euv174/20221024/coalign_shifts_step.h5","w") as f:
    #         f.create_dataset("x",data=eui_map_seq_coalign_shifts_x)
    #         f.create_dataset("y",data=eui_map_seq_coalign_shifts_y)

    #     eui_map_seq_coalign_shifts = {"x":eui_map_seq_coalign_shifts_x*u.arcsec,"y":eui_map_seq_coalign_shifts_y*u.arcsec}


    # eui_map_seq_coalign = coalignment.mapsequence_coalign_by_match_template(eui_map_seq,shift=eui_map_seq_coalign_shifts)
    # eui_map_seq_coalign.save("../../src/EUI/HRI/euv174/20221024/coalign_step/eui_map_seq_coalign_{index:03}.fits",overwrite=True)

    eui_files = sorted(glob("../../src/EUI/HRI/euv174/20221024/coalign_step/*.fits"))
    eui_map_seq_coalign = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)


    Txshift_hri, Tyshift_hri = (1.66986 + 2.49223)*u.arcsec,(7.60204 - 2.76366 )*u.arcsec

    map_181 = eui_map_seq_coalign[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0

    phi_los_map = sunpy.map.Map("../../src/coalign_map/20221024/phi_los_map_shifted.fits")
    phi_los_map_hrifov = phi_los_map.reproject_to(map_181.wcs)

    hmi_files = sorted(glob("../../src/HMI/20221024/lvl15_cutout/*.fits"))
    hmi_map_seq = sunpy.map.Map(hmi_files[:],sequence=True,memmap=True)
    hmi_date_obs = np.array([m.date.to_datetime() for m in hmi_map_seq])

    # hmi_map_repro_hrifov = []

    # for ii in range(len(hmi_map_seq)):
    #     with propagate_with_solar_surface(rotation_model="rigid"):
    #         hmi_map_repro_hrifov.append(hmi_map_seq[ii].reproject_to(map_181.wcs))

    # hmi_map_repro_hrifov = sunpy.map.Map(hmi_map_repro_hrifov,sequence=True)
    # hmi_map_repro_hrifov.save("../../src/HMI/20221024/lvl15_cutout_repro_hrifov/hmi_repro_hri_{index:03}.fits",overwrite=True)
    hmi_map_repro_hrifov_files = sorted(glob("../../src/HMI/20221024/lvl15_cutout_repro_hrifov/*.fits"))
    hmi_map_repro_hrifov = sunpy.map.Map(hmi_map_repro_hrifov_files,sequence=True,memmap=True)

    for ii, eui_map_ in enumerate(tqdm(eui_map_seq_coalign[:])):
        eui_map_date = eui_map_seq_coalign[ii].date.to_datetime()

        eui_map_fake = sunpy.map.Map(eui_map_.data,map_181.meta)

        hmi_map_index_ = get_nearest_hmi_index(eui_map_date,hmi_date_obs)

        hmi_map_ = hmi_map_repro_hrifov[hmi_map_index_]
        hmi_map_date = hmi_date_obs[hmi_map_index_]

        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [500,600]*u.pix,[670,760]*u.pix,
                            ImageNormalize(vmin=150,vmax=1.5e3,stretch=AsinhStretch(0.4)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi/east_1_wow/eui_hmi_cutout_{ii:03}.png")
        
        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [1800,200]*u.pix,[2000,500]*u.pix,
                            ImageNormalize(vmin=300,vmax=2.7e3,stretch=AsinhStretch(0.4)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi/west_2_wow/eui_hmi_cutout_{ii:03}.png")

        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [1700,500]*u.pix,[1970,800]*u.pix,
                            ImageNormalize(vmin=500,vmax=6e3,stretch=AsinhStretch(0.4)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi/west_1_wow/eui_hmi_cutout_{ii:03}.png")
        
        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_hhflare_195_velmap_derot_repro_hrifov,
                            [1120,740]*u.pix,[1260,900]*u.pix,
                            ImageNormalize(vmin=5e2,vmax=1.2e4,stretch=AsinhStretch(0.4)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi/center_1_wow/eui_hmi_cutout_{ii:03}.png")
        
        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [500,400]*u.pix,[710,600]*u.pix,
                            ImageNormalize(vmin=100,vmax=1.5e3,stretch=AsinhStretch(0.4)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi/east_2_wow/eui_hmi_cutout_{ii:03}.png")

        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [450,1100]*u.pix,[650,1400]*u.pix,
                            ImageNormalize(vmin=100,vmax=4e3,stretch=AsinhStretch(0.1)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi_noupflow/region_1_wow/eui_hmi_cutout_{ii:03}.png")

        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [500,800]*u.pix,[700,1000]*u.pix,
                            ImageNormalize(vmin=100,vmax=8e3,stretch=AsinhStretch(0.1)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi_noupflow/region_2_wow/eui_hmi_cutout_{ii:03}.png")
        
        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [850,700]*u.pix,[1050,1000]*u.pix,
                            ImageNormalize(vmin=100,vmax=7e3,stretch=AsinhStretch(0.1)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi_noupflow/region_3_wow/eui_hmi_cutout_{ii:03}.png")

        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [1100,1350]*u.pix,[1300,1600]*u.pix,
                            ImageNormalize(vmin=250,vmax=3e3,stretch=AsinhStretch(0.3)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi_noupflow/region_4_wow/eui_hmi_cutout_{ii:03}.png")

        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [1500,700]*u.pix,[1700,1000]*u.pix,
                            ImageNormalize(vmin=100,vmax=7e3,stretch=AsinhStretch(0.1)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi_noupflow/region_5_wow/eui_hmi_cutout_{ii:03}.png")

        plot_eui_hmi_cutout(eui_map_fake,hmi_map_,phi_los_map_hrifov,eis_195_velmap_derot_repro_shifted_hrifov,
                            [850,150]*u.pix,[1150,400]*u.pix,
                            ImageNormalize(vmin=150,vmax=2e3,stretch=AsinhStretch(0.15)),
                            ImageNormalize(vmin=-1000,vmax=1000),
                            eui_map_date, hmi_map_date,
                            f"../../figs/EUI/20221024/zoomin_hmi_noupflow/region_6_wow/eui_hmi_cutout_{ii:03}.png")

