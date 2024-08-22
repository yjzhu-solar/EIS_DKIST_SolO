import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.contour import QuadContourSet
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
from astropy.time import Time
import dkist

import cmcrameri.cm as cmcm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (AutoLocator, AutoMinorLocator, 
    FixedLocator, FixedFormatter, LogLocator, StrMethodFormatter)
from ipywidgets import interactive, widgets
from IPython.display import display, clear_output
from astropy.visualization import (AsinhStretch, LinearStretch,
        LogStretch, ImageNormalize)
import os
from copy import deepcopy   
from glob import glob
import h5py
from tqdm import tqdm
from watroo import wow

def get_nearest_dkist_index(eui_map_date,dkist_date_obs, time_shift=300*u.s):
    return (np.abs(eui_map_date + time_shift - dkist_date_obs)).argmin()

def reproject_to_above_surface(map,target_wcs,radius=None,height=None):
    if radius is None and height is not None:
        rsun_ref = map.meta["rsun_ref"] + height.to_value(u.m)
    elif radius is not None and height is None:
        rsun_ref = radius.to_value(u.m)
    else:
        rsun_ref = map.meta["rsun_ref"]

    map_new = deepcopy(map)
    target_wcs_new = deepcopy(target_wcs)

    map_new.meta["rsun_ref"] = rsun_ref
    target_wcs_new.wcs.aux.rsun_ref = rsun_ref

    return map_new.reproject_to(target_wcs_new)

def plot_eui_hmi_cutout(eui_map,vbi_dataset,vbi_coalign_dir,vbi_date_obs,
                        vbi_target_wcs,eui_norm,save_dir,figsize=(8,8)):

    eui_map_date = eui_map[0].date
    eui_map_fake = sunpy.map.Map(eui_map[0].data,map_181.meta)
    eui_map_fake = eui_map_fake.submap([1100,500]*u.pix, 
                                       top_right=[1360,750]*u.pix)
    eui_map_fake = sunpy.map.Map(wow(eui_map_fake.data, bilateral=1, denoise_coefficients=[5,3])[0],
                                    eui_map_fake.meta)
    
    vbi_map_index_ = get_nearest_dkist_index(eui_map_date,vbi_date_obs)

    vbi_data_ = np.load(f'{vbi_coalign_dir}/vbi_hbeta_coalign_of_{vbi_map_index_:03d}.npz')['img']
    vbi_data_ = vbi_data_[32:-32,32:-32]
    vbi_date_ = vbi_date_obs[vbi_map_index_] 
    vbi_norm_ = ImageNormalize(vmin=np.nanpercentile(vbi_data_, 0.1),
                                vmax=np.nanpercentile(vbi_data_, 99.9))

    eui_map_cutout = reproject_to_above_surface(eui_map_fake,target_wcs=vbi_target_wcs,
                                                height=2.8*u.Mm)
    
    fig = plt.figure(figsize=figsize,layout='constrained')

    ax1 = fig.add_subplot(121,projection=eui_map_cutout.wcs)
    ax2 = fig.add_subplot(122,projection=vbi_target_wcs)


    im1 = eui_map_cutout.plot(axes=ax1,norm=eui_norm,cmap="solar orbiterhri_euv174")
    ax1.set_title(f"Reprojected EUI/HRI 17.4 nm {eui_map_date.strftime('%Y-%m-%d %H:%M:%S')}")
    
    im2 = ax2.imshow(vbi_data_,cmap="gray",origin="lower",norm=vbi_norm_)
    ax2.set_title(f"VBI-B H-Beta {vbi_date_.strftime('%Y-%m-%d %H:%M:%S')}")


    ax2.set_ylabel(" ")
    ax2.set_xlabel(" ")

    # plt.show()

    def update_fig(ii, fig, ims, axes, eui_map, vbi_dataset, vbi_date_obs):

        eui_map_date = eui_map[ii].date
        eui_map_fake = sunpy.map.Map(eui_map[ii].data,map_181.meta)
        eui_map_fake = eui_map_fake.submap([1100,500]*u.pix, 
                                        top_right=[1360,750]*u.pix)
        eui_map_fake = sunpy.map.Map(wow(eui_map_fake.data, bilateral=1, denoise_coefficients=[5,3])[0],
                                        eui_map_fake.meta)
    
        vbi_map_index_ = get_nearest_dkist_index(eui_map_date,vbi_date_obs)

        vbi_data_ = np.load(f'{vbi_coalign_dir}/vbi_hbeta_coalign_of_{vbi_map_index_:03d}.npz')['img']
        vbi_data_ = vbi_data_[32:-32,32:-32]
        vbi_date_ = vbi_date_obs[vbi_map_index_] 
        vbi_norm_ = ImageNormalize(vmin=np.nanpercentile(vbi_data_, 0.1),
                                    vmax=np.nanpercentile(vbi_data_, 99.9))

        eui_map_cutout = reproject_to_above_surface(eui_map_fake,target_wcs=vbi_target_wcs,
                                                    height=2.8*u.Mm)
        

        ims[0].set_data(eui_map_cutout.data)
        ims[1].set_data(vbi_data_)
        ims[1].set_norm(vbi_norm_)

        axes[0].set_title(f"Reprojected EUI/HRI 17.4 nm {eui_map_date.strftime('%Y-%m-%d %H:%M:%S')}")
        axes[1].set_title(f"VBI-B H-Beta {vbi_date_.strftime('%Y-%m-%d %H:%M:%S')}")

    anim = animation.FuncAnimation(fig, update_fig, frames=tqdm(range(len(eui_map))), # frames=tqdm(range(20)),
                fargs=(fig, [im1,im2], [ax1,ax2], eui_map, vbi_dataset, vbi_date_obs), blit=False)
    
    if not os.path.exists(os.path.dirname(save_dir)):
        os.makedirs(os.path.dirname(save_dir))
    
    anim.save(save_dir, fps=30,dpi=300)
    

if __name__ == '__main__':
    eui_files = sorted(glob("../../../src/EUI/HRI/euv174/20221024/coalign_step_boxcar/*.fits"))
    eui_map_seq_coalign = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)
    
    Txshift_hri, Tyshift_hri = (1.66986 + 2.49223)*u.arcsec,(7.60204 - 2.76366 )*u.arcsec

    map_181 = eui_map_seq_coalign[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0

    map_181_wcs = map_181.wcs

    vbi_hmi_xshift, vbi_hmi_yshift = 4.40*u.arcsec, 1.46*u.arcsec

    dkist_vbi_target_map = sunpy.map.Map('../../../src/DKIST/vbi_1024/BJOLO/VBI_2022_10_24T18_59_13_686_00486136_I_BJOLO_L1.fits')
    dkist_vbi_target_map = dkist_vbi_target_map.shift_reference_coord(vbi_hmi_xshift,vbi_hmi_yshift)
    vbi_target_wcs = dkist_vbi_target_map.wcs[257:-256:4,257:-256:4]

    vbi_dataset = dkist.load_dataset('../../../src/DKIST/vbi_1024/BJOLO/')
    vbi_coalign_dir = '../../../sav/DKIST_of/BJOLO/33'

    plot_eui_hmi_cutout(eui_map_seq_coalign,vbi_dataset,vbi_coalign_dir,Time(vbi_dataset.headers['DATE-AVG']),
                        vbi_target_wcs,ImageNormalize(),
                        save_dir='../../../figs/DKIST/VBI_preview/vbi_1024_hbeta_eui.mp4',figsize=(10,4.5))




