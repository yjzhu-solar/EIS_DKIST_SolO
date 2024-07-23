import numpy as np
import sunpy 
import sunpy.map
from sunpy.coordinates import (propagate_with_solar_surface, 
                               Helioprojective, 
                               get_horizons_coord)
import eispac
from glob import glob
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.wcsapi import SlicedLowLevelWCS
from astropy.visualization import ImageNormalize, AsinhStretch
import astropy.constants as const
from astropy.io import fits
from sunraster.instr.spice import read_spice_l2_fits
from map_coalign import MapSequenceCoalign
import os
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib import animation

def interpolate_spice_map_to_target_wcs(spice_map, spice_coalign_wcs, spice_time, hri_map, target_wcs):
    spice_nx = spice_nt = spice_map.data.shape[1]
    spice_map = spice_map.submap([0, 120]*u.pix, top_right=[spice_nx, 699]*u.pix)
    spice_ny = spice_map.data.shape[0]
    spice_pix_t, spice_pix_y, spice_pix_x = np.indices((1,*spice_map.data.shape))
    spice_world_coords = spice_coalign_wcs.pixel_to_world(spice_pix_x, spice_pix_y, spice_pix_t)[0][0,:,:]

    solar_orbiter_loc = np.flip(get_horizons_coord('solar orbiter',
                                                {'start':spice_time[-1],
                                                'stop':spice_time[0],
                                                'step':f'{spice_nt}'}))
    
    spice_pix_y_in_target_wcs = np.zeros((spice_ny, spice_nx))
    spice_pix_x_in_target_wcs = np.zeros((spice_ny, spice_nx))

    for ii in range(spice_nt):
        spice_world_coord_t = SkyCoord(spice_world_coords[:,ii].Tx.to(u.arcsec), 
                                  spice_world_coords[:,ii].Ty.to(u.arcsec),
                                  frame='helioprojective',obstime=spice_time[ii], 
                                  observer=solar_orbiter_loc[ii], 
                                  rsun=hri_map.meta['rsun_ref']*u.m,)
        
        with propagate_with_solar_surface(rotation_model='rigid'):
            spice_pix_x_in_target_wcs[:,ii], spice_pix_y_in_target_wcs[:,ii] = target_wcs.world_to_pixel(spice_world_coord_t)

    hri_map_pix_y, hri_map_pix_x = np.indices(hri_map.data.shape)

    spice_map_interpolator = LinearNDInterpolator((spice_pix_x_in_target_wcs.flatten(), spice_pix_y_in_target_wcs.flatten()), spice_map.data.flatten())

    spice_map_interpolated = spice_map_interpolator(hri_map_pix_x, hri_map_pix_y)

    return sunpy.map.Map(spice_map_interpolated, hri_map.wcs)

def get_saffron_map(saffron_dir, saffron_filename, spice_cube_wcs, spice_time, hri_map, velmap=False):
    saffron_files = glob(os.path.join(saffron_dir, saffron_filename))

    saffron_intmaps = []
    if velmap:
        saffron_velmaps = []

    for saffron_file in saffron_files:
        saffron_map = sunpy.map.Map(saffron_file)
        saffron_intmap = saffron_map[0]
        saffron_intmap = interpolate_spice_map_to_target_wcs(saffron_intmap, spice_cube_wcs, spice_time, hri_map, hri_map.wcs)

        saffron_intmap.plot_settings['norm'] = ImageNormalize(vmin=np.nanpercentile(saffron_intmap.data, 0.5),
                                                                vmax=np.nanpercentile(saffron_intmap.data, 99.5),
                                                                stretch=AsinhStretch())
        saffron_intmaps.append(saffron_intmap)


        if velmap:
            saffron_velmap = saffron_map[1]

            saffron_velmap_data = saffron_velmap.data.copy()
            saffron_velmap_data = (saffron_velmap_data/np.nanmedian(saffron_velmap_data) - 1)*const.c.to_value(u.km/u.s)
            saffron_velmap_data = saffron_velmap_data - np.nanmedian(saffron_velmap_data[120:699,:], axis=0)
            saffron_velmap = sunpy.map.Map(saffron_velmap_data, saffron_map[1].meta)

            saffron_velmap = interpolate_spice_map_to_target_wcs(saffron_velmap, spice_cube_wcs, spice_time, hri_map, hri_map.wcs)
            saffron_velmap.plot_settings['norm'] = ImageNormalize(vmin=-40, vmax=40)
            saffron_velmaps.append(saffron_velmap)

    if velmap:
        if len(saffron_files) == 1:
            return saffron_intmaps[0], saffron_velmaps[0]
        else:
            return saffron_intmaps, saffron_velmaps
    else:
        if len(saffron_files) == 1:
            return saffron_intmaps[0]
        else:
            return saffron_intmaps

def plot_eui_cutout(eui_map_seq, spice_vel_map, bottom_left, top_right,  
                    eui_norm=None, save_dir=None, figsize=(4.5,4.5)):
    
    hri_map_seq_crop = eui_map_seq.submap(bottom_left, top_right=top_right)

    Txshift_hri, Tyshift_hri = (9.41462 - 20.8515)*u.arcsec, (7.05089-8.29747)*u.arcsec

    map_181 = hri_map_seq_crop[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0
    map_wcs = map_181.wcs

    fig = plt.figure(figsize=figsize,layout='constrained')
    ax = fig.add_subplot(111,projection=map_wcs)

    if eui_norm is None:
        eui_norm = ImageNormalize(vmin=np.nanpercentile(hri_map_seq_crop[181].data, 0.5),
                                vmax=np.nanpercentile(hri_map_seq_crop[181].data, 99.5),
                                stretch=AsinhStretch(0.15))

    im = hri_map_seq_crop[0].plot(axes=ax, norm=eui_norm)
    title = ax.set_title(f"EUI/HRI 17.4 nm {hri_map_seq_crop[0].date.strftime('%Y-%m-%d %H:%M:%S')}")

    bounds = ax.axis()
    spice_vel_map.draw_contours(levels=[-20],colors=["#005CAF",],alpha=0.8,
                                axes=ax)
    ax.axis(bounds)

    def update_plot(ii, im, title, hri_map_seq_crop):
        im.set_data(hri_map_seq_crop[ii].data)
        title.set_text(f"EUI/HRI 17.4 nm {hri_map_seq_crop[ii].date.strftime('%Y-%m-%d %H:%M:%S')}")

    anim = animation.FuncAnimation(fig, update_plot, frames=len(hri_map_seq_crop), 
                fargs=(im, title, hri_map_seq_crop), blit=False)

    if not os.path.exists(os.path.dirname(save_dir)):
        os.makedirs(os.path.dirname(save_dir))
    
    anim.save(save_dir, fps=30,dpi=400)



if __name__ == "__main__":
    eui_files = sorted(glob("../../../../src/EUI/HRI/euv174/20221020/coalign_step_boxcar/*.fits"))
    hri_map_seq = sunpy.map.Map(eui_files[:], memmap=True)
    hri_map_seq = MapSequenceCoalign(hri_map_seq)

    Txshift_hri, Tyshift_hri = (9.41462 - 20.8515)*u.arcsec, (7.05089-8.29747)*u.arcsec

    map_181 = hri_map_seq[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0

    spice_coalign_filename = '../../../../src/SPICE/20221020/solo_L2_spice-n-ras_20221020T231536_V06_150995364-000_coalign.fits'
    spice_coalign_cube = read_spice_l2_fits(spice_coalign_filename)

    spice_time = spice_coalign_cube['Ne VIII 770 - Peak'].time[0]

    with fits.open(spice_coalign_filename) as hdul:
        hdul.info()
        spice_coalign_header = hdul[3].header.copy()

    spice_coalign_header['CRVAL1'] = spice_coalign_header['CRVAL1'] - 14.2
    spice_coalign_header['CRVAL2'] = spice_coalign_header['CRVAL2'] - 4.5

    spice_coalign_wcs = WCS(spice_coalign_header).dropaxis(2)[:,120:700,:]

    saffron_dir = '../../../../src/SPICE/slimane/solo_L2.5_spice-n-ras_20221020T231536_V06_150995364-000/con-06'
    saffron_NeVIII_intmap, saffron_NeVIII_velmap = get_saffron_map(saffron_dir, '*770.42-ne_8*.fits',
                                                                spice_coalign_wcs, spice_time, map_181, velmap=True)
    
    # fig = plt.figure(figsize=(6,6),layout='constrained')
    # ax = fig.add_subplot(projection=map_181.wcs)
    # map_181.plot(axes=ax)
    # saffron_NeVIII_intmap.plot(axes=ax, alpha=0.5)
    # plt.show()
    
    plot_eui_cutout(hri_map_seq, saffron_NeVIII_velmap, [300, 700]*u.pix,
                    [500, 1000]*u.pix, 
                    save_dir='../../../../figs/EUI/20221020/upflows/20221020_EUI_SPICE_region1.mp4')
    

    


