import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patheffects as path_effects
import sunpy
import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface
import astropy.units as u
import astropy.constants as const
from astropy.visualization import ImageNormalize
from astropy.coordinates import SkyCoord
from astropy.time import Time
from watroo import wow
from glob import glob
from tqdm import tqdm
import os
from sjireader import read_iris_sji
from itertools import chain

def plot_iris_map(eui_files, iris_1400_map_dir, iris_2796_map_dir,
                  eui_bottom_left, eui_top_right, 
                  iris_1400_example_map_region, iris_2796_example_map_region,
                  eis_vel_map, figsize=(12.6, 6.5), save_dir=None):
    
    eui_map = sunpy.map.Map(eui_files[0])
    
    eui_earth_time = Time(eui_map.meta['date_ear'])

    eui_map_fake = sunpy.map.Map(wow(eui_map.data, bilateral=1, denoise_coefficients=[5,5])[0],
                                 map_181.meta)
    
    eui_map_fake = eui_map_fake.submap(eui_bottom_left,top_right=eui_top_right)

    eui_map_boundary = get_map_edge_coords(eui_map_fake, step=10)
    
    iris_1400_map_ = read_iris_sji(iris_1400_map_dir, index=eui_earth_time, sdo_rsun=True)
    iris_1400_map_ = sunpy.map.Map(wow(iris_1400_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_1400_map_.meta).rotate()
    iris_2796_map_ = read_iris_sji(iris_2796_map_dir, index=eui_earth_time, sdo_rsun=True)
    iris_2796_map_ = sunpy.map.Map(wow(iris_2796_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_2796_map_.meta).rotate()

    with propagate_with_solar_surface():
        iris_1400_map_r = iris_1400_map_.reproject_to(iris_1400_example_map_region.wcs)
        iris_2796_map_r = iris_2796_map_.reproject_to(iris_2796_example_map_region.wcs)


    fig = plt.figure(figsize=figsize,layout='constrained')

    gs = fig.add_gridspec(1, 3, width_ratios=[348/501*722/452,1,1])

    ax1 = fig.add_subplot(gs[0],projection=eui_map_fake.wcs)
    ax2 = fig.add_subplot(gs[1],projection=iris_1400_example_map_region.wcs)
    ax3 = fig.add_subplot(gs[2],projection=iris_2796_example_map_region.wcs)
    

    im1 = ax1.imshow(eui_map_fake.data, cmap='sdoaia171', norm=ImageNormalize())
    im2 = ax2.imshow(iris_1400_map_r.data, cmap='irissji1400', norm=ImageNormalize())
    im3 = ax3.imshow(iris_2796_map_r.data, cmap='irissji2796', norm=ImageNormalize())

    for ax_ in (ax2,ax3):
        bounds = ax_.axis()
        with propagate_with_solar_surface():
            ax_.plot_coord(eui_map_boundary, color='white', lw=1, ls='--', alpha=0.8)
        ax_.axis(bounds)
        ax_.coords[0].axislabels.set_visible(False)
        ax_.coords[1].axislabels.set_visible(False)

    # ax2.coords[1].set_ticklabel_visible(False)
    ax3.coords[1].set_ticklabel_visible(False)

    ax1.set_xlabel('Solar-X (arcsec)')
    ax1.set_ylabel('Solar-Y (arcsec)')

    text1 = ax1.text(0.05, 0.97, r'$\rm{HRI_{EUV}}$' + '\n' + r'$t_\oplus$' + f'{eui_earth_time.isot[:-4]}',
                     ha='left', va='top', color='white', transform=ax1.transAxes,
                     path_effects=[path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()],
                     fontsize=12)

    text2 = ax2.text(0.05, 0.97, f'IRIS/SJI 140.0 nm\n{iris_1400_map_.date.isot[:-4]}', ha='left', va='top', color='white',
                transform=ax2.transAxes,
                path_effects=[path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()],
                fontsize=12)
    
    text3 = ax3.text(0.05, 0.97, f'IRIS/SJI 279.6 nm\n{iris_2796_map_.date.isot[:-4]}', ha='left', va='top', color='white',
                     transform=ax3.transAxes,
                     path_effects=[path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()],
                     fontsize=12)
    
    for ax_ in (ax1,ax2,ax3):
        ax_.coords.grid(color='white', alpha=0.6, lw=0.5, linestyle=':')

        with propagate_with_solar_surface(rotation_model='rigid'):
            bounds = ax_.axis()
            eis_vel_map.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
                                    axes=ax_)
            ax_.axis(bounds)
    
    fig.canvas.draw()
    # plt.show()

    def update_fig(ii, ims, texts, eui_files, iris_1400_map_dir, iris_2796_map_dir, 
                  iris_1400_example_map_region, iris_2796_example_map_region,):
        
        eui_map = sunpy.map.Map(eui_files[ii])
        
        eui_earth_time = Time(eui_map.meta['date_ear'])

        eui_map_fake = sunpy.map.Map(wow(eui_map.data, bilateral=1, denoise_coefficients=[3,3])[0],
                                    map_181.meta)
        
        eui_map_fake = eui_map_fake.submap(eui_bottom_left,top_right=eui_top_right)

        iris_1400_map_ = read_iris_sji(iris_1400_map_dir, index=eui_earth_time, sdo_rsun=True)
        iris_1400_map_ = sunpy.map.Map(wow(iris_1400_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_1400_map_.meta).rotate()
        iris_2796_map_ = read_iris_sji(iris_2796_map_dir, index=eui_earth_time, sdo_rsun=True)
        iris_2796_map_ = sunpy.map.Map(wow(iris_2796_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_2796_map_.meta).rotate()

        with propagate_with_solar_surface():
            iris_1400_map_r = iris_1400_map_.reproject_to(iris_1400_example_map_region.wcs)
            iris_2796_map_r = iris_2796_map_.reproject_to(iris_2796_example_map_region.wcs)

        ims[0].set_array(eui_map_fake.data)
        ims[1].set_array(iris_1400_map_r.data)
        ims[2].set_array(iris_2796_map_r.data)

        texts[0].set_text(r'$\rm{HRI_{EUV}}$' + '\n' + r'$t_\oplus$' + f'{eui_earth_time.isot[:-4]}')
        texts[1].set_text(f'IRIS/SJI 140.0 nm\n{iris_1400_map_.date.isot[:-4]}')
        texts[2].set_text(f'IRIS/SJI 279.6 nm\n{iris_2796_map_.date.isot[:-4]}')

        return ims, texts
    
    anim = FuncAnimation(fig, update_fig, frames = tqdm(range(len(eui_files))), #tqdm(range(10)), 
                            fargs=([im1,im2,im3], [text1,text2,text3], eui_files, iris_1400_map_dir, iris_2796_map_dir,
                                    iris_1400_example_map_region, iris_2796_example_map_region),
                            blit=False)
    
    if not os.path.exists(os.path.dirname(save_dir)):
        os.makedirs(os.path.dirname(save_dir))
    
    anim.save(save_dir, fps=30, bitrate=-1, codec="libx264", extra_args=['-pix_fmt', 'yuv420p',
                                                                         '-filter:v', 'crop=1260:650'])

        
        # iris_1400_map_ = sunpy.map.Map(wow(iris_1400_map[ii].data, bilateral=1, denoise_coefficients=[3,3])[0], iris_1400_map[ii].meta).rotate()
        # iris_2796_map_ = read_iris_sji(iris_2796_map_dir, index=iris_1400_map_.date, sdo_rsun=True)
        # iris_2796_map_ = sunpy.map.Map(wow(iris_2796_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_2796_map_.meta).rotate()

        # with propagate_with_solar_surface():
        #     iris_1400_map_r = iris_1400_map_.reproject_to(iris_1400_example_map_region.wcs)
        #     iris_2796_map_r = iris_2796_map_.reproject_to(iris_2796_example_map_region.wcs)

        # ims[0].images[0].set_array(iris_1400_map_r.data)
        # ims[1].images[0].set_array(iris_2796_map_r.data)

        # texts[0].set_text(f'IRIS/SJI 140.0 nm\n{iris_1400_map_.date.isot[:-4]}')
        # texts[1].set_text(f'IRIS/SJI 279.6 nm\n{iris_2796_map_.date.isot[:-4]}')

        # return ims, texts
    
    # anim = FuncAnimation(fig, update_fig, frames = tqdm(range(len(iris_1400_map))),
    #                         fargs=(fig.axes, [text1,text2], iris_1400_map, iris_2796_map_dir,
    #                                 iris_1400_example_map_region, iris_2796_example_map_region),
    #                         blit=False)
    
    # if not os.path.exists(os.path.dirname(save_dir)):
    #     os.makedirs(os.path.dirname(save_dir))
    
    # anim.save(save_dir, fps=30,dpi=150, bitrate=-1, codec="libx264", 
    #           extra_args=['-pix_fmt', 'yuv420p', '-filter:v', 'crop=1200:974'])
        

def get_map_edge_coords(map, step=1):
    map_edges = sunpy.map.map_edges(map)

    x_pix = []
    y_pix = []

    if map_edges[1].shape[0] % step != 0:
        iter_1 = chain(range(0, map_edges[1].shape[0], step), [map_edges[1].shape[0]-1])
    else:
        iter_1 = range(0, map_edges[1].shape[0], step)
    for ii in iter_1:
        x_pix.append(map_edges[1][ii,0].value)
        y_pix.append(map_edges[1][ii,1].value)

    if map_edges[3].shape[0] % step != 0:
        iter_3 = chain(range(0, map_edges[3].shape[0], step), [map_edges[3].shape[0]-1])
    else:
        iter_3 = range(0, map_edges[3].shape[0], step)

    for ii in iter_3:
        x_pix.append(map_edges[3][ii,0].value)
        y_pix.append(map_edges[3][ii,1].value)

    if map_edges[0].shape[0] % step != 0:
        iter_0 = chain(range(map_edges[0].shape[0]-1, -1, -step), [0])
    else:
        iter_0 = range(map_edges[0].shape[0]-1, -1, -step)

    for ii in iter_0:
        x_pix.append(map_edges[0][ii,0].value)
        y_pix.append(map_edges[0][ii,1].value)

    if map_edges[2].shape[0] % step != 0:
        iter_2 = chain(range(map_edges[2].shape[0]-1, -1, -step), [0])
    else:
        iter_2 = range(map_edges[2].shape[0]-1, -1, -step)

    for ii in iter_2:
        x_pix.append(map_edges[2][ii,0].value)
        y_pix.append(map_edges[2][ii,1].value)
    
    return map.pixel_to_world(x_pix*u.pix,y_pix*u.pix)


if __name__ == '__main__':

    eis_195_velmap = sunpy.map.Map("../../../src/EIS/DHB_007_v2/20221020T2343/sunpymaps/eis_195_velmap_shift.fits")

    eui_files = sorted(glob("../../../src/EUI/HRI/euv174/20221020/coalign_step_boxcar/*.fits"))
    # eui_map_seq_coalign = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)

    Txshift_hri, Tyshift_hri = (9.41462 - 20.8515)*u.arcsec, (7.05089-8.29747)*u.arcsec

    map_181 = sunpy.map.Map(eui_files[195]).shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0


    iris_sji_1400_map_example = read_iris_sji('../../../src/IRIS/20221020/1905/iris_l2_20221020_190518_3640007428_SJI_1400_t000.fits',
                                              index=Time('2022-10-24T19:21:41'), sdo_rsun=True).rotate()
    iris_sji_2796_map_example = read_iris_sji('../../../src/IRIS/20221020/1905/iris_l2_20221020_190518_3640007428_SJI_2796_t000_deconvolved.fits',
                                              index=Time('2022-10-24T19:21:41'), sdo_rsun=True).rotate()

    iris_sji_1400_map_example = iris_sji_1400_map_example.submap(SkyCoord(-900*u.arcsec, 180*u.arcsec, frame=iris_sji_1400_map_example.coordinate_frame),
                                                    top_right=SkyCoord(-825*u.arcsec, 300*u.arcsec, frame=iris_sji_1400_map_example.coordinate_frame))
    
    iris_sji_2796_map_example = iris_sji_2796_map_example.submap(SkyCoord(-900*u.arcsec, 180*u.arcsec, frame=iris_sji_2796_map_example.coordinate_frame),
                                                    top_right=SkyCoord(-825*u.arcsec, 300*u.arcsec, frame=iris_sji_2796_map_example.coordinate_frame))
    
        
    iris_sji_1400_map_dir = '../../../src/IRIS/20221020/1905/iris_l2_20221020_190518_3640007428_SJI_2796_t000_deconvolved.fits'
    iris_sji_2796_map_dir = '../../../src/IRIS/20221020/1905/iris_l2_20221020_190518_3640007428_SJI_2796_t000_deconvolved.fits'


    plot_iris_map(eui_files, iris_sji_1400_map_dir, iris_sji_2796_map_dir, 
                  [1700,150]*u.pix, [2047,650]*u.pix,
                  iris_sji_1400_map_example, iris_sji_2796_map_example, eis_195_velmap,
                  figsize=(12.6, 6.5), save_dir="../../../figs/ms_eis_eui_upflow_movie/iris_west_low_res.mp4")



