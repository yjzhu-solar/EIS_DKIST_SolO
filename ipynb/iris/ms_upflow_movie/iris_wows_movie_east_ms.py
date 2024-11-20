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
from astropy.time import Time
from watroo import wow
from glob import glob
from tqdm import tqdm
import os
from sjireader import read_iris_sji

def plot_iris_map(iris_1330_map, iris_1400_map_dir, iris_2796_map_dir, iris_1330_example_map_region_1, 
                  iris_1400_example_map_region_1, iris_2796_example_map_region_1,
                  iris_1330_example_map_region_2, 
                  iris_1400_example_map_region_2, iris_2796_example_map_region_2,
                  eis_vel_map, figsize=(8,6), save_dir=None):
    
    iris_1330_map_ = sunpy.map.Map(wow(iris_1330_map[0].data, bilateral=1, denoise_coefficients=[3,3])[0], iris_1330_map[0].meta).rotate()

    iris_1400_map_ = read_iris_sji(iris_1400_map_dir, index=iris_1330_map_.date, sdo_rsun=True)
    iris_1400_map_ = sunpy.map.Map(wow(iris_1400_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_1400_map_.meta).rotate()
    iris_2796_map_ = read_iris_sji(iris_2796_map_dir, index=iris_1330_map_.date, sdo_rsun=True)
    iris_2796_map_ = sunpy.map.Map(wow(iris_2796_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_2796_map_.meta).rotate()

    with propagate_with_solar_surface():
        iris_1330_map_r1 = iris_1330_map_.reproject_to(iris_1330_example_map_region_1.wcs)
        iris_1400_map_r1 = iris_1400_map_.reproject_to(iris_1400_example_map_region_1.wcs)
        iris_2796_map_r1 = iris_2796_map_.reproject_to(iris_2796_example_map_region_1.wcs)  

        iris_1330_map_r2 = iris_1330_map_.reproject_to(iris_1330_example_map_region_2.wcs)
        iris_1400_map_r2 = iris_1400_map_.reproject_to(iris_1400_example_map_region_2.wcs)
        iris_2796_map_r2 = iris_2796_map_.reproject_to(iris_2796_example_map_region_2.wcs)
    
    fig = plt.figure(figsize=figsize,layout='constrained')

    ax1 = fig.add_subplot(231, projection=iris_1330_example_map_region_1.wcs)
    ax2 = fig.add_subplot(232, projection=iris_1400_example_map_region_1.wcs)
    ax3 = fig.add_subplot(233, projection=iris_2796_example_map_region_1.wcs)
    ax4 = fig.add_subplot(234, projection=iris_1330_example_map_region_2.wcs)
    ax5 = fig.add_subplot(235, projection=iris_1400_example_map_region_2.wcs)
    ax6 = fig.add_subplot(236, projection=iris_2796_example_map_region_2.wcs)

    ax1.imshow(iris_1400_map_r1.data, cmap='irissji1400', norm=ImageNormalize())
    ax2.imshow(iris_1330_map_r1.data, cmap='irissji1330', norm=ImageNormalize())
    ax3.imshow(iris_2796_map_r1.data, cmap='irissji2796', norm=ImageNormalize())
    ax4.imshow(iris_1400_map_r2.data, cmap='irissji1400', norm=ImageNormalize())
    ax5.imshow(iris_1330_map_r2.data, cmap='irissji1330', norm=ImageNormalize())
    ax6.imshow(iris_2796_map_r2.data, cmap='irissji2796', norm=ImageNormalize())

    for ax_ in (ax1,ax2,ax3,ax5,ax6):
        ax_.coords[0].axislabels.set_visible(False)
        ax_.coords[1].axislabels.set_visible(False)

    for ax_ in (ax2,ax3,ax5,ax6):
        ax_.coords[1].set_ticklabel_visible(False)

    ax4.set_xlabel('Solar-X (arcsec)')
    ax4.set_ylabel('Solar-Y (arcsec)')

    text1 = ax1.text(0.05, 0.95, f'IRIS/SJI 140.0 nm\n{iris_1400_map_.date.isot[:-4]}', ha='left', va='top', color='white', 
             transform=ax1.transAxes,
             path_effects=[path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()],
             fontsize=12)
    
    text2 = ax2.text(0.05, 0.95, f'IRIS/SJI 133.0 nm\n{iris_1330_map_.date.isot[:-4]}', ha='left', va='top', color='white', 
             transform=ax2.transAxes,
             path_effects=[path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()],
             fontsize=12)
    
    text3 = ax3.text(0.05, 0.95, f'IRIS/SJI 279.6 nm\n{iris_2796_map_.date.isot[:-4]}', ha='left', va='top', color='white',
            transform=ax3.transAxes,
            path_effects=[path_effects.Stroke(linewidth=2, foreground='black'),path_effects.Normal()],
            fontsize=12)
    
    for ax_ in (ax1,ax2,ax3,ax4,ax5,ax6):
        ax_.coords.grid(color='white', alpha=0.6, lw=0.5, linestyle=':')

        with propagate_with_solar_surface(rotation_model='rigid'):
            bounds = ax_.axis()
            eis_vel_map.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
                                    axes=ax_)
            ax_.axis(bounds)

    fig.canvas.draw()
    # plt.show()

    def update_fig(ii, axes, texts, iris_1330_map, iris_1400_map_dir, iris_2796_map_dir, iris_1330_example_map_region_1, 
                  iris_1400_example_map_region_1, iris_2796_example_map_region_1,
                  iris_1330_example_map_region_2, 
                  iris_1400_example_map_region_2, iris_2796_example_map_region_2,):

        iris_1330_map_ = sunpy.map.Map(wow(iris_1330_map[ii].data, bilateral=1, denoise_coefficients=[3,3])[0], iris_1330_map[ii].meta).rotate()

        iris_1400_map_ = read_iris_sji(iris_1400_map_dir, index=iris_1330_map_.date, sdo_rsun=True)
        iris_1400_map_ = sunpy.map.Map(wow(iris_1400_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_1400_map_.meta).rotate()
        iris_2796_map_ = read_iris_sji(iris_2796_map_dir, index=iris_1330_map_.date, sdo_rsun=True)
        iris_2796_map_ = sunpy.map.Map(wow(iris_2796_map_.data, bilateral=1, denoise_coefficients=[3,3])[0], iris_2796_map_.meta).rotate()

        with propagate_with_solar_surface():
            iris_1330_map_r1 = iris_1330_map_.reproject_to(iris_1330_example_map_region_1.wcs)
            iris_1400_map_r1 = iris_1400_map_.reproject_to(iris_1400_example_map_region_1.wcs)
            iris_2796_map_r1 = iris_2796_map_.reproject_to(iris_2796_example_map_region_1.wcs)  

            iris_1330_map_r2 = iris_1330_map_.reproject_to(iris_1330_example_map_region_2.wcs)
            iris_1400_map_r2 = iris_1400_map_.reproject_to(iris_1400_example_map_region_2.wcs)
            iris_2796_map_r2 = iris_2796_map_.reproject_to(iris_2796_example_map_region_2.wcs)

        axes[0].images[0].set_array(iris_1400_map_r1.data)
        axes[1].images[0].set_array(iris_1330_map_r1.data)
        axes[2].images[0].set_array(iris_2796_map_r1.data)
        axes[3].images[0].set_array(iris_1400_map_r2.data)
        axes[4].images[0].set_array(iris_1330_map_r2.data)
        axes[5].images[0].set_array(iris_2796_map_r2.data)

        texts[0].set_text(f'IRIS/SJI 140.0 nm\n{iris_1400_map_.date.isot[:-4]}')
        texts[1].set_text(f'IRIS/SJI 133.0 nm\n{iris_1330_map_.date.isot[:-4]}')
        texts[2].set_text(f'IRIS/SJI 279.6 nm\n{iris_2796_map_.date.isot[:-4]}')

        return axes, texts
    
    anim = FuncAnimation(fig, update_fig, frames = tqdm(range(len(iris_1330_map))),  #tqdm(range(10))
                         fargs=(fig.axes, [text1,text2,text3], iris_1330_map, iris_1400_map_dir, iris_2796_map_dir, 
                                iris_1330_example_map_region_1, iris_1400_example_map_region_1, iris_2796_example_map_region_1,
                                iris_1330_example_map_region_2, iris_1400_example_map_region_2, iris_2796_example_map_region_2),
                         blit=False)

    if not os.path.exists(os.path.dirname(save_dir)):
        os.makedirs(os.path.dirname(save_dir))
    
    anim.save(save_dir, fps=5,dpi=150, bitrate=-1, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'])
        
        
        


# def plot_eui_cutout(eui_map, eis_vel_map, bottom_left, top_right,
#                     save_dir, figsize=(5,7)):
    
#     eui_earth_time = Time(eui_map[0].meta['date_ear'])

#     eui_map_fake = sunpy.map.Map(wow(eui_map[0].data, bilateral=1, denoise_coefficients=[3,3])[0],
#                                  map_181.meta)
    
#     eui_map_fake = eui_map_fake.submap(bottom_left,top_right=top_right)

#     fig = plt.figure(figsize=figsize,layout='constrained')

#     ax = fig.add_subplot(111,projection=eui_map_fake)
#     im = eui_map_fake.plot(axes=ax, norm=ImageNormalize(), cmap='sdoaia171', title=None)
#     ax.set_title(r"WOW $\rm{{HRI_{{EUV}}}}$ $t_{{\oplus}}$ {}".format(eui_earth_time.strftime("%Y-%m-%d %H:%M:%S")))

#     bounds = ax.axis()
#     bounds = ax.axis()
#     eis_vel_map.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
#                             axes=ax)
#     ax.axis(bounds)

#     fig.canvas.draw()
#     # plt.show()    

#     def update_fig(ii, im, ax):

#         eui_earth_time = Time(eui_map[ii].meta['date_ear'])
#         eui_map_fake = sunpy.map.Map(wow(eui_map[ii].data, bilateral=1, denoise_coefficients=[3,3])[0],
#                                      map_181.meta)
#         eui_map_fake = eui_map_fake.submap(bottom_left,top_right=top_right)
#         im.set_data(eui_map_fake.data)
#         ax.set_title(r"WOW $\rm{{HRI_{{EUV}}}}$ $t_{{\oplus}}$ {}".format(eui_earth_time.strftime("%Y-%m-%d %H:%M:%S")))
#         return im, ax
    
#     anim = FuncAnimation(fig, update_fig, frames = tqdm(range(len(eui_map))), #tqdm(range(10)),
#                          fargs=(im, ax), blit=False)

#     if not os.path.exists(os.path.dirname(save_dir)):
#         os.makedirs(os.path.dirname(save_dir))
    
#     anim.save(save_dir, fps=30,dpi=200, bitrate=-1, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'])

if __name__ == '__main__':

    eis_195_velmap = sunpy.map.Map("../../../src/EIS/DHB_007_v2/20221025T0023/sunpymaps/eis_195_velmap_derot.fits")

    eui_files = sorted(glob("../../../src/EUI/HRI/euv174/20221024/coalign_step_boxcar/*.fits"))
    eui_map_seq_coalign = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)

    Txshift_hri, Tyshift_hri = (1.66986 + 2.49223)*u.arcsec,(7.60204 - 2.76366 )*u.arcsec

    map_181 = eui_map_seq_coalign[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0
    hri_map_1024 = map_181.submap([300,200]*u.pix, top_right=[2000,1900]*u.pix)
    hri_map_1024_zoomin_1 = hri_map_1024.submap([210,400]*u.pix, top_right=[370,560]*u.pix)
    hri_map_1024_zoomin_2 = hri_map_1024.submap([230,200]*u.pix, top_right=[430,400]*u.pix)

    iris_sji_1330_map_example = read_iris_sji('../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_1330_t000.fits', 
                                              index=1, sdo_rsun=True).rotate()
    iris_sji_1400_map_example = read_iris_sji('../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_1400_t000.fits',
                                              index=1, sdo_rsun=True).rotate()
    iris_sji_2796_map_example = read_iris_sji('../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_2796_t000_deconvolved.fits',
                                              index=1, sdo_rsun=True).rotate()

    with propagate_with_solar_surface(rotation_model='rigid'):
        iris_sji_1330_map_example_region_1 = iris_sji_1330_map_example.submap(hri_map_1024_zoomin_1.bottom_left_coord,
                                                            top_right=hri_map_1024_zoomin_1.top_right_coord).rotate()
        iris_sji_1400_map_example_region_1 = iris_sji_1400_map_example.submap(hri_map_1024_zoomin_1.bottom_left_coord,
                                                            top_right=hri_map_1024_zoomin_1.top_right_coord).rotate()
        iris_sji_2796_map_example_region_1 = iris_sji_2796_map_example.submap(hri_map_1024_zoomin_1.bottom_left_coord,
                                                            top_right=hri_map_1024_zoomin_1.top_right_coord).rotate()
        
        iris_sji_1330_map_example_region_2 = iris_sji_1330_map_example.submap(hri_map_1024_zoomin_2.bottom_left_coord,
                                                            top_right=hri_map_1024_zoomin_2.top_right_coord).rotate()
        iris_sji_1400_map_example_region_2 = iris_sji_1400_map_example.submap(hri_map_1024_zoomin_2.bottom_left_coord,
                                                            top_right=hri_map_1024_zoomin_2.top_right_coord).rotate()
        iris_sji_2796_map_example_region_2 = iris_sji_2796_map_example.submap(hri_map_1024_zoomin_2.bottom_left_coord,
                                                            top_right=hri_map_1024_zoomin_2.top_right_coord).rotate()
    
        
    iris_sji_1330_map_all = read_iris_sji('../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_1330_t000.fits', sdo_rsun=True)
    iris_sji_1400_map_dir = '../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_1400_t000.fits'
    iris_sji_2796_map_dir = '../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_2796_t000_deconvolved.fits'
    # iris_sji_1400_map_all = read_iris_sji('../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_1400_t000.fits', sdo_rsun=True)
    # iris_sji_2796_map_all = read_iris_sji('../../../src/IRIS/20221024/2219/iris_l2_20221024_221954_3620511149_SJI_2796_t000_deconvolved.fits', sdo_rsun=True)

    plot_iris_map(iris_sji_1330_map_all, iris_sji_1400_map_dir, iris_sji_2796_map_dir,
                    iris_sji_1330_map_example_region_1, iris_sji_1400_map_example_region_1, iris_sji_2796_map_example_region_1,
                    iris_sji_1330_map_example_region_2, iris_sji_1400_map_example_region_2, iris_sji_2796_map_example_region_2,
                    eis_195_velmap,figsize=(8,6), save_dir="../../../figs/ms_eis_eui_upflow_movie/iris_east_low_res.mp4")



