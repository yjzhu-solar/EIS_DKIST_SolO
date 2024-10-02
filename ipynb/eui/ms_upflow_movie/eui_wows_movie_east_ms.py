import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sunpy
import sunpy.map
import astropy.units as u
import astropy.constants as const
from astropy.visualization import ImageNormalize
from astropy.time import Time
from watroo import wow
from glob import glob
from tqdm import tqdm
import os

def plot_eui_cutout(eui_map, eis_vel_map, bottom_left, top_right,
                    save_dir, figsize=(5,7)):
    
    eui_earth_time = Time(eui_map[0].meta['date_ear'])

    eui_map_fake = sunpy.map.Map(wow(eui_map[0].data, bilateral=1, denoise_coefficients=[3,3])[0],
                                 map_181.meta)
    
    eui_map_fake = eui_map_fake.submap(bottom_left,top_right=top_right)

    fig = plt.figure(figsize=figsize,layout='constrained')

    ax = fig.add_subplot(111,projection=eui_map_fake)
    im = eui_map_fake.plot(axes=ax, norm=ImageNormalize(), cmap='sdoaia171', title=None)
    ax.set_title(r"WOW $\rm{{HRI_{{EUV}}}}$ $t_{{\oplus}}$ {}".format(eui_earth_time.strftime("%Y-%m-%d %H:%M:%S")))

    bounds = ax.axis()
    bounds = ax.axis()
    eis_vel_map.draw_contours(levels=[-10,-5,5,10],colors=["#005CAF","#58B2DC","#F05E1C","#E83015"],alpha=0.8,
                            axes=ax)
    ax.axis(bounds)

    fig.canvas.draw()
    # plt.show()    

    def update_fig(ii, im, ax):

        eui_earth_time = Time(eui_map[ii].meta['date_ear'])
        eui_map_fake = sunpy.map.Map(wow(eui_map[ii].data, bilateral=1, denoise_coefficients=[3,3])[0],
                                     map_181.meta)
        eui_map_fake = eui_map_fake.submap(bottom_left,top_right=top_right)
        im.set_data(eui_map_fake.data)
        ax.set_title(r"WOW $\rm{{HRI_{{EUV}}}}$ $t_{{\oplus}}$ {}".format(eui_earth_time.strftime("%Y-%m-%d %H:%M:%S")))
        return im, ax
    
    anim = FuncAnimation(fig, update_fig, frames = tqdm(range(len(eui_map))), #tqdm(range(10)),
                         fargs=(im, ax), blit=False)

    if not os.path.exists(os.path.dirname(save_dir)):
        os.makedirs(os.path.dirname(save_dir))
    
    anim.save(save_dir, fps=30,dpi=200, bitrate=-1, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'])

if __name__ == '__main__':

    eis_195_velmap_derot_repro_shifted_hrifov = sunpy.map.Map("../../../src/EIS/DHB_007_v2/20221025T0023/sunpymaps/eis_195_velmap_derot_repro_hrifov.fits")
    eis_hhflare_195_velmap_derot_repro_hrifov = sunpy.map.Map("../../../src/coalign_map/20221024/eis_hhflare_195_velmap_derot_repro_hrifov.fits")

    eui_files = sorted(glob("../../../src/EUI/HRI/euv174/20221024/coalign_step_boxcar/*.fits"))
    eui_map_seq_coalign = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)

    Txshift_hri, Tyshift_hri = (1.66986 + 2.49223)*u.arcsec,(7.60204 - 2.76366 )*u.arcsec

    map_181 = eui_map_seq_coalign[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0

    plot_eui_cutout(eui_map_seq_coalign,
                    eis_195_velmap_derot_repro_shifted_hrifov,
                    [510,340]*u.pix,
                    [820,800]*u.pix,
                    "../../../figs/ms_eis_eui_upflow_movie/hri_east_low_res.mp4",)



