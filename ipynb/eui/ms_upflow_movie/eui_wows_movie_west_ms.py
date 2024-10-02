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

def plot_eui_cutout(eui_map, bottom_left, top_right,
                    save_dir, figsize=(5,6)):
    
    eui_earth_time = Time(eui_map[0].meta['date_ear'])

    eui_map_fake = sunpy.map.Map(wow(eui_map[0].data, bilateral=1, denoise_coefficients=[5,5])[0],
                                 map_181.meta)
    
    eui_map_fake = eui_map_fake.submap(bottom_left,top_right=top_right)

    fig = plt.figure(figsize=figsize,layout='constrained')

    ax = fig.add_subplot(111,projection=eui_map_fake)
    im = eui_map_fake.plot(axes=ax, norm=ImageNormalize(), cmap='sdoaia171', title=None)
    ax.set_title(r"WOW $\rm{{HRI_{{EUV}}}}$ $t_{{\oplus}}$ {}".format(eui_earth_time.strftime("%Y-%m-%d %H:%M:%S")))

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
    
    anim.save(save_dir, fps=30, bitrate=-1, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'])

if __name__ == '__main__':

    # eis_195_velmap_derot_repro_shifted_hrifov = sunpy.map.Map("../../../src/EIS/DHB_007_v2/20221025T0023/sunpymaps/eis_195_velmap_derot_repro_hrifov.fits")
    # eis_hhflare_195_velmap_derot_repro_hrifov = sunpy.map.Map("../../../src/coalign_map/20221024/eis_hhflare_195_velmap_derot_repro_hrifov.fits")

    eui_files = sorted(glob("../../../src/EUI/HRI/euv174/20221026/coalign_step_boxcar/*.fits"))
    eui_map_seq_coalign = sunpy.map.Map(eui_files[:],sequence=True,memmap=True)

    Txshift_hri, Tyshift_hri = (-0.0235313 - 6.3736)*u.arcsec, (7.82867 - 0.685765)*u.arcsec

    map_181 = eui_map_seq_coalign[181].shift_reference_coord(Txshift_hri,Tyshift_hri)
    map_181.meta["rsun_ref"] = 696000000.0

    plot_eui_cutout(eui_map_seq_coalign,
                    [1700,300]*u.pix,
                    [2020,700]*u.pix,
                    "../../../figs/ms_eis_eui_upflow_movie/hri_west_low_res.mp4",)



