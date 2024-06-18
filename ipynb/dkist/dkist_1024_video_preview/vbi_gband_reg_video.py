import dkist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.visualization import ImageNormalize
from tqdm import tqdm
from scipy.io import readsav
from glob import glob

vbi_gband_dir = '../../../src/DKIST/vbi_1024/AEZDV/'
vbi_gband_dataset = dkist.load_dataset(vbi_gband_dir)
reg_save_files = sorted(glob('../../../sav/DKIST_reg/AEZDV/*.sav'))

fig, ax = plt.subplots(figsize=(10, 10),constrained_layout=True)

dkist_save = readsav(reg_save_files[0])

norm = ImageNormalize(vmin=np.nanpercentile(dkist_save['im_ref'], 0.1),
                        vmax=np.nanpercentile(dkist_save['im_ref'], 99.9))

im = ax.imshow(dkist_save['im_ref'], cmap='gray', norm=norm, origin='lower')
ax.set_title(f'VBI-B H-beta {vbi_gband_dataset[0].meta['headers'][0]['DATE-AVG']}')

def update_fig(ii, fig, ax, im, dataset):
    dkist_save = readsav(reg_save_files[ii])
    im.set_data(dkist_save['im_ref'])
    im.set_norm(ImageNormalize(vmin=np.nanpercentile(dkist_save['im_ref'], 0.1),
                                 vmax=np.nanpercentile(dkist_save['im_ref'], 99.9)))
    ax.set_title(f'VBI-B H-beta {dataset[ii].meta['headers'][0]['DATE-AVG']}')

anim = animation.FuncAnimation(fig, update_fig, frames=tqdm(range(vbi_gband_dataset.dimensions[0].value.astype(int)),), 
            fargs=(fig, ax, im, vbi_gband_dataset), blit=False)

anim.save('../../../figs/DKIST/VBI_preview/vbi_1024_gband_reg.mp4', fps=30,dpi=400)
