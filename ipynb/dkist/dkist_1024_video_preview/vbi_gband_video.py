import dkist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.visualization import ImageNormalize
from tqdm import tqdm

vbi_gband_dir = '../../../src/DKIST/vbi_1024/AEZDV/'
vbi_gband_dataset = dkist.load_dataset(vbi_gband_dir)

fig, ax = plt.subplots(figsize=(10, 10),constrained_layout=True)

norm = ImageNormalize(vmin=np.nanpercentile(vbi_gband_dataset[0].data, 0.1),
                        vmax=np.nanpercentile(vbi_gband_dataset[0].data, 99.9))

im = ax.imshow(vbi_gband_dataset[0].data, cmap='gray', norm=norm, origin='lower')
ax.set_title(f'VBI-B G-band {vbi_gband_dataset[0].meta['headers'][0]['DATE-AVG']}')


def update_fig(ii, fig, ax, im, dataset):
    im.set_data(dataset[ii].data)
    im.set_norm(ImageNormalize(vmin=np.nanpercentile(dataset[ii].data, 0.1),
                                 vmax=np.nanpercentile(dataset[ii].data, 99.9)))
    ax.set_title(f'VBI-B G-band {dataset[ii].meta['headers'][0]['DATE-AVG']}')

anim = animation.FuncAnimation(fig, update_fig, frames=tqdm(range(vbi_gband_dataset.dimensions[0].value.astype(int)),), 
            fargs=(fig, ax, im, vbi_gband_dataset), blit=False)

anim.save('../../../figs/DKIST/VBI_preview/vbi_1024_gband.mp4', fps=30,dpi=400)





