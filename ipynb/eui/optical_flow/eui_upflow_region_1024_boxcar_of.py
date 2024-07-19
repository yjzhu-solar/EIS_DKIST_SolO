import numpy as np
import sunpy 
import sunpy.map
from map_coalign import MapSequenceCoalign
from watroo import AtrousTransform, Triangle,denoise
import cv2
from glob import glob
import astropy.units as u
import h5py
from tqdm import tqdm


def calc_boxcar_of(eui_map_seq, bottom_left=None, top_right=None,
                   denoise_sigma=[5,3], of_params=[0.5, 3, 15, 3, 5, 1.2, 0]):
    eui_map_seq_crop = eui_map_seq.submap(bottom_left, top_right=top_right)

    array_of = np.zeros([*eui_map_seq_crop[0].data.shape,2,len(eui_map_seq_crop)])

    img_prev = denoise(eui_map_seq_crop[0].data, denoise_sigma, Triangle)

    for ii, map in enumerate(tqdm(eui_map_seq_crop)):
        if ii == 0:
            pass
        else:
            img_next = denoise(map.data, denoise_sigma, Triangle)

            prev_8bit = cv2.normalize(img_prev, None, 0, 255, cv2.NORM_MINMAX).astype('uint8')
            next_8bit = cv2.normalize(img_next, None, 0, 255, cv2.NORM_MINMAX).astype('uint8')
            array_of[:,:,:,ii] = cv2.calcOpticalFlowFarneback(prev_8bit, next_8bit, None, *of_params)

            img_prev = img_next
            
    return array_of

if __name__ == '__main__':
    eui_files = sorted(glob("/home/yjzhu/Solar/EIS_DKIST_SolO/src/EUI/HRI/euv174/20221024/coalign_step_boxcar/*.fits"))
    eui_map_seq_coalign = sunpy.map.Map(eui_files[:])
    eui_map_seq_coalign = MapSequenceCoalign(eui_map_seq_coalign)

    # of_east_1 = calc_boxcar_of(eui_map_seq_coalign, bottom_left=[500,600]*u.pix,top_right=[670,760]*u.pix)

    # with h5py.File("/home/yjzhu/Solar/EIS_DKIST_SolO/sav/optical_flow/of_east_1.h5", "w") as f:
    #     f.create_dataset("of_east_1", data=of_east_1)

    of_west_1 = calc_boxcar_of(eui_map_seq_coalign, bottom_left=[1700,500]*u.pix,top_right=[1970,800]*u.pix)

    with h5py.File("/home/yjzhu/Solar/EIS_DKIST_SolO/sav/optical_flow/of_west_1.h5", "w") as f:
        f.create_dataset("of_west_1", data=of_west_1)




