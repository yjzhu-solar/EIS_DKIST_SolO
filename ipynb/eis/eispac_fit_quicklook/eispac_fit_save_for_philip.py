import eispac 
import h5py
import sunpy
import sunpy.map
from glob import glob
import os 


file_dirs = ["../../../src/EIS/HH_Flare/20221021T0315/",
             "../../../src/EIS/HH_Flare/20221021T0612/",
             "../../../src/EIS/HH_Flare/20221021T1814/",
             "../../../src/EIS/HH_Flare/20221022T0808/",
             "../../../src/EIS/HH_Flare/20221022T1118/",
             "../../../src/EIS/HH_Flare/20221022T1834/"]

path_to_save = "/home/yjzhu/Downloads/eis_int_vel_maps_for_philip/"

if not os.path.exists(path_to_save):
    os.makedirs(path_to_save)

for file_dir in file_dirs:
    fitres_files = sorted(glob(os.path.join(file_dir,"fitres","*.fe_12_195_119.1c-0.fit.h5")))

    for jj, fitres_file in enumerate(fitres_files[:]):

        fit_res = eispac.read_fit(fitres_file)
        int_map = fit_res.get_map(component=0, measurement="intensity")
        vel_map = fit_res.get_map(component=0, measurement="vel")

        int_meta_dict = {}

        for keys in ('cdelt1', 'cdelt2', 'crpix1', 'crpix2', 'crval1', 'crval2', 'date_obs'):
            int_meta_dict[keys] = int_map.meta[keys]

        filename_to_save = os.path.join(path_to_save,
                                        os.path.basename(fitres_file).replace(".fit.h5",".fit_int_vel_maps.h5"))

        with h5py.File(filename_to_save, "w") as hf:
            hf.create_dataset("intensity", data=int_map.data)
            hf.create_dataset("velocity", data=vel_map.data)
            hf.attrs.update(int_meta_dict)

        print(f"Saved {filename_to_save}")



