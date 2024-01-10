from eispac_fit_plot_el_dhb import eispac_fit_plot
import os 
from glob import glob   



file_dirs = ["../../../src/EIS/EL_DHB_01/20221023T0044/"]

for ii, file_dir in enumerate(file_dirs):
    fitres_files = sorted(glob(os.path.join(file_dir,"fitres","*.fe_12_195_119.1c-0.fit.h5")))

    for jj, fitres_file in enumerate(fitres_files[:]):
        eispac_fit_plot(fitres_file,save_dir="../../../figs/EIS/quicklook_videos/EL_DHB_01/20221023/",index=jj,
                        figsize=(8,6))