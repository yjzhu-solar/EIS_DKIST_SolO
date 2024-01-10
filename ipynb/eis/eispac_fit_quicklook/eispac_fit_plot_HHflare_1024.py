from eispac_fit_plot import eispac_fit_plot
import os 
from glob import glob   



file_dirs = ["../../../src/EIS/HH_Flare/20221024T0800/","../../../src/EIS/HH_Flare/20221024T1808/"]

for ii, file_dir in enumerate(file_dirs):
    fitres_files = sorted(glob(os.path.join(file_dir,"fitres","*.fe_12_195_119.1c-0.fit.h5")))

    for jj, fitres_file in enumerate(fitres_files[:]):
        eispac_fit_plot(fitres_file,save_dir="../../../figs/EIS/quicklook_videos/HH_flare/20221024/",index=jj)