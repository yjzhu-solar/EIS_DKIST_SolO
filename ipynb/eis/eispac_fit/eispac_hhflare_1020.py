import eispac 
from eispac_fit import eispac_fit
import os
from glob import glob

if __name__ == "__main__":

    file_dirs = ["../../../src/EIS/HH_Flare/20221020T0705/","../../../src/EIS/HH_Flare/20221020T1854/"]

    fe_12_195_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.1c.template.h5"

    for ii, file_dir in enumerate(file_dirs):
        data_files = sorted(glob(os.path.join(file_dir,"data_cube","*.data.h5")))

        for data_file in data_files:
            eispac_fit(data_file, fe_12_195_1c_template, save_dir=os.path.join(file_dir,"fitres"))

    


