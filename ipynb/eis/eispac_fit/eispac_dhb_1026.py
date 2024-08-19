import eispac 
from eispac_fit import eispac_fit
from glob import glob
import os

# data_filepaths = ["../../../src/EIS/DHB_007_v2/20221020T0343/",
#                   "../../../src/EIS/DHB_007_v2/20221020T2343/",
#                   "../../../src/EIS/DHB_007_v2/20221022T0630/",
#                   "../../../src/EIS/DHB_007_v2/20221022T2140/",
#                   "../../../src/EIS/DHB_007_v2/20221023T2038/",
#                   "../../../src/EIS/DHB_007_v2/20221025T0023/"]        

data_filepaths = ["../../../src/EIS/DHB_007_v2/20221026T0713/"] 

fe_12_195_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.1c.template.h5"
fe_12_195_2c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.2c.template.h5"
fe_13_202_1c_template = "../../../src/EIS/EISPAC_templates/fe_13_202_044.1c.template.h5"         

if __name__ == "__main__":

    for data_filepath in data_filepaths:
        data_files = sorted(glob(os.path.join(data_filepath,"*.data.h5"))) 

        for data_file in data_files:
            eispac_fit(data_file, fe_12_195_1c_template)
            # eispac_fit(data_filepath, fe_12_195_2c_template)
            eispac_fit(data_file, fe_13_202_1c_template)

    

    
