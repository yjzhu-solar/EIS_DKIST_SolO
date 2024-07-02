import eispac 
from eispac_fit import eispac_fit

if __name__ == "__main__":
    fe_08_186_1c_template = "../../../src/EIS/EISPAC_templates/fe_08_186_601.1c.template.h5"


    data_filepaths = ["../../../src/EIS/DHB_007_v2/20221025T0023/eis_20221025_014811.data.h5",
                      "../../../src/EIS/DHB_007_v2/20221020T2343/eis_20221021_010842.data.h5"]

    for data_filepath in data_filepaths:
        eispac_fit(data_filepath, fe_08_186_1c_template)

