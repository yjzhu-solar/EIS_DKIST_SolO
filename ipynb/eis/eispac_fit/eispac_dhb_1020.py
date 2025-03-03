import eispac 
from eispac_fit import eispac_fit

if __name__ == "__main__":

    data_filepath = "../../../src/EIS/DHB_007_v2/20221020T0343/eis_20221020_034349.data.h5"

    # fe_12_195_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.1c.template.h5"
    # fe_12_195_2c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.2c.template.h5"
    # fe_13_202_1c_template = "../../../src/EIS/EISPAC_templates/fe_13_202_044.1c.template.h5"
    # fe_10_184_1c_template = "../../../src/EIS/EISPAC_templates/fe_10_184_536.1c.template.h5"
    # fe_12_186_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_186_880.1c.template.h5"
    # si_10_258_1c_template = "../../../src/EIS/EISPAC_templates/si_10_258_375.1c.template.h5"
    # si_10_261_1c_template = "../../../src/EIS/EISPAC_templates/si_10_261_058.1c.template.h5"
    fe_13_203_2c_template = "../../../src/EIS/EISPAC_templates/fe_13_203_826.2c.template.h5"

    # eispac_fit(data_filepath, fe_12_195_1c_template)
    # eispac_fit(data_filepath, fe_12_195_2c_template)
    # eispac_fit(data_filepath, fe_13_202_1c_template)
    # eispac_fit(data_filepath, fe_10_184_1c_template)
    # eispac_fit(data_filepath, fe_12_186_1c_template)
    # eispac_fit(data_filepath, si_10_258_1c_template)
    # eispac_fit(data_filepath, si_10_261_1c_template)

    eispac_fit(data_filepath, fe_13_203_2c_template)

    

    
