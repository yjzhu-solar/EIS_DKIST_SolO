import eispac 
from eispac_fit import eispac_fit

if __name__ == "__main__":

    data_filepath = "../../../src/EIS/20230330T1203/eis_20230330_120319.data.h5"

    fe_12_195_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.1c.template.h5"
    fe_10_184_1c_template = "../../../src/EIS/EISPAC_templates/fe_10_184_536.1c.template.h5"
    fe_08_185_1c_template = "../../../src/EIS/EISPAC_templates/fe_08_185_213.1c.template.h5"
    si_07_275_1c_template = "../../../src/EIS/EISPAC_templates/si_07_275_368.1c.template.h5"
    fe_09_188_1c_template = "../../../src/EIS/EISPAC_templates/fe_09_188_497.1c.template.h5"
    fe_09_197_1c_template = "../../../src/EIS/EISPAC_templates/fe_09_197_862.1c.template.h5"


    eispac_fit(data_filepath, fe_12_195_1c_template)
    eispac_fit(data_filepath, fe_10_184_1c_template, smooth_width=3)
    eispac_fit(data_filepath, fe_08_185_1c_template, smooth_width=3)
    eispac_fit(data_filepath, si_07_275_1c_template, smooth_width=3)
    eispac_fit(data_filepath, fe_09_188_1c_template, smooth_width=3)
    eispac_fit(data_filepath, fe_09_197_1c_template, smooth_width=3)


    

    
