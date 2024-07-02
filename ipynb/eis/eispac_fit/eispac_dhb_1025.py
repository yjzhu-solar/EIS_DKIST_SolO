import eispac 
from eispac_fit import eispac_fit

if __name__ == "__main__":

    data_filepath = "../../../src/EIS/DHB_007_v2/20221025T0023/eis_20221025_014811.data.h5"

    fe_08_185_1c_template = "../../../src/EIS/EISPAC_templates/fe_08_185_213.1c.template.h5"
    fe_09_188_1c_template = "../../../src/EIS/EISPAC_templates/fe_09_188_497.1c.template.h5"
    fe_09_197_1c_template = "../../../src/EIS/EISPAC_templates/fe_09_197_862.1c.template.h5"
    fe_10_184_1c_template = "../../../src/EIS/EISPAC_templates/fe_10_184_536.1c.template.h5"
    fe_11_188_2c_template = "../../../src/EIS/EISPAC_templates/fe_11_188_216.2c.template.h5"
    fe_12_186_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_186_880.1c.template.h5"
    fe_12_192_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_192_394.1c.template.h5"
    fe_12_195_1c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.1c.template.h5"
    fe_12_195_2c_template = "../../../src/EIS/EISPAC_templates/fe_12_195_119.2c.template.h5"
    fe_13_202_1c_template = "../../../src/EIS/EISPAC_templates/fe_13_202_044.1c.template.h5"
    fe_13_203_2c_template = "../../../src/EIS/EISPAC_templates/fe_13_203_826.2c.template.h5"
    fe_14_264_1c_template = "../../../src/EIS/EISPAC_templates/fe_14_264_787.1c.template.h5"
    fe_15_284_1c_template = "../../../src/EIS/EISPAC_templates/fe_15_284_160.1c.template.h5"
    fe_15_284_2c_template = "../../../src/EIS/EISPAC_templates/fe_15_284_160.2c.template.h5"
    fe_16_262_1c_tempalte = "../../../src/EIS/EISPAC_templates/fe_16_262_984.1c.template.h5"

    o__04_279_1c_template = "../../../src/EIS/EISPAC_templates/o__04_279_933.1c.template.h5"
    o__06_184_1c_template = "../../../src/EIS/EISPAC_templates/o__06_184_117.1c.template.h5"
    mg_06_270_2c_template = "../../../src/EIS/EISPAC_templates/mg_06_270_394.2c.template.h5"
    si_07_275_1c_template = "../../../src/EIS/EISPAC_templates/si_07_275_368.1c.template.h5"
    mg_07_276_1c_template = "../../../src/EIS/EISPAC_templates/mg_07_276_153.1c.template.h5"
    mg_07_280_1c_template = "../../../src/EIS/EISPAC_templates/mg_07_280_737.1c.template.h5"

    si_10_258_1c_template = "../../../src/EIS/EISPAC_templates/si_10_258_375.1c.template.h5"
    si_10_261_1c_template = "../../../src/EIS/EISPAC_templates/si_10_261_058.1c.template.h5"

    templates = [fe_08_185_1c_template, fe_09_188_1c_template,
                fe_09_197_1c_template, fe_10_184_1c_template,
                fe_11_188_2c_template, fe_12_186_1c_template, 
                fe_12_192_1c_template, fe_12_195_1c_template, 
                fe_12_195_2c_template, fe_13_202_1c_template, 
                fe_13_203_2c_template, fe_14_264_1c_template, 
                fe_15_284_1c_template, fe_15_284_2c_template, 
                fe_16_262_1c_tempalte, o__04_279_1c_template, 
                o__06_184_1c_template, mg_06_270_2c_template, 
                si_07_275_1c_template, mg_07_276_1c_template, 
                mg_07_280_1c_template, si_10_258_1c_template, 
                si_10_261_1c_template]
    
    for template in templates:
        eispac_fit(data_filepath, template)


    

    
