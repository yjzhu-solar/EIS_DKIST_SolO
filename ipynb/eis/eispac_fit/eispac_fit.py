import eispac
import os
from glob import glob

def eispac_fit(datafile, template_path, save_dir=None, ncpu="max", smooth_width=1,
               overwrite=False):
    fit_template = eispac.read_template(template_path)
    datacube = eispac.read_cube(datafile, fit_template.central_wave)
    if smooth_width > 1:
        datacube = datacube.smooth_cube(width=smooth_width)

    if datacube is None:
        raise Exception("The wavelength of the spectral line is not included in this dataset.")
    else:
        if save_dir is None:
            save_dir = os.path.dirname(datafile)

        template_basename = os.path.basename(template_path)
        data_basename = os.path.basename(datafile)
        fit_basename = data_basename.replace("data", ".".join(template_basename.split(".")[0:2]) + "-*.fit")
        
        fit_file_list = glob(os.path.join(save_dir, fit_basename))

        if len(fit_file_list) > 0 and not overwrite:
            print(f"Fit file {fit_file_list[0]} already exists. Skipping...")
        else:
            fit_res = eispac.fit_spectra(datacube, fit_template, ncpu=ncpu)
            eispac.save_fit(fit_res,save_dir=save_dir)


    


