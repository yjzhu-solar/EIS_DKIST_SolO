import eispac
import os

def eispac_fit(datafile, template_path, save_dir=None, ncpu="max", smooth_width=1):
    fit_template = eispac.read_template(template_path)
    datacube = eispac.read_cube(datafile, fit_template.central_wave)
    if smooth_width > 1:
        datacube = datacube.smooth_cube(width=smooth_width)

    if datacube is None:
        raise Exception("The wavelength of the spectral line is not included in this dataset.")
    else:
        fit_res = eispac.fit_spectra(datacube, fit_template, ncpu=ncpu)

        if save_dir is None:
            save_dir = os.path.dirname(datafile)

        eispac.save_fit(fit_res,save_dir=save_dir)


    


