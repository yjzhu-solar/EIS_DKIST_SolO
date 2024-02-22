pro iris_1024_sgf_c_ii_1334, wvl_shift=wvl_shift, sav=sav

    iris_raster_file = "../src/IRIS/20221024/2322/iris_l2_20221024_232249_3600609177_raster_t000_r00000.fits"

    wd = iris_getwindata(iris_raster_file, 1334,/calib,wrange=[1333.5,1335.5])

    if n_elements(wvl_shift) eq 0 then begin
        wave_corr = iris_prep_wavecorr_l2(iris_raster_file)
        wave_corr_fuv = wave_corr.corr_fuv
        wvl_shift = 0d
    endif 

    eis_fit_template,wd,template
    iris_auto_fit, wd, fit, template=template
    print,"Old refwvl: ",fit.refwvl
    fit.refwvl = 1334.5323 + wvl_shift
    refwvl = fit.refwvl
    print,"New refwvl: ",fit.refwvl
    iris_fit_viewer, wd, fit

    if n_elements(sav) ne 0 then begin
        int = eis_get_fitdata(fit, /int,err=int_err)
        vel = eis_get_fitdata(fit, /vel,err=vel_err)
        wid = eis_get_fitdata(fit, /wid,err=wid_err)
        chi2 = fit.chi2
        save,filename="../src/IRIS/20221024/2322/fit_res/CII_1334_raster0.sav",int,int_err,vel,vel_err,wid,wid_err,chi2,refwvl,wave_corr_fuv
    endif
end 
