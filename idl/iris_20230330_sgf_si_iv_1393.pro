pro iris_20230330_sgf_si_iv_1393, wvl_shift=wvl_shift, sav=sav

    iris_raster_file = "../src/IRIS/20230330/iris_l2_20230330_090346_3400109477_raster_t000_r00003.fits"

    wd = iris_getwindata(iris_raster_file, 1393,/calib,wrange=[1393,1394.5])
    wd = eis_bin_windata(wd, xbin=2, ybin=2)

    if n_elements(wvl_shift) eq 0 then begin
        wave_corr = iris_prep_wavecorr_l2(iris_raster_file)
        wave_corr_fuv = wave_corr.corr_fuv
        wvl_shift = 0d
    endif 

    eis_fit_template,wd,template
    iris_auto_fit, wd, fit, template=template
    print,"Old refwvl: ",fit.refwvl
    fit.refwvl = 1393.755 + wvl_shift
    refwvl = fit.refwvl
    print,"New refwvl: ",fit.refwvl
    iris_fit_viewer, wd, fit

    if n_elements(sav) ne 0 then begin
        int = eis_get_fitdata(fit, /int,err=int_err)
        vel = eis_get_fitdata(fit, /vel,err=vel_err)
        wid = eis_get_fitdata(fit, /wid,err=wid_err)
        chi2 = fit.chi2
        save,filename="../src/IRIS/20230330/fit_res/SiIV_1393_raster0_rebin.sav",int,int_err,vel,vel_err,wid,wid_err,chi2,refwvl,wave_corr_fuv
    endif

end 
