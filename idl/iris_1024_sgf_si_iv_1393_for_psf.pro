pro iris_1024_sgf_si_iv_1393_for_psf, wvl_shift=wvl_shift, sav=sav

    iris_raster_files = file_search("../src/IRIS/20221024/2219/*raster*")
    print, "processing files:"
    print, iris_raster_files
    n_files = n_elements(iris_raster_files)

    for ii = 0, n_files - 1 do begin
        wd = iris_getwindata(iris_raster_files[ii], 1393,/calib,wrange=[1393,1394.5])

        if n_elements(wvl_shift) eq 0 then begin
            wave_corr = iris_prep_wavecorr_l2(iris_raster_files[ii])
            wave_corr_fuv = wave_corr.corr_fuv
            wvl_shift = 0d
        endif 

        if ii eq 0 then begin 
            eis_fit_template,wd,template
        endif 

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
            save,filename="../src/IRIS/20221024/2219/fit_res/SiIV_1393_raster_" + STRING(ii, FORMAT='(I1)') + ".sav",int,int_err,vel,vel_err,wid,wid_err,chi2,refwvl,wave_corr_fuv
        endif
    endfor 

end 
