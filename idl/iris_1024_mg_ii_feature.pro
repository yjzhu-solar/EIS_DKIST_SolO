pro iris_1024_mg_ii_feature

    iris_raster_file = "../src/IRIS/20221024/2322/iris_l2_20221024_232249_3600609177_raster_t000_r00000.fits"

    wave_corr = iris_prep_wavecorr_l2(iris_raster_file)
    wave_corr_nuv = wave_corr.corr_nuv
    
    iris_get_mg_features_lev2, iris_raster_file, 6, [-40, 40], lc, rp, bp

    save,filename="../src/IRIS/20221024/2322/fit_res/MgII_raster0.sav",lc,rp,bp,wave_corr_nuv
    
end