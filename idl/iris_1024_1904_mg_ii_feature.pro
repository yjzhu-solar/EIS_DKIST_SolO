pro iris_1024_1904_mg_ii_feature

    iris_raster_file = "../src/IRIS/20221024/1904/iris_l2_20221024_190447_3643101203_raster_t000_r00000.fits"

    wave_corr = iris_prep_wavecorr_l2(iris_raster_file)
    wave_corr_nuv = wave_corr.corr_nuv
    
    iris_get_mg_features_lev2, iris_raster_file, 3, [-40, 40], lc, rp, bp, /onlyk

    save,filename="../src/IRIS/20221024/1904/fit_res/MgII_raster0.sav",lc,rp,bp,wave_corr_nuv
    
end