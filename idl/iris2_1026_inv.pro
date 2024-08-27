pro iris2_1026_inv, show=show

    raster_file = "../src/IRIS/20221026/0026/iris_l2_20221026_002630_3600609177_raster_t000_r00000.fits"
    iris2model = iris2(raster_file, delta_mu=0.2, weights_windows=[1.,1/2.,1.,1/3.], $
                    dir_save="../src/IRIS/20221026/0026/fit_res/")


    if n_elements(show) ne 0 then begin
        sel = show_iris2model(iris2model)
    endif
            

end