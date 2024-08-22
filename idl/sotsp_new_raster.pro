pro sotsp_new_raster

sotsp_file = "../src/SOTSP/20221026/lvl2/20221026_184140.fits"

read_sotsp,sotsp_file,index,data,scan_info=scan_info,/xycoord
slitmap=sotsp_rasterize(scan_info)

slitmap_dims=n_elements(slitmap)
new_dims=size(data,/dim)
new_dims[0]=slitmap_dims  ;  sets the correct size of the new data array
data_new=make_array(dim=new_dims,type=size(data,/type))
for ii=0,slitmap_dims[0]-1 do if slitmap[ii] ge 0 then data_new[ii,*,*]=data[slitmap[ii],*,*]


save,filename="../src/SOTSP/20221026/lvl2/sotsp_lvl2_missing_col_corrected.sav", data_new, slitmap

end