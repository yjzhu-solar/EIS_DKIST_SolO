pro sotsp_new_raster

sotsp_files = ["../src/SOTSP/20221026/lvl2/20221026_184140/20221026_184140.fits", $
                "../src/SOTSP/20221025/lvl2/20221025_184140/20221025_184140.fits", $
                "../src/SOTSP/20221024/lvl2/20221024_184138/20221024_184138.fits", $
                "../src/SOTSP/20221023/lvl2/20221023_184139/20221023_184139.fits"]

for ii = 0, n_elements(sotsp_files)-1 do begin
  sotsp_file = sotsp_files[ii]
  ; sotsp_file = "../src/SOTSP/20221023/lvl2/20221023_184139/20221023_184139.fits"
  ; sotsp_file = "../src/SOTSP/20221025/lvl2/20221025_184140/20221025_184140.fits"
  ; sotsp_file = "../src/SOTSP/20221026/lvl2/20221026_184140/20221026_184140.fits"
  print, sotsp_file

  read_sotsp,sotsp_file,index,data,scan_info=scan_info, phead=phead,/xycoord
  slitmap=sotsp_rasterize(scan_info)

  slitmap_dims=n_elements(slitmap)
  new_dims=size(data,/dim)
  new_dims[0]=slitmap_dims  ;  sets the correct size of the new data array
  data_new=make_array(dim=new_dims,type=size(data,/type))
  for jj=0,slitmap_dims[0]-1 do if slitmap[jj] ge 0 then data_new[jj,*,*]=data[slitmap[jj],*,*]

  save_filename = FILE_DIRNAME(sotsp_file) + '/' + $
    'sotsp_lvl2_missing_col_corrected.sav'

  save,filename=save_filename, data_new, slitmap, scan_info, index, data, phead

endfor
end