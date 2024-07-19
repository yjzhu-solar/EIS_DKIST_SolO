pro xrt_calib

xrt_files = file_search("/home/yjzhu/Solar/EIS_DKIST_SolO/src/XRT/20221024/l1/L1*.fits")

n_files = n_elements(xrt_files)

for ii = 0, n_files - 1 do begin

read_xrt, xrt_files[ii], index, data
xrt_synleaksub, index, data, index_ll_corr, data_ll_corr
xrt_spotcor, index_ll_corr, data_ll_corr, index_refilled, data_refilled
xrt_deconvolve,index_refilled, data_refilled,index_deconv,data_deconv
write_xrt, index_deconv, data_deconv, outdir='/home/yjzhu/Solar/EIS_DKIST_SolO/src/XRT/20221024/l1_deconv/'

endfor

end 