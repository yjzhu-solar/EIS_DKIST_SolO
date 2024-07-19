pro eis_1021_auto_fit

windata = eis_getwindata('../src/EIS/DHB_007_v2/20221020T2343/idl_l0_l1/eis_l1_20221021_010842.fits', 185.2, /refill)
windata = eis_bin_windata(windata,ybin=4)
eis_fit_template, windata, template
eis_wvl_select, windata, wvl_select
eis_auto_fit,windata, fit_data, template=template, wvl_select=wvl_select 
vel_feviii_185 = eis_get_fitdata(fit_data,/vel)

windata = eis_getwindata('../src/EIS/DHB_007_v2/20221020T2343/idl_l0_l1/eis_l1_20221021_010842.fits', 188.495, /refill)
windata = eis_bin_windata(windata, ybin=4)
eis_fit_template, windata, template
eis_wvl_select, windata, wvl_select
eis_auto_fit,windata, fit_data, template=template, wvl_select=wvl_select
vel_feix_188 = eis_get_fitdata(fit_data,/vel)

windata = eis_getwindata('../src/EIS/DHB_007_v2/20221020T2343/idl_l0_l1/eis_l1_20221021_010842.fits', 197.862, /refill)
windata = eis_bin_windata(windata, ybin=4)
eis_fit_template, windata, template
eis_wvl_select, windata, wvl_select
eis_auto_fit,windata, fit_data, template=template, wvl_select=wvl_select
vel_feix_197 = eis_get_fitdata(fit_data,/vel)

windata = eis_getwindata('../src/EIS/DHB_007_v2/20221020T2343/idl_l0_l1/eis_l1_20221021_010842.fits', 184.6, /refill)
windata = eis_bin_windata(windata, ybin=4)
eis_fit_template, windata, template
eis_wvl_select, windata, wvl_select
eis_auto_fit,windata, fit_data, template=template, wvl_select=wvl_select
vel_fex_184 = eis_get_fitdata(fit_data,/vel)

save, filename='../src/EIS/DHB_007_v2/20221020T2343/idl_l0_l1/eis_auto_fit_res_fe08_09_10.sav', vel_feviii_185,$
 vel_feix_188, vel_feix_197, vel_fex_184

end