pro gen_eis_radcal

date_obs = '25-Oct-2022'

sw_wvl = findgen(2048)*0.022 + 165.1
lw_wvl = findgen(2048)*0.022 + 245.1

gdz_sw = eis_recalibrate_intensity_new(date_obs, sw_wvl, 1, /gdz)
gdz_lw = eis_recalibrate_intensity_new(date_obs, lw_wvl, 1, /gdz)

hpw_sw = eis_recalibrate_intensity_new(date_obs, sw_wvl, 1, /wul)
hpw_lw = eis_recalibrate_intensity_new(date_obs, lw_wvl, 1, /wul)

new_sw = dblarr(2048)
new_lw = dblarr(2048)

for i=0, 2047 do begin
    new_sw[i] = eis_recalibrate_intensity_new(date_obs, sw_wvl[i], 1)
    new_lw[i] = eis_recalibrate_intensity_new(date_obs, lw_wvl[i], 1)
endfor


save,filename='../sav/eis_radcal_20221025.sav', sw_wvl, lw_wvl, gdz_sw, gdz_lw, hpw_sw, hpw_lw, new_sw, new_lw

end