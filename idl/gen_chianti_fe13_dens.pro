pro gen_chianti_fe13_dens

    n_dens = 26
    dens = 6 + findgen(n_dens)*0.2d
    n_heights = 26
    heights = 1 + findgen(n_heights)*0.02d

    FeXIII_1074_1079 = dblarr(n_dens, n_heights)
    FeXIII_202_203 = dblarr(n_dens, n_heights)

    for ii = 0, n_heights - 1 do begin

        FeXIII_data = emiss_calc("fe_13",6.25,dens=dens,radtemp=5770d,rphot=heights[ii])

        FeXIII_1074_id = where(FeXIII_data.lambda eq 10749.000)
        FeXIII_1079_id = where(FeXIII_data.lambda eq 10801.000)
        FeXIII_202_id = where(FeXIII_data.lambda eq 202.044)
        FeXIII_203_7_id = where(FeXIII_data.lambda eq 203.795)
        FeXIII_203_8_id = where(FeXIII_data.lambda eq 203.826)

        FeXIII_1074_1079[*,ii] = reform((FeXIII_data.em)[0,*,FeXIII_1074_id]/(FeXIII_data.em)[0,*,FeXIII_1079_id])
        FeXIII_202_203[*,ii] = reform((FeXIII_data.em)[0,*,FeXIII_202_id]/((FeXIII_data.em)[0,*,FeXIII_203_7_id] + (FeXIII_data.em)[0,*,FeXIII_203_8_id]))

    endfor


    save,filename="../sav/CHIANTI/FeXIII_dens_diag.sav",FeXIII_1074_1079,FeXIII_202_203,dens,heights


end 

