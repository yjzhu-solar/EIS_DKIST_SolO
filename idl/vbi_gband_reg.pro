pro vbi_gband_reg, kernel_size = kernel_size

vbi_gband_dir = '../src/DKIST/vbi_1024/AEZDV/'
vbi_files = file_search(vbi_gband_dir + '*.fits')
n_files = n_elements(vbi_files)
vbi_gband_sav_dir = '../sav/DKIST_reg/AEZDV/'

if n_elements(kernel_size) eq 0 then begin
    kernel_size = 32
endif

for ii = 0, n_files - 2 do begin

    if ii eq 0 then begin
        read_sdo, vbi_files[ii], index_ref, im_ref, /USE_SHARED_LIB
        im_ref = im_ref[100:-100, 100:-100]
    endif else begin 
        im_ref = im_reg
    endelse

    
    read_sdo, vbi_files[ii+1], index_orig, im_orig, /USE_SHARED_LIB
    im_orig = im_orig[100:-100, 100:-100]

    kernel = bytarr(kernel_size, kernel_size)
    im_reg = reg(im_orig, im_ref, kernel)


    save, filename=vbi_gband_sav_dir + file_basename(vbi_files[ii], '.fits') + '.sav', im_ref

    if ii eq n_files - 2 then begin
        im_ref = im_reg
        save, filename=vbi_gband_sav_dir + file_basename(vbi_files[ii+1], '.fits') + '.sav', im_ref
    endif

    print, 'Done: ' + strtrim(ii,2) + ' of ' + strtrim(n_files - 2,2)
endfor

end