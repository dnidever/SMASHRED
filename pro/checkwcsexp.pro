pro checkwcsexp,expnum

nexp = n_elements(expnum)
for e=0,nexp-1 do begin

  ;files = file_search(expnum[e]+'_??.fits',count=nfiles)
  files = file_search(expnum[e]+'_[0-9][0-9].fits',count=nfiles)
  if nfiles eq 0 then files = file_search(expnum[e]+'_[0-9].fits',count=nfiles)
  bd = where(stregex(files,'a.fits$',/boolean) eq 1 or stregex(files,'s.fits$',/boolean) eq 1,nbd)
  if nbd gt 0 then remove,bd,files
  if nfiles eq 0 then begin
    print,'NO files for ',expnum[e]
    continue
  endif

  ;print,transpose(files)
  for i=0,nfiles-1 do begin
    dum = strsplit(file_basename(files[i],'.fits'),'_',/extract)
    chip = reform(dum[1])
    hdr=headfits(files[i])
    xsz=sxpar(hdr,'NAXIS1')
    ysz=sxpar(hdr,'NAXIS2')
    xyad,hdr,0,0,ra00,dec00
    xyad,hdr,xsz,0,ra10,dec10
    xyad,hdr,xsz,ysz,ra11,dec11
    xyad,hdr,0,ysz,ra01,dec01
    xyad,hdr,xsz/2.,ysz/2.,rac,decc
    ; Get WCSFIT RMS
    indrms = where(stregex(hdr,'HISTORY WCSFIT: RMS',/boolean) eq 1,nindrms)
    if nindrms gt 0 then begin
      dum = strsplit(hdr[indrms[0]],' ',/extract)
      rms = float((strsplit(dum[2],'=',/extract))[1])
    endif else rms = 9999.99
    print,files[i],rms

    ra0=ten(sxpar(hdr,'TELRA',count=nra))*360./24
    if nra eq 0 then ra0=ten(sxpar(hdr,'RA',count=nra))
    dec0=ten(sxpar(hdr,'TELDEC',count=ndec))
    if ndec eq 0 then dec0=ten(sxpar(hdr,'DEC',count=nra))

    if (i eq 0) then plot,[ra00,ra10,ra11,ra01,ra00]-ra0,[dec00,dec10,dec11,dec01,dec00]-dec0,xr=[-1.2,1.2]/cos(dec0*!pi/180),yr=[-1.2,1.2],xsty=1,ysty=1,tit=files[i] else $
      oplot,[ra00,ra10,ra11,ra01,ra00]-ra0,[dec00,dec10,dec11,dec01,dec00]-dec0
    xyouts,rac-ra0,decc-dec0,chip
  endfor
  if nexp gt 1 then foo=get_kbrd()
endfor

end

