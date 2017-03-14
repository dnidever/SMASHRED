pro checkwcsexp,expnum

nexp = n_elements(expnum)
for e=0,nexp-1 do begin

  files = file_search(expnum[e]+'_??.fits',count=nfiles)
  bd = where(stregex(files,'a.fits$',/boolean) eq 1 or stregex(files,'s.fits$',/boolean) eq 1,nbd)
  if nbd gt 0 then remove,bd,files
  if nfiles eq 0 then begin
    print,'NO files for ',expnum[e]
    continue
  endif

  print,transpose(files)
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

    ra0=ten(sxpar(hdr,'TELRA'))*360./24
    dec0=ten(sxpar(hdr,'TELDEC'))

    if (i eq 0) then plot,[ra00,ra10,ra11,ra01,ra00]-ra0,[dec00,dec10,dec11,dec01,dec00]-dec0,xr=[-1.2,1.2]/cos(dec0*!pi/180),yr=[-1.2,1.2],xsty=1,ysty=1,tit=files[i] else $
      oplot,[ra00,ra10,ra11,ra01,ra00]-ra0,[dec00,dec10,dec11,dec01,dec00]-dec0
    xyouts,rac-ra0,decc-dec0,chip
  endfor
  if nexp gt 1 then foo=get_kbrd()
endfor

end

