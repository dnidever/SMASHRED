pro checkwcs,field

files=file_search(field+'/'+field+'*[0-9].fits')
dpos1=strpos(files,field+'-')
expnum=files
for i=0,n_elements(files)-1 do expnum[i]=strmid(files[i],dpos1[i]+strlen(field)+1,8)
expu=expnum[uniq(expnum,sort(expnum))]
chip=expnum
for i=0,n_elements(files)-1 do begin
   dpos2=strpos(files[i],expnum[i]+'_')
   chip[i]=strmid(files[i],dpos2+strlen(expnum[i])+1,2)
endfor

for j=0,n_elements(expu)-1 do begin
  pk=where(expnum eq expu[j])
  print,transpose(files[pk])
  for i=0,n_elements(pk)-1 do begin
    hdr=headfits(files[pk[i]])
    xsz=sxpar(hdr,'NAXIS1')
    ysz=sxpar(hdr,'NAXIS2')
    xyad,hdr,0,0,ra00,dec00
    xyad,hdr,xsz,0,ra10,dec10
    xyad,hdr,xsz,ysz,ra11,dec11
    xyad,hdr,0,ysz,ra01,dec01
    xyad,hdr,xsz/2.,ysz/2.,rac,decc

    ra0=ten(sxpar(hdr,'TELRA'))*360./24
    dec0=ten(sxpar(hdr,'TELDEC'))

    if (i eq 0) then plot,[ra00,ra10,ra11,ra01,ra00]-ra0,[dec00,dec10,dec11,dec01,dec00]-dec0,xr=[-1.2,1.2]/cos(dec0*!pi/180),yr=[-1.2,1.2],xsty=1,ysty=1,tit=files[pk[i]] else oplot,[ra00,ra10,ra11,ra01,ra00]-ra0,[dec00,dec10,dec11,dec01,dec00]-dec0
    xyouts,rac-ra0,decc-dec0,chip[pk[i]]
  endfor
  foo=get_kbrd()
endfor

end

