pro fixwcs,flist

readcol,flist,files,f='a'
dpos1=strpos(files,'_')
chip=files
for i=0,n_elements(files)-1 do chip[i]=strmid(files[i],dpos1[i]+1,2)

for i=0,n_elements(files)-1 do begin
  dpos=strpos(files[i],'/F')
  rdir=strmid(files[i],0,dpos+1)
  froot=strmid(files[i],dpos+1,3)
  img=mrdfits(files[i],0,hdr)
  extast,hdr,ain
  fref=file_search(rdir+froot+'*'+chip[i]+'.fits')
  match,fref,files,ela,elb
  fref[ela]='-1'
  use=where(fref ne '-1')
  if (use[0] ne -1) then begin
    href=headfits(fref[use[0]])
    extast,href,aout
    ra0=ten(sxpar(href,'TELRA'))*360./24
    dec0=ten(sxpar(href,'TELDEC'))
    aout.crval=aout.crval-[ra0,dec0]
  endif else restore,file='/d0/smash/data/cal/refast_'+chip[i]+'.sav'
  afix=aout
  ra0=ten(sxpar(hdr,'TELRA'))*360./24
  dec0=ten(sxpar(hdr,'TELDEC'))
  afix.crval=aout.crval+[ra0,dec0]
  afix.dateobs=ain.dateobs
  afix.mjdobs=ain.mjdobs
  hdr2=hdr
  putast,hdr2,afix
  writefits,files[i],img,hdr2
endfor



end
