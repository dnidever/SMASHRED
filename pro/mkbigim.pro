pro mkbigim,file,clobber=clobber

if file_test(file) eq 0 then begin
  print,file,' NOT FOUND'
  return
endif

base = file_basename(file,'.fits')
outfile = base+'_bin.fits'

if file_test(outfile) eq 1 and not keyword_set(clobber) then begin
  print,outfile,' EXISTS and /clobber not set'
  return
endif

print,'Binning ',file

nxarr = lonarr(9)
nyarr = lonarr(9)
for i=0,8 do begin
  head = headfits(file,exten=i+1)
  nxarr[i] = sxpar(head,'naxis1')
  nyarr[i] = sxpar(head,'naxis2')
endfor

; Stitch everything together in ONE BIG IMAGE
; before binning otherwise there can be offsets
nx = long(total(nxarr[0:2]))
ny = nyarr[0]+nyarr[3]+nyarr[6]

xoff = [0,nxarr[0],nxarr[0]+nxarr[1]]
yoff = [0,nyarr[0],nyarr[0]+nyarr[3]]

;nx = min(nxarr)/10
;ny = min(nyarr)/10
print,'NX=',strtrim(nx,2),' NY=',strtrim(ny,2)

;nx = 1094L
;ny = 943L

;imarr = fltarr(9,nx,ny)
bigim = fltarr(nx,ny)
for i=0,8 do begin
  print,'Reading extension ',strtrim(i+1,2)
  fits_read,file,im,head,exten=i+1
  col = i mod 3
  row = i/3
  x0 = xoff[col]
  x1 = x0+nxarr[i]-1
  y0 = yoff[row]
  y1 = y0+nyarr[i]-1
  ;print,i,x0,x1,y0,y1
  bigim[x0:x1,y0:y1] = im
  ;imarr[i,*,*] = rebin(im[0:10*nx-1,0:10*ny-1],nx,ny)
endfor

;bigim = fltarr(3*nx,3*ny)   
;bigim[0:nx-1,0:ny-1] = reform(imarr[0,*,*])
;bigim[nx:2*nx-1,0:ny-1] = reform(imarr[1,*,*])
;bigim[2*nx:3*nx-1,0:ny-1] = reform(imarr[2,*,*])
;
;bigim[0:nx-1,ny:2*ny-1] = reform(imarr[3,*,*])
;bigim[nx:2*nx-1,ny:2*ny-1] = reform(imarr[4,*,*])
;bigim[2*nx:3*nx-1,ny:2*ny-1] = reform(imarr[5,*,*])
;
;bigim[0:nx-1,2*ny:3*ny-1] = reform(imarr[6,*,*])
;bigim[nx:2*nx-1,2*ny:3*ny-1] = reform(imarr[7,*,*])
;bigim[2*nx:3*nx-1,2*ny:3*ny-1] = reform(imarr[8,*,*])

nx2 = nx/10
ny2 = ny/10
bigim2 = rebin(bigim[0:nx2*10-1,0:ny2*10-1],nx2,ny2)

;displayc,bigim,/z

print,'Writing to ',outfile
MWRFITS,bigim2,outfile,/create
;MWRFITS,bigim,outfile,/create
;MWRFITS,imarr,outfile

end
