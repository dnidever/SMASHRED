pro sim_astromerr

dir = '/data/smash/cp/red/photred/simim/'
cd,dir

inpbase = 'test'
loadals,inpbase+'.als',als
back = 1727.10
backsig = 24.1

; load the als input file
readline,inpbase+'.als.inp',inplines

; Redo the coordinates to put them on a grid
hd = headfits(inpbase+'.fits')
nx = sxpar(hd,'naxis1')
ny = sxpar(hd,'naxis2')
step = 100
nxg = nx/step
nyg = ny/step
x = lindgen(nxg)*step+step/2
y = lindgen(nyg)*step+step/2
xx = x#replicate(1,nyg)
yy = replicate(1,nxg)#y
als0 = als
als = als[0:nxg*nyg-1]
als.x = (xx)(*)
als.y = (yy)(*)
nals = n_elements(als)
; flat distribution in mag
coef = robust_poly_fitq(als0.mag,alog10(als0.err),2)
;als.mag = randomu(seed,nals)*8+10
als.mag = randomu(seed,nals)*10+10
als.err = 10^poly(als.mag,coef)

nmocks = 100 ;50
diff = fltarr(2,nals,nmocks)+!values.f_nan
sharp = fltarr(nals,nmocks)+!values.f_nan
undefine,allalsout,allalsin
for i=0,nmocks-1 do begin
  print,'Mock=',strtrim(i+1,2)
  base = 'test2'

  ; Copy files to temporary mock files
  file_copy,inpbase+['.psf','.opt','.als.opt'],base+['.psf','.opt','.als.opt'],/allow,/over

  ; Add random offset to X/Y positions
  als1 = als
  als1.x += randomu(seed,nals)*4.-2
  als1.y += randomu(seed,nals)*4.-2

  ; Make the image with 
  print,'Making synthetic image'
  file_delete,base+'.fits',/allow
  daophot_addstar,inpbase,als1,base+'.fits',/blank

  ; Add background level and noise
  print,'Adding background'
  fits_read,base+'.fits',im,head
  sz = size(im)
  im += back
  im += randomn(seed,sz[1],sz[2])*backsig
  ;modfits,base+'.fits',im,head
  mwrfits,im,base+'.fits',head,/create

  ; run find and aper
  print,'Running find/aper'
  file_delete,base+['.log','.coo','.ap','.psf'],/allow
  undefine,lines
  push,lines,'#!/bin/sh'
  push,lines,'daophot << END_DAOPHOT >> '+base+'.log'
  push,lines,'options'
  push,lines,base+'.opt'
  push,lines,' '
  push,lines,'attach '+base+'.fits'
  push,lines,'find'
  push,lines,'1,1'
  push,lines,base+'.coo'
  push,lines,'y'
  push,lines,'photometry'
  push,lines,' '
  push,lines,' '
  push,lines,base+'.coo'
  push,lines,base+'.ap'
  push,lines,'exit'
  push,lines,'END_DAOPHOT'
  writeline,base+'.dao.inp',lines
  file_chmod,base+'.dao.inp','755'o
  spawn,'./'+base+'.dao.inp'

  ; Make the allstar input file
  inplines1 = inplines  
  inplines1[11] = base+'.fits'
  inplines1[12] = base+'.psf'
  inplines1[13] = base+'.ap'
  inplines1[14] = base+'.als'
  inplines1[15] = base+'s.fits'

  ; run allstar
  print,'Running allstar'
  file_copy,inpbase+'.psf',base+'.psf',/allow,/over
  file_delete,base+'.als',/allow
  file_delete,base+'s.fits',/allow
  spawn,'allstar < '+base+'.als.inp >> '+base+'.log'

  ; Load the output als file
  loadals,base+'.als',alsout

  ; match to original catalogs
  srcmatch,als1.x,als1.y,alsout.x,alsout.y,10,ind1,ind2,count=nmatch
  print,strtrim(nmatch,2),' matches'

  ; save the differences
  diff[0,ind1,i] = als1[ind1].x-alsout[ind2].x
  diff[1,ind1,i] = als1[ind1].y-alsout[ind2].y
  sharp[ind1,i] = alsout[ind2].sharp

  ; save the input and output structures
  add_tag,als1,'mock',i+1,als1
  push,allalsin,als1
  add_tag,alsout,'mock',i+1,alsout
  push,allalsout,alsout

  ;stop

endfor

xsig = mad(reform(diff[0,*,*]),dim=2,/zero)
ysig = mad(reform(diff[1,*,*]),dim=2,/zero)
gd = where(finite(xsig) eq 1,ngd)
plot,als[gd].mag,xsig[gd],ps=1

; save,als,diff,xsig,ysig,file='sim_astromerr.dat'

;bindata,alog10(1.0/als[gd].err),alog10(xsig[gd]),xbin,ybin,binsize=0.2,min=0.5,max=3.0,gdind=gdind
bindata,alog10(1.0/als[gd].err),alog10(xsig[gd]),xbin,ybin,binsize=0.2,min=-0.4,max=3.0,gdind=gdind
plot,1.0/als[gd].err,xsig[gd],ps=1,yr=[0.0001,2],/ylog,/xlog
oplot,10^xbin,10^ybin,ps=-1,co=250,sym=3
coef = robust_poly_fit(xbin[gdind],ybin[gdind],1)
oplot,1.0/als.err,10^poly(alog10(1.0/als.err),coef),ps=3,co=200
print,coef
;     0.505401     -1.08208
fwhm = 3.774
gain = 3.91
;astrometric error should be ~FWHM/SNR
oplot,1.0/als.err,fwhm*als.err,ps=1,co=80

; the slope is a bit different, probably
; because of the background

; From Mendez et al. (2014), Eqn. 17
; source dominated (strong) case
; 0.644*FWHM/SNR
oplot,1.0/als.err,0.644*fwhm*als.err,ps=3,co=100
; Works well!!
; background dominated (weak) case
; 0.73685 * sqrt(sky)*FWHM^1.5 / Flux
flux = 10^((25-als.mag)*0.4)
oplot,1.0/als.err,0.73685 * sqrt(back)*FWHM^1.5 / flux,ps=3,co=150
; doesn't work well, low by ~2.7x, maybe I'm not probing the weak
; regime very well

; This works well at the low S/N end!!!
oplot,1.0/als.err,2*fwhm^1.5 * sqrt(back)/flux,ps=3,co=250

stop

end
