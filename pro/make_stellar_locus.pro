;+
;
; MAKE_STELLAR_LOCUS
;
; Define the stellar locus in the various bands.
;
; INPUTS:
;  str      A structure with SMASH photometry.
;  =gmax    The maximum g-value to use of the sources that define the
;             stellar locus.  The default is gmax=21.0.
;  /pl      Make some plots.
;
; OUTPUTS:
;  tstr     The structure of the stellar locus.
;
; USAGE:
;  IDL>make_stellar_locus,str,tstr
; 
; By D.Nidever  Nov 2016
;-

pro make_stellar_locus,str,tstr,gmax=gmax,pl=pl

; Make stellar locus templates

; Not enough inputs
if n_elements(str) eq 0 then begin
  print,'Syntax - make_stellar_locus,str,tstr,gmax=gmax,pl=pl'
endif

;if n_elements(str) eq 0 then str=mrdfits('~/decam/FieldB/FieldB_alf_calib_shapecuts_ebv.fits',1)
;if n_elements(str) eq 0 then str=mrdfits('~/smash/reduction/catalogs/final/v3/Field160_combined_allobj.fits.gz',1)  
;;g = where(str.imag lt 21 and str.gmag lt 50 and str.rmag lt 50)
;g = where(str.gmag lt 22,ng)    ; and str.gmag lt 50 and str.rmag lt 50)
if n_elements(gmax) eq 0 then gmax = 21  ; 22
;g = where(str.g lt gmax,ng)
g = where(str.g lt gmax and (str.r lt 50 or str.i lt 50) and abs(str.sharp) lt 1,ng)
print,'Max g = ',gmax

;  needs to be detected and gri, removes crap
str1 = str[g]
nstr = ng


; DEREDDEN the PHOTOMETRY
;  Extinction coefficients
;  u  4.239
;  g  3.303
;  r  2.285
;  i  1.698
;  z  1.263
;print,'DEREDDENING the photometry!!!!'
;str1.u -= str1.ebv*4.239
;str1.g -= str1.ebv*3.303
;str1.r -= str1.ebv*2.285
;str1.i -= str1.ebv*1.698
;str1.z -= str1.ebv*1.263

; Different ways to fit the stellar locus ridgeline
; Contour ridge line
; Poly
; Medfilt
; Bindata
; Medfilt in sorted list


dx = 0.01
bin = 0.10
gi0 = -1.0  ;0.30
gi1 = 3.50
;nbin = (gi1-gi0)/bin
nbin = (gi1-gi0)/dx + 1
tstr = replicate({gibin:0.0,ui:9999.0,uisig:9999.0,uiscatter:9999.0,gi:0000.0,gisig:9999.0,giscatter:9999.0,$
                  ri:9999.0,risig:9999.0,riscatter:9999.0,zi:9999.0,zisig:9999.0,ziscatter:9999.0,$
                  cols:fltarr(5)+9999,sigs:fltarr(5)+9999,scatters:fltarr(5)+9999,mederr:fltarr(5)+9999,$
                  n:0L},nbin)
gitemp = lindgen(nbin)*dx+gi0+0.5*dx
tstr.gibin = gitemp
tags = tag_names(tstr)

; Fit Bsplines to the three 2CDs

; u, g, r, i, z
mag = fltarr(ng,5)
mag[*,0] = str1.u
mag[*,1] = str1.g
mag[*,2] = str1.r
mag[*,3] = str1.i
mag[*,4] = str1.z
err = fltarr(ng,5)
err[*,0] = str1.uerr
err[*,1] = str1.gerr
err[*,2] = str1.rerr
err[*,3] = str1.ierr
err[*,4] = str1.zerr

bands = ['u','g','r','i','z']
mind = [0,2,4]  ; u, r, z
for i=0,2 do begin
  nord = 3 ;4
  bkspace = 0.3 ; 50
  gd = where( mag[*,1] lt 50 and mag[*,3] lt 50 and mag[*,mind[i]] lt 50 and $
              mag[*,1]-mag[*,3] gt gi0 and mag[*,1]-mag[*,3] lt gi1,ngd)
  if ngd eq 0 then goto,BOMB
  xx = mag[gd,1]-mag[gd,3]
  yy = mag[gd,mind[i]]-mag[gd,3]
  colerr = sqrt( err[gd,mind[i]]^2 + err[gd,3]^2 )
  invvar = 1.0/colerr^2
  dum = bspline_iterfit(xx,yy,invvar=invvar,nord=nord,bkspace=bkspace,yfit=model)
  diff = yy-model
  sig = mad(diff,/zero)
  if keyword_set(pl) then begin
    plot,xx,diff,ps=3,xr=[-1,4],yr=[-4,4]*sig
    oplot,[-10,10],[0,0],linestyle=2
  endif

  ; Measure sigs and scatters
  si = sort(xx)
  nbinsigs = 100 ;200
  nsigs = floor(ngd/nbinsigs)
  sigs = fltarr(nsigs)
  scatters = sigs
  mngi = sigs
  mederr = sigs
  for j=0L,nsigs-1 do begin
    lo = nbinsigs*j  &  hi=(lo+nbinsigs-1) < (ngd-1)
    sigs[j] = mad(diff[si[lo:hi]],/zero)
    scatters[j] = sqrt( sqrt(median(diff[si[lo:hi]]^2))^2 - sqrt(median(colerr[si[lo:hi]]^2))^2 )
    mngi[j] = mean(xx[si[lo:hi]])
    mederr[j] = median(colerr[si[lo:hi]])
  endfor

  ; Smooth
  smsigs = gsmooth(sigs,3)
  smscatters = gsmooth(scatters,3)
  smmederr = gsmooth(mederr,3)
  
  ; Make the templates
  interp,xx[si],model[si],gitemp,coltemp
  interp,mngi,smsigs,gitemp,sigstemp
  interp,mngi,smscatters,gitemp,scatterstemp
  interp,mngi,smmederr,gitemp,mederrtemp

  tind1 = where(tags eq strupcase(bands[mind[i]])+'I',ntind1)
  tstr.(tind1[0]) = coltemp
  tind2 = where(tags eq strupcase(bands[mind[i]])+'ISIG',ntind2)
  tstr.(tind2[0]) = sigstemp
  tind3 = where(tags eq strupcase(bands[mind[i]])+'ISCATTER',ntind3)
  tstr.(tind3[0]) = scatterstemp
  ;tstr.(mind[i]+1) = coltemp
  ;tstr.(mind[i]+2) = sigstemp
  ;tstr.(mind[i]+3) = scatterstemp

  tstr.cols[mind[i]] = coltemp
  tstr.sigs[mind[i]] = sigstemp
  tstr.scatters[mind[i]] = scatterstemp
  tstr.mederr[mind[i]] = mederrtemp
  tstr.n = ngd

  ;stop
  BOMB:
endfor

tstr.cols[1] = gitemp ; fill in g-i color
tstr.cols[3] = 0
tstr.gi = gitemp

; Interpolate scatter and sigs for g
tstr.sigs[1] = 0.5*(tstr.sigs[0]+tstr.sigs[2])
tstr.scatters[1] = 0.5*(tstr.scatters[0]+tstr.scatters[2])
tstr.gisig = tstr.sigs[1]
tstr.giscatter = tstr.scatters[1]
; Interpolate scatter and sigs for i
;tstr.cols[3] = 0.5*(tstr.cols[2]+tstr.cols[4])
tstr.sigs[3] = 0.5*(tstr.sigs[2]+tstr.sigs[4])
tstr.scatters[3] = 0.5*(tstr.scatters[2]+tstr.scatters[4])

; Make sure all scatters are non-negative
tstr.scatters = tstr.scatters > 0.0

;stop

; Plots
if keyword_set(pl) then begin
  plotc,(indgen(5)+1)#replicate(1,nbin),tstr.cols,replicate(1,5)#tstr.gibin,ps=3,xr=[0,6],yr=[-1,6]
  ;plot,[0],[0],/nodata,xr=[0,6],yr=[-5,5]
  ind = indgen(5)+1
  bottom = 50
  ncolors = 200
  color = SCALE( tstr.gibin, [min(tstr.gibin),max(tstr.gibin)], [bottom, bottom+ncolors-1])
  for i=0,nbin-1 do oplot,ind,tstr[i].cols,co=color[i]

  plotc,replicate(1,5)#indgen(nbin),tstr.sigs,ps=3,yr=[-0.1,0.5]
  oplot,indgen(nbin),tstr.uisig,co=80
  ;oplot,indgen(nbin),tstr.mederr[0],co=80
  oplot,indgen(nbin),tstr.scatters[0],co=80
  oplot,indgen(nbin),tstr.gisig,co=150
  ;oplot,indgen(nbin),tstr.mederr[1],co=150
  oplot,indgen(nbin),tstr.scatters[1],co=150
  oplot,indgen(nbin),tstr.risig,co=200
  ;oplot,indgen(nbin),tstr.mederr[2],co=200
  oplot,indgen(nbin),tstr.scatters[2],co=200
  oplot,indgen(nbin),tstr.zisig,co=250
  ;oplot,indgen(nbin),tstr.mederr[4],co=250
  oplot,indgen(nbin),tstr.scatters[4],co=250
endif

;stop

;setdisp
;!p.font = 0
;
;ps_open,'smash_stellar_templates_bspline_seds',/color,thick=4,/encap
;plotc,(indgen(5)+1)#replicate(1,nbin),tstr.cols,replicate(1,5)#tstr.gibin,ps=3,xr=[0.5,5.5],yr=[-1,6],xs=1,ys=1,$
;       xtickn=['u','g','r','i','z'],xtickv=[1,2,3,4,5],$
;      xtit='Band',ytit='Color (X-i)',tit='Empirical Stellar Templates (color-coded by g-i)'
;;plot,[0],[0],/nodata,xr=[0,6],yr=[-5,5]
;ind = indgen(5)+1
;bottom = 50
;ncolors = 200
;color = SCALE( tstr.gibin, [min(tstr.gibin),max(tstr.gibin)], [bottom, bottom+ncolors-1])
;for i=0,nbin-1,5 do oplot,ind,tstr[i].cols,co=color[i]
;ps_close
;ps2jpg,'smash_stellar_templates_bspline_seds.eps',/eps
;
;
;ps_open,'smash_stellar_templates_bspline_internalscatter',/color,thick=4,/encap
;plotc,[0],[0],/nodata,ps=3,xr=[-1,3.5],yr=[0.0,0.17],xs=1,ys=1,xtit='g-i',$
;      ytit='Internal Scatter',tit='Internal Scatter'   ;yr=[-0.1,1]
;;oplot,[0,1e4],[0,0],linestyle=2
;;oplot,tstr.gi,tstr.uisig,co=80
;oplot,tstr.gi,tstr.uiscatter,co=80
;;oplot,indgen(nbin),tstr.mederr[0],co=80
;;oplot,tstr.gi,tstr.gisig,co=150
;oplot,tstr.gi,tstr.giscatter,co=150
;;oplot,indgen(nbin),tstr.mederr[1],co=150
;;oplot,tstr.gi,tstr.risig,co=200
;oplot,tstr.gi,tstr.riscatter,co=200
;;oplot,indgen(nbin),tstr.mederr[2],co=200
;;oplot,tstr.gi,tstr.zisig,co=250
;oplot,tstr.gi,tstr.ziscatter,co=250
;;oplot,indgen(nbin),tstr.mederr[4],co=250
;al_legend,['u-i','g-i','r-i','z-i'],textcolor=[80,150,200,250],/top,/right,charsize=1.3
;ps_close
;ps2jpg,'smash_stellar_templates_bspline_internalscatter.eps',/eps
;

; There DEFINITELY a magnitude dependence in the U colors in the blue
; because you go from MW halo to LMC stellar pops.  Not sure if this
; is reddening or age/metallicity.  I think it's reddening and
; not the LMC.

; I don't think it's reddening.  It just shifts everything.

if keyword_set(pl) then begin

  ; Make 3-panel plot showing the stellar distribution in the four
  ; colors vs. g-i and overplotting the stellar template values
  ; Also see what the templates look like in the other color-color
  ; diagrams

  pos1 = [0,0,1,0.33]
  pos2 = [0,0.33,1,0.66]
  pos3 = [0,0.66,1,0.99]
  
  ;save = 1
  if keyword_set(save) then begin
    setdisp
    !p.font = 0
    file='smash_stellar_locus'
    ps_open,file,thick=4,/color,/encap
    device,/inches,xsize=8.0,ysize=10.0
    charsize=1.3
    x0 = 0.10
    x1 = 0.98
    yoff = 0.06
    dy = 0.27
    posim1 = [x0,yoff,x1,dy]
    poscol1 = [x0,0.33-0.02,x1,0.33-0.01]
    posim2 = [x0,0.33+yoff,x1,0.33+dy]
    poscol2 = [x0,0.66-0.02,x1,0.66-0.01]
    posim3 = [x0,0.66+yoff,x1,0.66+dy]
    poscol3 = [x0,0.99-0.02,x1,0.99-0.01]
 endif
   
  ; u-i vs. g-i
  undefine,dum,im
  hess,str1.g-str1.i,str1.u-str1.i,dum,im,dx=0.05,dy=0.1,xr=[-1,3.5],yr=[-1,6],xarr=xarr,yarr=yarr,/noplot
  displayc,im,xarr,yarr,posim=posim1,poscol=poscol1,xtit='g-i',ytit='u-i',/log,charsize=charsize,maskv=0,maskc=255
  oplot,tstr.gi,tstr.ui,co=250

  ; r-i vs. g-i
  undefine,dum,im
  hess,str1.g-str1.i,str1.r-str1.i,dum,im,dx=0.05,dy=0.1,xr=[-1,3.5],yr=[-1,3],xarr=xarr,yarr=yarr,/noplot
  displayc,im,xarr,yarr,posim=posim2,poscol=poscol2,xtit='g-i',ytit='r-i',/log,charsize=charsize,/noerase,maskv=0,maskc=255
  oplot,tstr.gi,tstr.ri,co=250
  ; MISSING THE GIANTS IN THE RED!!!

  ; z-i vs. g-i
  undefine,dum,im
  hess,str1.g-str1.i,str1.z-str1.i,dum,im,dx=0.05,dy=0.07,xr=[-1,3.5],yr=[-2,1],xarr=xarr,yarr=yarr,/noplot
  displayc,im,xarr,yarr,posim=posim3,poscol=poscol3,xtit='g-i',ytit='z-i',/log,charsize=charsize,/noerase,$
           xticklen=0.04,maskv=0,maskc=255
  oplot,tstr.gi,tstr.zi,co=250

  if keyword_set(save) then begin
    ps_close
    ps2jpg,file+'.eps',/eps
  endif
  
endif


;; Other 2CDS
;
;; g-z vs. g-i
;hess,str1.g-str1.i,str1.g-str1.z,dx=0.01,dy=0.01,/log,yr=[-2,6]
;oplot,tstr.gi,tstr.gi-tstr.zi,ps=8,co=250
;
;; u-i vs. u-g
;hess,str1.u-str1.g,str1.u-str1.i,dx=0.01,dy=0.01,/log,xr=[-1,3],yr=[-2,6]
;oplot,tstr.ui-tstr.gi,tstr.ui,ps=8,co=250
;; works well
;
;; g-r vs. u-g
;hess,str1.u-str1.g,str1.g-str1.r,dx=0.02,dy=0.03,/log,xr=[0,3],yr=[-1,2]
;oplot,tstr.ui-tstr.gi,tstr.gi-tstr.ri,ps=8,co=250
;; it doesn't capture the spread in u-g at high g-r, maybe it's okay
;
;; r-i vs. u-g
;hess,str1.u-str1.g,str1.r-str1.i,dx=0.02,dy=0.03,/log,xr=[0,3],yr=[-1,2]
;oplot,tstr.ui-tstr.gi,tstr.ri,ps=8,co=250
;; works well
;
;; g-i vs. g-r
;hess,str1.g-str1.r,str1.g-str1.i,dx=0.02,dy=0.03,/log,xr=[-1,3],yr=[-1,4]
;oplot,tstr.gi-tstr.ri,tstr.gi,ps=8,co=250
; MISSES giants in the red!!!

; Overlap, except for the red giants in some colors, the stellar
; templates work and capture the multi-band stellar locus!!

;g = where(str.i lt 23)
;plotc,str[g].g-str[g].i,str[g].u-str[g].i,str[g].g,ps=3,xr=[-1,3.5],yr=[0,6]  
;plotc,str[g].g-str[g].ebv*3.303-(str[g].i-str[g].ebv*1.698),str[g].u-str[g].ebv*4.239-(str[g].i-str[;g].ebv*1.698),str[g].g,ps=3,xr=[-1,3.5],yr=[0,6]

if keyword_set(stp) then stop

end
