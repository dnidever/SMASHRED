pro compare_trilegal_2cds_corrected

; Compare TRILEGAL simulation to our data in color-color diagrams
; with the corrections

field = 'Field140'

plotsdir = '/data/smash/cp/red/photred/trilegal/plots/'

tri = mrdfits('/data/smash/cp/red/photred/trilegal/'+field+'_trilegal_decam.fits.gz',1)
gdtri = where(tri.des_g le 22.0,ngdtri)
tri0 = tri
tri = tri[gdtri]
ttags = tag_names(tri)
ntri = n_elements(tri)
obj = mrdfits('/data/smash/cp/red/photred/catalogs/final/v3/'+field+'_combined_allobj_bright.fits',1)
otags = tag_names(obj)
nobj = n_elements(obj)

; Apply the corrections
otri = tri
; u=0.03094 -0.02137 (normal), u-g color
ucoef = [0.03094d0, -0.02137d0]
tri.decam_u -= poly(otri.decam_u-otri.des_g,ucoef)
; g=0.00045 -0.00311 (normal), g-r color
;gcoef = [0.00045d0, -0.00311d0]
gcoef = [0.20d0, -0.20d0]  ; this makes g-z vs. r-i work
tri.des_g -= poly(otri.des_g-otri.des_r,gcoef)
; r=0.00219 -0.06581 (normal), g-r color
rcoef = [0.00219d0, -0.06581d0]
tri.des_r -= poly(otri.des_g-otri.des_r,rcoef)
; i=-0.00426 -0.31326 (normal), i-z color
icoef = [-0.00426d0, -0.31326d0]
tri.des_i -= poly(otri.des_i-otri.des_z,icoef)
; z=0.00029 -0.18814 (normal), i-z color
zcoef = [0.00029d0, -0.18814d0]
tri.des_z -= poly(otri.des_i-otri.des_z,zcoef)

color1 = ['g-i','g-i','g-i','g-i','r-i']
color2 = ['u-g','g-r','r-i','i-z','g-z']
;color1 = ['r-i']
;color2 = ['g-z']
for i=0,n_elements(color1)-1 do begin

  col1band1 = first_el(strsplit(color1[i],'-',/extract))
  col1band2 = first_el(strsplit(color1[i],'-',/extract),/last)
  col2band1 = first_el(strsplit(color2[i],'-',/extract))
  col2band2 = first_el(strsplit(color2[i],'-',/extract),/last)
  indc1b1 = where(otags eq strupcase(col1band1),nindc1b1)
  indc1b2 = where(otags eq strupcase(col1band2),nindc1b2)
  indc2b1 = where(otags eq strupcase(col2band1),nindc2b1)
  indc2b2 = where(otags eq strupcase(col2band2),nindc2b2)

  tindc1b1 = where(ttags eq 'DES_'+strupcase(col1band1),ntindc1b1)
  if col1band1 eq 'u' then tindc1b1 = where(ttags eq 'DECAM_U',ntindc1b1)
  tindc1b2 = where(ttags eq 'DES_'+strupcase(col1band2),ntindc1b2)
  if col1band2 eq 'u' then tindc1b2 = where(ttags eq 'DECAM_U',ntindc1b2)
  tindc2b1 = where(ttags eq 'DES_'+strupcase(col2band1),ntindc2b1)
  if col2band1 eq 'u' then tindc2b1 = where(ttags eq 'DECAM_U',ntindc2b1)
  tindc2b2 = where(ttags eq 'DES_'+strupcase(col2band2),ntindc2b2)
  if col2band2 eq 'u' then tindc2b2 = where(ttags eq 'DECAM_U',ntindc2b2)

  col1 = obj.(indc1b1)-obj.(indc1b2)
  col2 = obj.(indc2b1)-obj.(indc2b2)
  si1 = sort(col1) & si2 = sort(col2)
  tcol1 = tri.(tindc1b1)-tri.(tindc1b2)
  tcol2 = tri.(tindc2b1)-tri.(tindc2b2)
  tsi1 = sort(tcol1) & tsi2 = sort(tcol2)
  xr = [ col1[si1[round(0.01*nobj)]] < tcol1[tsi1[round(0.01*ntri)]],$
         col1[si1[round(0.99*nobj)]] > tcol1[tsi1[round(0.99*ntri)]] ]
  xr = [xr[0]-0.1, xr[1]+0.1]
  ;xr = [-1.5 > min([col1,tcol1]), 6.0 < max([col1,tcol1])]
  xr = round(10*xr)/10.  ; round to nearest 10th
  dx = round( 100.*range(xr)/100. )/100.  ; round to nearest 100th
  yr = [ col2[si2[round(0.01*nobj)]] < tcol2[tsi2[round(0.01*ntri)]],$
         col2[si2[round(0.99*nobj)]] > tcol2[tsi2[round(0.99*ntri)]] ]
  yr = [yr[0]-0.1, yr[1]+0.1]
  ;yr = [-1.5 > min([col2,tcol2]), 6.0 < max([col2,tcol2])]
  yr = round(10*yr)/10.  ; round to nearest 10th
  dy = round( 100.*range(yr)/100. )/100.  ; round to nearest 100th
  charsize = 1.5

  hess,col1,col2,dum,oim,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
  oim2 = gsmooth(oim,1)
  hess,tcol1,tcol2,dum,tim,dx=dx,dy=dy,xr=xr,yr=yr,/noplot

  ; Make the plot
  !p.font = 0
  setdisp
  file = plotsdir+'trilegal_2cds_'+field+'_'+repstr(color1[i],'-','')+'_'+repstr(color2[i],'-','')+'_corrg2'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=10,ysize=11
  displayc,oim,xarr,yarr,xtit=color1[i],ytit=color2[i],tit='SMASH data',/log,$
           posim=[0.08,0.58,0.98,0.90],poscol=[0.08,0.97,0.98,0.99],charsize=charsize
  ; this contour is just so it will overplot properly below
  contour,alog10(oim2>1),xarr,yarr,nlevels=5,color=255,position=[0.08,0.08,0.98,0.40],/noerase
  displayc,tim,xarr,yarr,xtit=color1[i],ytit=color2[i],tit='TRILEGAL Simulation Corrected',/log,/noerase,$
           posim=[0.08,0.08,0.98,0.40],poscol=[0.08,0.47,0.98,0.49],charsize=charsize
  contour,alog10(oim2>1),xarr,yarr,nlevels=5,/over,color=255
  ps_close
  ps2jpg,file+'.eps',/eps

  ;stop
endfor


stop

end
