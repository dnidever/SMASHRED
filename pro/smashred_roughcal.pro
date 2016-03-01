pro smashred_roughcal,file

; Perform rough calibration using SLR and APASS
; input the "final" catalog of merged sources

dir = file_dirname(file)
base = file_basename(file,'.fits')
field = first_el(strsplit(base,'_',/extract))
if strmid(file,2,/reverse) eq '.gz' then base=file_basename(file,'.fits.gz')
print,'Calibrating ',file
cat = MRDFITS(file,1)

rar = minmax(cat.ra)
decr = minmax(cat.dec)

offset = {field:'',goff:99.99,gofferr:99.99,gnused:-1L,ioff:99.99,iofferr:99.99,inused:-1L,$
          uoff:99.99,uofferr:99.99,unused:-1L,roff:99.99,rofferr:99.99,rnused:-1L,$
          zoff:99.99,zofferr:99.99,znused:-1L}
offset.field = field

; Get APASS data for this field 
;-------------------------------
; get from disk or VizieR?
; /data/smash/apass, in 30 deg blocks of RA
apasslo = lindgen(12)*30
apasshi = lindgen(12)*30+30
lo = first_el(where(apasslo le rar[0]),/last)
hi = first_el(where(apasshi ge rar[1]))
n = hi-lo+1
ind = indgen(n)+lo
print,strtrim(n,2),' APASS files to load'
undefine,apass
for i=0,n-1 do begin
  apass1 = mrdfits('/data/smash/apass/apassdr7_ra'+strtrim(apasslo[ind[i]],2)+'-'+strtrim(apasshi[ind[i]],2)+'.fits.gz',1)
  push,apass,apass1
endfor
undefine,apass1

; Find the matches
dcr = 0.5
srcmatch,apass.ra,apass.dec,cat.ra,cat.dec,dcr,ind1,ind2,count=nmatch,/sph
print,strtrim(nmatch,2),' matches'
apass2 = apass[ind1]
cat2 = cat[ind2]

; Measure photometric offset in g-band
; V, B_V, B, GP, RP, IP, and "ERR" for each
; All negative errors means that the star was observed on one night
; only for that particular magnitude.  Therefore, the error is purely Poisson rather than nightly scatter.
sig = sqrt(apass2.gperr^2+cat2.gerr^2)
g_gd = where(cat2.g lt 50 and cat2.gerr lt 0.05 and abs(cat2.sharp) lt 1 and cat2.chi lt 3 and abs(apass2.gperr) lt 0.2,ng_gd)
if ng_gd gt 0 then begin
  robust_mean,apass2[g_gd].gp-cat2[g_gd].g,grobmn,grobsig,sig=sig[g_gd]
  print,'g-band mag offset (APASS-CAT)=',stringize(grobmn,ndec=4),'+/-',stringize(grobsig/sqrt(ng_gd),ndec=4)
  ; Applying the offset
  g_ind = where(cat.g lt 50,ng_ind)
  cat[g_ind].g += grobmn
  ; Save the values
  offset.goff = grobmn
  offset.gofferr = grobsig/sqrt(ng_gd)
  offset.gnused = ng_gd
  ; Plotting
  if keyword_set(pl) then begin
    xr = minmax(apass2[g_gd].gp)
    plot,apass2[g_gd].gp,apass2[g_gd].gp-(cat2[g_gd].g+grobmn),ps=3,xr=xr,yr=[-0.5,0.5],xs=1,ys=1,xtit='APASS gp',ytit='APASS gp - g'
    oplot,xr,[0,0],co=250
  endif
endif

; Do i-band as well
sig = sqrt(apass2.iperr^2+cat2.ierr^2)
i_gd = where(cat2.i lt 50 and cat2.ierr lt 0.05 and abs(cat2.sharp) lt 1 and cat2.chi lt 3 and abs(apass2.iperr) lt 0.2,ni_gd)
if ni_gd gt 0 then begin
  robust_mean,apass2[i_gd].ip-cat2[i_gd].i,irobmn,irobsig,sig=sig[i_gd]
  print,'i-band mag offset (APASS-CAT)=',stringize(irobmn,ndec=4),'+/-',stringize(irobsig/sqrt(ni_gd),ndec=4)
  ; Applying the offset
  i_ind = where(cat.i lt 50,ni_ind)
  cat[i_ind].i += irobmn
  ; Save the values
  offset.ioff = irobmn
  offset.iofferr = irobsig/sqrt(ni_gd)
  offset.inused = ni_gd
  ; Plotting
  if keyword_set(pl) then begin
    xr = minmax(apass2[i_gd].ip)
    plot,apass2[i_gd].ip,apass2[i_gd].ip-(cat2[i_gd].i+irobmn),ps=3,xr=xr,yr=[-0.5,0.5],xs=1,ys=1,xtit='APASS ip',ytit='APASS ip - i'
    oplot,xr,[0,0],co=250
  endif
endif

; Apply color correction (to other bands) using SLR
;---------------------------------------------------
gd = where(cat.g lt 22 and cat.gerr lt 0.05 and cat.i lt 50 and abs(cat.sharp) lt 1 and cat.chi lt 3,ngd)  ; get good sources
if ngd eq 0 then begin
  print,'No good g and i-band sources to calibrate'
  return
endif

; g-i is the fiducial color
; u-g
; g-r
; g-z

; polynomial coefficients derived from FieldB
gzcoef = [ -0.12239349d,  1.18763673d,  0.08909598d, -0.08486944d,  0.02408244d]  ; -0.35<g-i<2.9
grcoef = [0.01719993d,  0.63184690d,  0.01000558d,  0.13702503d, -0.09827074d,  0.01533803] ; -0.5<g-i<2.7
gucoef = [ 0.15109622d, -2.51319671d,  0.63049340d, -0.01046708]  ; 0.35<g-i<2.5

; how do we know which color to shift, horizontal or vertical?
; maybe need to calibrate i-band as well?

; g-z, -0.35<g-i<2.9, consistent with PARSEC DES isochrones
giz_gd = where(cat.g lt 22 and cat.gerr lt 0.05 and cat.i lt 50 and cat.z lt 50 and $
               ;cat.g-cat.i gt -0.3 and cat.g-cat.i lt 3.0 and $
               cat.g-cat.i gt 0.0 and cat.g-cat.i lt 2.0 and $  ; this region works better
               abs(cat.sharp) lt 1 and cat.chi lt 3,ngiz_gd)  ; get good sources
if ngiz_gd gt 10 then begin
  ;gzcoef2 = robust_poly_fit(cat[giz_gd].g-cat[giz_gd].i,cat[giz_gd].g-cat[giz_gd].z,4)
  sig = sqrt(cat[giz_gd].gerr^2+cat[giz_gd].ierr^2)
  gzmod = poly(cat[giz_gd].g-cat[giz_gd].i,gzcoef)
  ;gzoff = median(gzmod-(cat[giz_gd].g-cat[giz_gd].z))
  robust_mean,gzmod-(cat[giz_gd].g-cat[giz_gd].z),zrobmn,zrobsig,sig=sig
  print,'g-z color offset (FIDUCIAL-CAT)=',stringize(zrobmn,ndec=4),'+/-',stringize(zrobsig/sqrt(ngiz_gd),ndec=4)
  ; Applying the offset                                                                                              
  z_ind = where(cat.z lt 50,nz_ind)
  cat[z_ind].z -= zrobmn
  ; Save the values
  offset.zoff = -zrobmn  ; save the ADDITIVE offset
  offset.zofferr = zrobsig/sqrt(ngiz_gd)
  offset.znused = ngiz_gd
  ; Plotting
  if keyword_set(pl) then begin
    plot,cat[giz_gd].g-cat[giz_gd].i,cat[giz_gd].g-cat[giz_gd].z,ps=3,xr=[-0.5,2.5],yr=[-0.5,3],xs=1,ys=1,xtit='g-i',ytit='g-z'
    oplot,cat[giz_gd].g-cat[giz_gd].i,gzmod,ps=3,co=250
  endif
endif

; g-r, -0.5<g-i<2.7, offset from PARSEC DES isochrone by -0.77
gir_gd = where(cat.g lt 22 and cat.gerr lt 0.05 and cat.i lt 50 and cat.r lt 50 and $
               cat.g-cat.i ge 0.0 and cat.g-cat.i le 2.0 and $
               abs(cat.sharp) lt 1 and cat.chi lt 3,ngir_gd)  ; get good sources
if ngir_gd gt 10 then begin
  ;grcoef2 = robust_poly_fit(cat[gir_gd].g-cat[gir_gd].i,cat[gir_gd].g-cat[gir_gd].r,5)
  sig = sqrt(cat[gir_gd].gerr^2+cat[gir_gd].rerr^2)
  grmod = poly(cat[gir_gd].g-cat[gir_gd].i,grcoef)
  robust_mean,grmod-(cat[gir_gd].g-cat[gir_gd].r),rrobmn,rrobsig,sig=sig
  print,'g-r color offset (FIDUCIAL-CAT)=',stringize(rrobmn,ndec=4),'+/-',stringize(rrobsig/sqrt(ngir_gd),ndec=4)
  ; Applying the offset                                                                                              
  r_ind = where(cat.r lt 50,nr_ind)
  cat[r_ind].r -= rrobmn
  ; Save the values
  offset.roff = -rrobmn  ; save the ADDITIVE offset
  offset.rofferr = rrobsig/sqrt(ngir_gd)
  offset.rnused = ngir_gd
  ; Plotting
  if keyword_set(pl) then begin
    plot,cat[gir_gd].g-cat[gir_gd].i,cat[gir_gd].g-cat[gir_gd].r,ps=3,xr=[-0.5,2.5],yr=[-0.5,2],xs=1,ys=1,xtit='g-i',ytit='g-r'
    oplot,cat[gir_gd].g-cat[gir_gd].i,grmod,ps=3,co=250
  endif
endif

; g-u, 0.35<g-i<2.5, 
giu_gd = where(cat.g lt 22 and cat.gerr lt 0.05 and cat.i lt 50 and cat.u lt 50 and $
;               cat.g-cat.i ge 0.35 and cat.g-cat.i le 2.5 and $
               cat.g-cat.i ge 0.40 and cat.g-cat.i le 2.0 and $
              abs(cat.sharp) lt 1 and cat.chi lt 3,ngiu_gd)  ; get good sources
if ngiu_gd gt 10 then begin
  ;gucoef2 = robust_poly_fit(cat[giu_gd].g-cat[giu_gd].i,cat[giu_gd].g-cat[giu_gd].u,4)
  sig = sqrt(cat[giu_gd].gerr^2+cat[giu_gd].uerr^2)
  gumod = poly(cat[giu_gd].g-cat[giu_gd].i,gucoef)
  ;guoff = median(gumod-(cat[giu_gd].g-cat[giu_gd].u))
  robust_mean,gumod-(cat[giu_gd].g-cat[giu_gd].u),urobmn,urobsig,sig=sig
  print,'g-u color offset (FIDUCIAL-CAT)=',stringize(urobmn,ndec=4),'+/-',stringize(urobsig/sqrt(ngiu_gd),ndec=4)
  ; Applying the offset                                                                                              
  u_ind = where(cat.u lt 50,nu_ind)
  cat[u_ind].u -= urobmn
  ; Save the values
  offset.uoff = -urobmn  ; save the ADDITIVE offset
  offset.uofferr = urobsig/sqrt(ngiu_gd)
  offset.unused = ngiu_gd
  ; Plotting
  if keyword_set(pl) then begin
    plot,cat[giu_gd].g-cat[giu_gd].i,cat[giu_gd].g-cat[giu_gd].u,ps=3,xr=[0.0,2.5],yr=[-4,0],xs=1,ys=1,xtit='g-i',ytit='g-u'
    oplot,cat[giu_gd].g-cat[giu_gd].i,gumod,ps=3,co=250
  endif
endif

; There are some residual structure in all of these, but u-band is
; a VERY BAD FIT!
; Maybe refit things with the real data!


; Save the results
; catalog
outfile = dir+'/'+base+'_roughcal.fits'
print,'Writing roughly calibrated catalog to ',outfile
MWRFITS,cat,outfile,/create
; offset values
offsetfile = dir+'/'+base+'_roughcal_offsets.fits'
print,'Writing calibration offsets to ',offsetfile
MWRFITS,offset,offsetfile,/create

;stop

end
