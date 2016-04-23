pro get_transphot_night,arr,mfitstr,verbose=verbose,silent=silent,pl=pl

; Fit color and zpterm for each night+chip and airmass term
; per night. itertate until convergence

undefine,mfitstr

narr = n_elements(arr)
imjd = arr[0].mjd

; Get unique chip tags
ui = uniq(arr.chip,sort(arr.chip))
uchips = arr[ui].chip
nchips = n_elements(uchips)

; Initialize the fitting structure
mfitstr = replicate({mjd:-1L,chip:-1L,nightnum:0,nstars:-1L,zpterm:999.0,zptermerr:999.0,$
                     colterm:999.0,coltermerr:999.0,rms:999.0,sig:999.0,chisq:999.0,$
                     nbrt:0L,brtrms:999.0,brtsig:999.0,brtchisq:999.0,status:-1,amterm:0.0,amtermerr:999.0},nchips)

; Iterate until converge
endflag = 0
count = 0L
lastmfitstr = mfitstr
WHILE (endflag eq 0) do begin

  ; Fit color term and zpterm for each chip separately
  ;---------------------------------------------------

  resid = fltarr(narr)+999.0
  for i=0,nchips-1 do begin
    ichip = uchips[i]
    MATCH,arr.chip,uchips[i],ind1,ind2,count=nmatch,/sort

    mfitstr[i].mjd = imjd
    mfitstr[i].chip = ichip
    mfitstr[i].nstars = nmatch

    if nmatch gt 2 then begin

      ; Fit the color term
      ; CMAG - 6
      ; COL - 8
      x = arr[ind1].(8)
      mag = arr[ind1].mag
      y = arr[ind1].mag-arr[ind1].(6)
      ywam = y   ; with extinction still in 
      y -= arr[ind1].airmass*mfitstr[i].amterm ; subtract AIRMASS term
      err = arr[ind1].err

      ;coef0 = poly_fit(x,y,1,measure_errors=err,sigma=sigma0,yerror=yerror0,status=status0,yfit=yfit0)
      coef0 = robust_poly_fitq(x,y,1)
      yfit0 = poly(x,coef0)
      yerror0 = sqrt(mean((y-yfit0)^2))
      ; reject outliers and refit
      diff0 = y-yfit0
      gd = where(abs(diff0) lt 3*yerror0 and err lt 0.05,ngd)
      coef = robust_poly_fitq(x[gd],y[gd],1)  ; use robust fit for final values
      yfit = poly(x,coef)
      ydiff = y-yfit
      ; use poly_fit to get uncertainties
      coef1 = poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1)

      resid[ind1] = ywam-yfit  ; save the residuals WITH extinction portion 

      chisq = mean((ydiff/err)^2)

      mfitstr[i].zpterm = coef[0]
      mfitstr[i].zptermerr = coeferr[0]
      mfitstr[i].colterm = coef[1]
      mfitstr[i].coltermerr = coeferr[1]
      mfitstr[i].rms = yerror
      mfitstr[i].sig = mad(ydiff)
      mfitstr[i].chisq = chisq
      mfitstr[i].status = status

      ; RMS of bright stars
      gbright = where(err lt 0.05 and mag lt min(mag)+3,nbright)
      mfitstr[i].nbrt = nbright
      if nbright gt 1 then begin
        brtrms = sqrt(mean(ydiff[gbright]^2))
        brtchisq = mean((ydiff[gbright]/err[gbright])^2)
        mfitstr[i].brtrms = brtrms
        mfitstr[i].brtsig = mad(ydiff[gbright])
        mfitstr[i].brtchisq = brtchisq
      endif

      if keyword_set(verbose) then print,strtrim(i+1,2),' ',ichip,' ',nmatch,' ',coef[0],' ',coef[1],' ',yerror,' ',brtrms
    endif
  endfor  ; chip loop

  ; Fit airmass term for all data
  ;------------------------------
  x = arr.airmass
  y = resid
  mag = arr.mag
  err = arr.err
  ;mfitstr[i].minam = min(x)
  ;mfitstr[i].maxam = max(x)

  ;coef0 = poly_fit(x,y,1,measure_errors=err,sigma=sigma0,yerror=yerror0,status=status0,yfit=yfit0)
  coef0 = robust_poly_fitq(x,y,1)
  yfit0 = poly(x,coef0)
  yerror0 = sqrt(mean((y-yfit0)^2))
  ; reject outliers and refit
  ydiff0 = y-yfit0
  gd = where(abs(ydiff0) lt 3*yerror0 and err lt 0.05,ngd)
  if yerror0 lt 0.0001 then gd = where(abs(ydiff0) lt 3*0.01 and err lt 0.05,ngd)
  coef = robust_poly_fitq(x[gd],y[gd],1)
  yfit = poly(x,coef)
  ydiff = y-yfit
  ; use poly_fit to get uncertainties
  coef1 = poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1)

  mfitstr.amterm = coef[1]
  mfitstr.amtermerr = coeferr[1]

  ;resid2[ind1] = diff  ; save the residuals

  nrms = sqrt(mean((ydiff/err)^2))

; CHECK THAT I'M USING THE CORRECT SIGN FOR THESE TERMS!!!!!!

  ; Have we converged
  ; difference in zpterm, colterm, airmass terms
  maxdiff_thresh = 0.0001
  maxiter = 10
  diff = [ abs(mfitstr.zpterm-lastmfitstr.zpterm), $
           abs(mfitstr.colterm-lastmfitstr.colterm), $
           abs(mfitstr[0].amterm-lastmfitstr[0].amterm) ]
  if max(diff) lt maxdiff_thresh or count ge maxiter then endflag=1
  if not keyword_set(silent) then print,count,median(mfitstr.zpterm),median(mfitstr.colterm),mfitstr[0].amterm,$
          max(diff[0:nchips-1]),max(diff[nchips:2*nchips-1]),max(diff[2*nchips]),format='(I5,6F11.5)'

  ; save last values
  lastmfitstr = mfitstr

  count++

ENDWHILE

if keyword_set(pl) then begin
  plotc,x,ydiff,err,ps=3,yr=[-2,2],xtit='Airmass',ytit='Residuals',tit='Color-coded by error'
  oplot,[0,10],[0,0],linestyle=2
endif

;stop

end

;------

pro get_transphot,file

; Check the STDRED data to see if we need separate color terms
; for each chip

if n_elements(file) eq 0 then begin
  print,'Syntax - get_transphot,file'
  return
endif

arr = MRDFITS(file,1)
narr = n_elements(arr)

; Get unique nights
ui = uniq(arr.mjd,sort(arr.mjd))
umjds = arr[ui].mjd
nmjds = n_elements(umjds)

; Loop through nights and fit the color and airmass terms
mstr = replicate({mjd:0L,nightnum:0,date:'',nchips:0L,nstars:0L,rms:99.0,brtrms:99.0,zpterm:99.0,zptermsig:99.99,$
                  colterm:99.0,coltermsig:99.99,amterm:99.0,amtermsig:99.99,badsoln:-1,photometric:-1},nmjds)
observatory,'ctio',obs
undefine,fitstr
print,'  NUM    MJD     DATE     NSTARS    ZPTERM    COLTERM     AMTERM       RMS      BRTRMS    MEDSIG'
for i=0,nmjds-1 do begin

  MATCH,arr.mjd,umjds[i],ind1,ind2,count=nmatch,/sort
  arr1 = arr[ind1]
  GET_TRANSPHOT_NIGHT,arr1,mfitstr,/silent
  mfitstr.nightnum = i+1

  ; Measure total RMS and NRMS per night
  totresid = mfitstr.nstars*mfitstr.rms^2
  rms = sqrt(mean(total(totresid)/nmatch))  ; rms for all stars of this night
  totbrtresid = mfitstr.nbrt*mfitstr.brtrms^2
  brtrms = sqrt(mean(total(totbrtresid)/total(mfitstr.nbrt)))
  medsig = median(mfitstr.sig)

  ; Save information to MJD structure
  mstr[i].mjd = umjds[i]
  ui = uniq(mfitstr.chip,sort(mfitstr.chip))
  mstr[i].nchips = n_elements(ui)
  mstr[i].nightnum = i+1
  mstr[i].nstars = nmatch
  mstr[i].rms = rms
  mstr[i].brtrms = brtrms
  mstr[i].zpterm = median(mfitstr.zpterm)
  mstr[i].zptermsig = mad(mfitstr.zpterm)
  mstr[i].colterm = median(mfitstr.colterm)
  mstr[i].coltermsig = mad(mfitstr.colterm)
  mstr[i].amterm = median(mfitstr.amterm)
  mstr[i].amtermsig = mad(mfitstr.amterm)

  ; Convert MJD to YYYYMMDD date
  jd = umjds[i]+2400000.5d0-0.5+obs.tz/24.  
  caldat,jd,month,day,year,hour,min,sec
  date = strtrim(year,2)+string(month,format='(I02)')+string(day,format='(I02)')
  mstr[i].date = date

  ; Print summary info to screen
  print,i+1,umjds[i],date,nmatch,median(mfitstr.zpterm),median(mfitstr.colterm),mfitstr[0].amterm,rms,brtrms,medsig,format='(I5,I8,A10,I8,6F11.4)'
  ; Save structure
  push,fitstr,mfitstr

  ;stop

endfor
nfitstr = n_elements(fitstr)

rnd = randomu(seed,nfitstr)*0.6-0.3

plotc,fitstr.nightnum+rnd,fitstr.zpterm,fitstr.sig,ps=3,xr=[0,nmjds+1],yr=[-1,1],xs=1,ys=1,max=0.1,$
      xtit='Night Number',ytit='ZP term',tit='Zero-point term (color-coded by zptermerr)'

; Photometric night, based on zpterm scatter
bdmjd = where(mstr.brtrms gt median(mstr.brtrms)+2.5*mad(mstr.brtrms),nbdmjd)
print,strtrim(nbdmjd,2),' nights have high RMS'
if nbdmjd gt 0 then print,umjds[bdmjd]

; 8   56664  20140106, only 16 stars,  SO FEW STARS, DON'T TRUST ANY OF THE RESULTS
;      maybe determine single global zpterm for the night
; 9   56665  20140107, only 51 stars,  not much better
bdmjd = where(mstr.nstars lt 100,nbdmjd)
for i=0,nbdmjd-1 do begin

  ; Fit global zpterm, use colterm and amterm from neighbors

  stop

endfor


stop

; Unique CHIPS
ui = uniq(fitstr.chip,sort(fitstr.chip))
uchips = fitstr[ui].chip
nuchips = n_elements(uchips)

; Measure median and RMS per chip
chipmedzpterm = fltarr(nuchips)+99.99
chipzptermsig = fltarr(nuchips)+99.99
chipzptermerr = fltarr(nuchips)+99.99
chiprelzpterm = fltarr(nfitstr)+99.99
chipmedcolterm = fltarr(nuchips)+99.99
chipcoltermsig = fltarr(nuchips)+99.99
chipcoltermerr = fltarr(nuchips)+99.99
chiprelcolterm = fltarr(nfitstr)+99.99
chipmedrms = fltarr(nuchips)+99.99

for i=0,nuchips-1 do begin
  ind = where(fitstr.chip eq uchips[i],nind)
  if nind gt 0 then begin
    chipmedzpterm[i] = median([fitstr[ind].zpterm])
    chipzptermsig[i] = mad([fitstr[ind].zpterm])
    chipzptermerr[i] = chipzptermsig[i] / sqrt(nind)
    chiprelzpterm[ind] = fitstr[ind].zpterm - chipmedzpterm[i]
    chipmedcolterm[i] = median([fitstr[ind].colterm])
    chipcoltermsig[i] = mad([fitstr[ind].colterm])
    chipcoltermerr[i] = chipcoltermsig[i] / sqrt(nind)
    chiprelcolterm[ind] = fitstr[ind].colterm - chipmedcolterm[i]
    chipmedrms[i] = median([fitstr[ind].rms])
  endif
endfor

plotc,fitstr.chip+rnd,fitstr.zpterm,fitstr.mjd,ps=3,yr=[-0.5,0.5],xs=1
plotc,fitstr.chip+rnd,fitstr.zpterm,fitstr.zptermerr,ps=3,yr=[-0.5,0.5],xs=1,max=0.01
oplot,uchips,chipmedzpterm,ps=-1
plotc,fitstr.chip+rnd,fitstr.colterm,fitstr.mjd,ps=3,yr=[-0.3,0.1],xs=1
plotc,fitstr.chip+rnd,fitstr.colterm,fitstr.coltermerr,ps=3,yr=[-0.3,0.1],xs=1,max=0.01
oplot,uchips,chipmedcolterm,ps=-1


; Measure median relative chip zpterm/colterm per night
mjdchipmedzpterm = fltarr(nmjds)
mjdchipzptermsig = fltarr(nmjds)
mjdchipzptermerr = fltarr(nmjds)
mjdchipmedcolterm = fltarr(nmjds)
mjdchipcoltermsig = fltarr(nmjds)
mjdchipcoltermerr = fltarr(nmjds)
nightnum = lonarr(nfitstr)
for i=0,nmjds-1 do begin
  ind = where(fitstr.mjd eq umjds[i],nind)
  ; ZP term
  mjdchipmedzpterm[i] = median( chiprelzpterm[ind] )
  mjdchipzptermsig[i] = mad(chiprelzpterm[ind])
  mjdchipzptermerr[i] = mjdchipzptermsig[i] / sqrt(nind)
  ; Color term
  mjdchipmedcolterm[i] = median( chiprelcolterm[ind] )
  mjdchipcoltermsig[i] = mad(chiprelcolterm[ind])
  mjdchipcoltermerr[i] = mjdchipcoltermsig[i] / sqrt(nind)
  nightnum[ind] = i+1   ; night number
endfor

; Zero-point
plotc,fitstr.mjd+rnd,chiprelzpterm,fitstr.coltermerr,ps=3,yr=[-0.1,0.1],max=0.01,xtit='MJD',ytit='Relative ZP Term',tit='ZP relative to median per chip'
oplot,umjds,mjdchipmedzpterm,ps=1,sym=2
; Color term vs. night
;plotc,fitstr.mjd+rnd,chiprelcolterm,fitstr.coltermerr,ps=3,yr=[-0.02,0.02],max=0.01,xtit='MJD',ytit='Relative Color Term',tit='Color term relative to median per chip'
;oplot,umjds,mjdchipmedcolterm,ps=1,sym=2
z = fitstr.zpterm ; fitstr.coltermerr
;zr = minmax(z)
zr = [-0.5,0.5]
ztit = 'zpterm'
plotc,nightnum+rnd,chiprelcolterm,z,ps=8,sym=0.6,xs=1,ys=1,xr=[0,max(nightnum)+1],yr=[-0.06,0.06],min=zr[0],max=zr[1],$
      xtit='Night Number',ytit='Relative Color Term',tit='Color term relative to median per chip (color-coded by '+ztit+')'
;col = scale(fitstr.coltermerr,[min(fitstr.coltermerr),0.01],[50,250])<250
col = 50 > scale(z,zr,[50,250]) < 250
for i=0,nfitstr-1 do oplot,[0,0]+nightnum[i]+rnd[i],fitstr[i].coltermerr*[-1,1]+chiprelcolterm[i],co=col[i]
oplot,indgen(nmjds)+1,mjdchipmedcolterm,ps=1,sym=2,thick=2
oplot,[0,100],[-1,-1]*mad(mjdchipmedcolterm),linestyle=2
oplot,[0,100],[1,1]*mad(mjdchipmedcolterm),linestyle=2


; There definitively are mean color term variations for a few
; nights, it's about a ~3% variation
; COULD THESE BE NON-PHOTOMETRIC NIGHTS??


; Color term vs. chip
z = nightnum
zr = minmax(z)
ztit = 'nightnum'
plotc,fitstr.chip+rnd,chiprelcolterm,z,ps=8,sym=0.6,xs=1,ys=1,xr=[0,max(uchips)+1],yr=[-0.06,0.06],min=zr[0],max=zr[1],$
      xtit='Chip Number',ytit='Relative Color Term',tit='Color term relative to median per chip (color-coded by '+ztit+')'
col = 50 > scale(z,zr,[50,250]) < 250
for i=0,nfitstr-1 do oplot,[0,0]+fitstr[i].chip+rnd[i],fitstr[i].coltermerr*[-1,1]+chiprelcolterm[i],co=col[i]
oplot,indgen(nmjds)+1,chipmedcolterm,ps=1,sym=2,thick=2
oplot,[0,100],[0,0],linestyle=1
oplot,[0,100],[-1,-1]*mad(chiprelcolterm),linestyle=2
oplot,[0,100],[1,1]*mad(chiprelcolterm),linestyle=2
; I don't see any real systematics here


; Discrepant nights
bd = where(abs(mjdchipmedcolterm) gt 0.04,nbd)
print,strtrim(nbd,2),' nights with discrepant median color terms'
if nbd gt 0 then print,umjds[bd]

stop

; some chip-to-chip variations
; some night-to-night variations, especially between two years
;   could some of those be from using SDSS and the Smith standards?
;  maybe it was a non-photometric night???

; These groups of days have lower color terms then the ones after that
; 56369-56372 and 56509+56510
; March 2013 and Aug 2013
; you see this in g-band and r-band but not i-band or u-band
; z-band doesn't show a constant offset, but the chip-to-chip patterns
; seem to change

; 56677-56686 color-terms are a bit high, but these are the Jan 21-28, 2014
; days when only small-field Smith/DES fields were observed and so
; only a handful of chips have observations
; the one day in that range that has low RMS seems to have a similar
; color-term in the r-band.

plotc,mchfitstr.chip+rnd,mchfitstr.rms,mchfitstr.mjd,ps=3,yr=[-2,2],xs=1  


; airmass terms
ploterror,mfitstr.mjd-min(mfitstr.mjd),mfitstr.coef[0],mfitstr.mjd*0,mfitstr.coeferr[1],ps=1,yr=[-1,1]

; 56678 and 56680 have less than 50 stars
; 56677, 56678, 56680, 56683-56686, less than 600 stars, not good fits
; maybe use the "standard" value for these


rndmjd = randomu(seed,numjdchip)*0.8-0.4
plotc,mchfitstr.mjd+rndmjd,mchfitstr.coef[0],mchfitstr.chip,ps=3,yr=[-1,1],xs=1
plotc,mchfitstr.mjd+rndmjd,mchfitstr.coef[1],mchfitstr.chip,ps=3,yr=[-0.3,0.1],xs=1
; plot using unique night numbers



; Measure median and RMS per chip
chipmed = fltarr(2,nuchips)+99.99
chiprms = fltarr(nuchips)+99.99
for i=0,nuchips-1 do begin
  ind = where(mchfitstr2.chip eq uchips[i],nind)
  if nind gt 0 then begin
    chipmed[0,i] = median([mchfitstr2[ind].coef[0]])
    chipmed[1,i] = median([mchfitstr2[ind].coef[1]])
    chiprms[i] = median([mchfitstr2[ind].rms])
  endif
endfor


plothist,mfitstr.medrms,bin=0.04,xr=[0,2]

bdmjd = where(mfitstr.rms gt median(mfitstr.rms)+2.5*mad(mfitstr.rms),nbdmjd)
print,strtrim(nbdmjd,2),' nights have high RMS'
if nbdmjd gt 0 then print,umjds[bdmjd]

; There are definitely some outlier stars, should remove them
plot,arr.mag,resid3,ps=3,yr=[-1,1]

; Compare the two RMS values, should improve slightly with airmass
; term removed
plot,mchfitstr.rms,ps=1,yr=[0,1]
oplot,mchfitstr2.rms,ps=4,co=250

stop

end
