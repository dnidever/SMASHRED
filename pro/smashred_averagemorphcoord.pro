;+
;
; SMASHRED_AVERAGEMORPHCOORD
;
; Calculate average morphological parameters (chi, sharp, flag, prob)
; and RA/DEC coordinates using ALLSRC and ALLOBJ.
;
; INPUTS:
;  fstr    The structure with information for each exposure.
;  chstr   The structure with information for each chip.
;  allsrc  The structure with information for each source detection.
;  allobj  The structure with information for each unique object.
;  /deeponly   Only use the deep exposures.
;  /shortonly  Only use the short exposures.
;  /alslowsnrcut   Downweight ALLSTAR S/N<5 sharp (short exposures)
;                    because are "bad".  Will keep the sharp if it
;                    is the only detection for that object.
;
; OUTPUTS:
;  The morphological parameters (chi, sharp, flag, prob) and
;  coordinate parameters (ra, dec, rascatter, decscatter) are
;  updated in ALLOBJ.
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>smashred_averagemorphcoord,fstr,chstr,allsrc,allobj
;
; By D.Nidever  March 2016
;-

pro smashred_averagemorphcoord,fstr,chstr,allsrc,allobj,error=error,alslowsnrcut=alslowsnrcut,deeponly=deeponly,$
                               shortonly=shortonly

; Not enough inputs
if n_elements(fstr) eq 0 or n_elements(chstr) eq 0 or n_elements(allsrc) eq 0 or n_elements(allobj) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_averagemorphcoord,fstr,chstr,allsrc,allobj,alslowsnrcut=alslowsnrcut,deeponly=deeponly,shortonly=shortonly'
  return
endif

; Get average chi, sharp, flag, prob
nallobj = n_elements(allobj)
totchi = fltarr(nallobj) & numchi = lon64arr(nallobj)
totsharp = fltarr(nallobj) & totwtsharp = fltarr(nallobj)
totprob = fltarr(nallobj) & numprob = lon64arr(nallobj)
flag = intarr(nallobj) & numflag = lon64arr(nallobj)
; Exposures to use
nfstr = n_elements(fstr)
; All exposures
if not keyword_set(deeponly) and not keyword_set(shortonly) then begin
  nexpind = nfstr
  expind = lindgen(nexpind)
endif
;  Only deep exposures
if keyword_set(deeponly) then expind = where(fstr.exptime gt 100,nexpind)
;  Only short exposures
if keyword_set(deeponly) then expind = where(fstr.exptime lt 100,nexpind)

; Exposure loop
for i=0,nexpind-1 do begin
  gd = where(allobj.srcfindx[expind[i]] ge 0,ngd)
  ; CHI
  chi1 = fltarr(nallobj)+!values.f_nan
  chi1[gd] = allsrc[allobj[gd].srcfindx[expind[i]]].chi
  gdchi = where(finite(chi1) eq 1 and chi1 lt 1e5,ngdchi)
  if ngdchi gt 0 then begin
    totchi[gdchi] += chi1[gdchi]
    numchi[gdchi]++
  endif
  ;; SHARP
  ;; use "weights" with wt=1 for normal values
  ;; and wt=0.0001 for short, allstar and S/N<5
  ;; Using weights will still use the low S/N
  ;; sharp value if it's the ONLY detection, which
  ;  is what I want.
  sharp1 = fltarr(nallobj)+!values.f_nan
  sharp1[gd] = allsrc[allobj[gd].srcfindx[expind[i]]].sharp
  gdsharp = where(finite(sharp1) eq 1 and sharp1 lt 1e5,ngdsharp)
  if ngdsharp gt 0 then begin
    wtsharp1 = fltarr(nallobj)
    wtsharp1[gdsharp] = 1.0
    ; Set wt to 1e-4 for ALLSTAR short exposures with S/N<5
    if fstr[expind[i]].exptime lt 100 and keyword_set(alslowsnrcut) then begin
      snr = fltarr(nallobj)
      snr[gd] = 1.087/allsrc[allobj[gd].srcfindx[expind[i]]].err
      lowsnrind = where(finite(sharp1) eq 1 and sharp1 lt 1e5 and snr lt 5,nlowsnrind)
      if nlowsnrind gt 0 then wtsharp1[lowsnrind]=1e-4
    endif
    totsharp[gdsharp] += wtsharp1[gdsharp]*sharp1[gdsharp]
    totwtsharp[gdsharp] += wtsharp1[gdsharp]
  endif
  ; PROB
  prob1 = fltarr(nallobj)+!values.f_nan
  prob1[gd] = allsrc[allobj[gd].srcfindx[expind[i]]].prob
  gdprob = where(finite(prob1) eq 1 and prob1 lt 50 and prob1 gt -0.5,ngdprob)
  if ngdprob gt 0 then begin
    totprob[gdprob] += prob1[gdprob]
    numprob[gdprob]++
  endif
  ; FLAG,  logical OR across all detections
  ;    0s are essentially ignored
  flag1 = intarr(nallobj)-1
  flag1[gd] = allsrc[allobj[gd].srcfindx[expind[i]]].flag
  gdflag = where(flag1 ne -1,ngdflag)
  if ngdflag gt 0 then begin
    flag[gdflag] OR= flag1[gdflag]
    numflag[gdflag]++
  endif
endfor
; Make average CHI
gdchi = where(numchi gt 0,ngdchi)
avgchi = fltarr(nallobj)+99.99
if ngdchi gt 0 then avgchi[gdchi]=totchi[gdchi]/numchi[gdchi]
allobj.chi = avgchi
; Make average SHARP
gdsharp = where(totwtsharp gt 0,ngdsharp)
avgsharp = fltarr(nallobj)+99.99
if ngdsharp gt 0 then avgsharp[gdsharp]=totsharp[gdsharp]/totwtsharp[gdsharp]
allobj.sharp = avgsharp
; Make average PROB
gdprob = where(numprob gt 0,ngdprob)
avgprob = fltarr(nallobj)+99.99
if ngdprob gt 0 then avgprob[gdprob]=totprob[gdprob]/numprob[gdprob]
allobj.prob = avgprob
; Stuff FLAG in
bdflag = where(numflag eq 0,nbdflag)
if nbdflag gt 0 then flag[bdflag]=-1
allobj.flag = flag


; Measure astrometric median scatter
;-----------------------------------

; ra/dec scatter
totalra = dblarr(nallobj)
totaldec = dblarr(nallobj)
numobs = lonarr(nallobj)
nchstr = n_elements(chstr)
for k=0,nchstr-1 do begin
  if chstr[k].nsrc gt 0 then begin
    ind = lindgen(chstr[k].nsrc)+chstr[k].allsrcindx
    totalra[allsrc[ind].cmbindx] += allsrc[ind].ra
    totaldec[allsrc[ind].cmbindx] += allsrc[ind].dec
    numobs[allsrc[ind].cmbindx]++
  endif
endfor
newra = totalra/(numobs>1)
newdec = totaldec/(numobs>1)
bd = where(numobs eq 0,nbd)
if nbd gt 0 then newra[bd]=999999.0
if nbd gt 0 then newdec[bd]=999999.0

; measure scatter, RMS
;  sqrt(mean(diff^2))
totalradiff = dblarr(nallobj)
totaldecdiff = dblarr(nallobj)
for k=0,nchstr-1 do begin
  if chstr[k].nsrc gt 0 then begin
    ind = lindgen(chstr[k].nsrc)+chstr[k].allsrcindx
    totalradiff[allsrc[ind].cmbindx] += (newra[allsrc[ind].cmbindx] - allsrc[ind].ra)^2
    totaldecdiff[allsrc[ind].cmbindx] += (newdec[allsrc[ind].cmbindx] - allsrc[ind].dec)^2
  endif
endfor
newrascatter = sqrt( totalradiff/(numobs>1) ) * 3600 * cos(newdec/!radeg)
newdecscatter = sqrt( totaldecdiff/(numobs>1) ) * 3600
if nbd gt 0 then newrascatter[bd]=99.99
if nbd gt 0 then newdecscatter[bd]=99.99

; Set scatter=99.99 for numobs=1
oneobs = where(numobs eq 1,noneobs)
if noneobs gt 0 then newrascatter[oneobs]=99.99
if noneobs gt 0 then newdecscatter[oneobs]=99.99

allobj.ra = newra
allobj.dec = newdec
allobj.rascatter = newrascatter
allobj.decscatter = newdecscatter

; Calculate uncertainty in average RA/DEC
;  add uncertainties in quadrature
totalraerr = dblarr(nallobj)
totaldecerr = dblarr(nallobj)
for k=0,nchstr-1 do begin
  if chstr[k].nsrc gt 0 then begin
    ind = lindgen(chstr[k].nsrc)+chstr[k].allsrcindx
    totalraerr[allsrc[ind].cmbindx] += allsrc[ind].raerr^2
    totaldecerr[allsrc[ind].cmbindx] += allsrc[ind].decerr^2
  endif
endfor
; take sqrt to finish the quadrature and 
; divide by N since we want the
; uncertainty of the average
raerr = sqrt(totalraerr)/(numobs>1)
decerr = sqrt(totaldecerr)/(numobs>1)
bd = where(numobs eq 0,nbd)
if nbd gt 0 then raerr[bd]=999999.0
if nbd gt 0 then decerr[bd]=999999.0
allobj.raerr = raerr
allobj.decerr = decerr

;stop

end
