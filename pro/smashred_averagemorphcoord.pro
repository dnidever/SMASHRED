;+
;
; SMASHRED_AVERAGEMORPH
;
; Calculate average morphological parameters (chi, sharp, flag, prob)
; using ALLSRC and ALLOBJ.
;
; INPUTS:
;  fstr    The structure with information for each exposure.
;  chstr   The structure with information for each chip.
;  allstr  The structure with information for each source detection.
;  allobj  The structure with information for each unique object.
;
; OUTPUTS:
;  The morphological parameters (chi, sharp, flag, prob) are updated
;  in ALLOBJ.
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>smashred_averagemorph,fstr,chstr,allsrc,allobj
;
; By D.Nidever  March 2016
;-

pro smashred_averagemorph,fstr,chstr,allsrc,allobj,error=error

; Not enough inputs
if n_elements(fstr) eq 0 or n_elements(chstr) eq 0 or n_elements(allsrc) eq 0 or n_elements(allobj) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_averagemorph,fstr,chstr,allstr,allobj'
  return
endif

; Get average chi, sharp, flag, prob
nallobj = n_elements(allobj)
totchi = fltarr(nallobj) & numchi = lon64arr(nallobj)
totsharp = fltarr(nallobj) & numsharp = lon64arr(nallobj)
totprob = fltarr(nallobj) & numprob = lon64arr(nallobj)
flag = intarr(nallobj) & numflag = lon64arr(nallobj)
nfstr = n_elements(fstr)
for i=0,nfstr-1 do begin
  gd = where(allobj.srcfindx[i] ge 0,ngd)
  ; CHI
  chi1 = fltarr(nallobj)+!values.f_nan
  chi1[gd] = allsrc[allobj[gd].srcfindx[i]].chi
  gdchi = where(finite(chi1) eq 1 and chi1 lt 1e5,ngdchi)
  if ngdchi gt 0 then begin
    totchi[gdchi] += chi1[gdchi]
    numchi[gdchi]++
  endif
  ; SHARP
  sharp1 = fltarr(nallobj)+!values.f_nan
  sharp1[gd] = allsrc[allobj[gd].srcfindx[i]].sharp
  gdsharp = where(finite(sharp1) eq 1 and sharp1 lt 1e5,ngdsharp)
  if ngdsharp gt 0 then begin
    totsharp[gdsharp] += sharp1[gdsharp]
    numsharp[gdsharp]++
  endif
  ; PROB
  prob1 = fltarr(nallobj)+!values.f_nan
  prob1[gd] = allsrc[allobj[gd].srcfindx[i]].prob
  gdprob = where(finite(prob1) eq 1 and prob1 lt 50 and prob1 gt -0.5,ngdprob)
  if ngdprob gt 0 then begin
    totprob[gdprob] += prob1[gdprob]
    numprob[gdprob]++
  endif
  ; FLAG,  logical OR across all detections
  ;    0s are essentially ignored
  flag1 = intarr(nallobj)-1
  flag1[gd] = allsrc[allobj[gd].srcfindx[i]].flag
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
gdsharp = where(numsharp gt 0,ngdsharp)
avgsharp = fltarr(nallobj)+99.99
if ngdsharp gt 0 then avgsharp[gdsharp]=totsharp[gdsharp]/numsharp[gdsharp]
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

;stop

end
