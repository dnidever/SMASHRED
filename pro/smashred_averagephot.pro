;+
;
; SMASHRED_AVERAGEPHOT
;
; Calculate average photometry for each unique object
; and filter using ALLSRC and ALLOBJT
;
; INPUTS:
;  fstr       The structure with information for each exposure.
;  chstr      The structure with information for each chip.
;  allstr     The structure with information for each source detection.
;  allobj     The structure with information for each unique object.
;  /usecalib  Average the calibrated photometry (CMAG/CERR).  The default
;               is to use the instrumental photometry (MAG/ERR).
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  The photometric magnitude and error columns will be updated in ALLOBJ.
;  =error     The error message if one occurred.
;
; USAGE:
;  IDL>smashred_averagephot,fstr,chstr,allsrc,allobj
;
; By D.Nidever  March 2016
;-

pro smashred_averagephot,fstr,chstr,allsrc,allobj,usecalib=usecalib,error=error,silent=silent

; Not enough inputs
if n_elements(fstr) eq 0 or n_elements(chstr) eq 0 or n_elements(allsrc) eq 0 or n_elements(allobj) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_averagephot,fstr,chstr,allstr,allobj,usecalib=usecalib,error=error,silent=silent'
  return
endif

; Get unique filters
ui = uniq(fstr.filter,sort(fstr.filter))
ufilter = fstr[ui].filter
nufilter = n_elements(ufilter)

; Use calibrated magnitudes, check that we have them
if keyword_set(usecalib) then begin
  srctags = tag_names(allsrc)
  dum = where(srctags eq 'CMAG',ncmag)
  dum = where(srctags eq 'CERR',ncerr)
  if ncmag eq 0 or ncerr eq 0 then begin
    error = '/USECALIB set but NO CMAG/CERR columsn in ALLSRC'
    if not keyword_set(silent) then print,error
    return
  endif
endif

; Combine photometry from same filter
print,'Combing all of the photometry'
for i=0,n_elements(ufilter)-1 do begin
  filtind = where(fstr.filter eq ufilter[i],nfiltind)
  print,ufilter[i]

  ; Indices for the magnitude and errors in ALLOBJ
  magind = where(lallobjtags eq ufilter[i])
  errind = where(lallobjtags eq ufilter[i]+'err')

  ; Only one exposure for this filter, copy
  if nfiltind eq 1 then begin
    ; Get stars that have detections in this frame
    gd = where(allobj.srcfindx[filtind] ge 0,ngd)
    if keyword_set(usecalib) then begin  ; calibrated phot
      allobj[gd].(magind) = allsrc[allobj[gd].srcfindx[filtind[0]]].cmag
      allobj[gd].(errind) = allsrc[allobj[gd].srcfindx[filtind[0]]].cerr
    endif else begin                ; instrumental phot
      allobj[gd].(magind) = allsrc[allobj[gd].srcfindx[filtind[0]]].mag
      allobj[gd].(errind) = allsrc[allobj[gd].srcfindx[filtind[0]]].err
    endelse
    bd = where(allobj.(magind) gt 50,nbd)
    if nbd gt 0 then begin
      allobj[bd].(magind) = 99.99
      allobj[bd].(errind) = 9.99
    endif

  ; Multiple exposures for this filter to average
  endif else begin

    mag = fltarr(n_elements(allobj),nfiltind)+99.99
    err = mag*0+9.99
    for k=0,nfiltind-1 do begin
      gd = where(allobj.srcfindx[filtind[k]] ge 0,ngd)
      if keyword_set(usecalib) then begin  ; calibrated phot
        mag[gd,k] = allsrc[allobj[gd].srcfindx[filtind[k]]].cmag
        err[gd,k] = allsrc[allobj[gd].srcfindx[filtind[k]]].cerr
      endif else begin                ; instrumental phot 
        mag[gd,k] = allsrc[allobj[gd].srcfindx[filtind[k]]].mag
        err[gd,k] = allsrc[allobj[gd].srcfindx[filtind[k]]].err
      endelse
    endfor
    ; Ignore non-detections
    bd = where(mag eq 0.0 or mag gt 50,nbd)
    if nbd gt 0 then begin
      mag[bd] = !values.f_nan
      err[bd] = !values.f_nan
    endif

    ; copied from phot_overlap.pro
    flux = 2.511864d^mag
    wt = 1.0d0/err^2
    totalwt = total(wt,2,/nan)
    totalflux = total(flux*wt,2,/nan)
    totalerr = total((err^2)*wt,2,/nan) 
    newflux = totalflux/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bd = where(finite(newmag) eq 0,nbd)
    if nbd gt 0 then begin
      newmag[bd] = 99.99
      newerr[bd] = 9.99
    endif

    allobj.(magind) = newmag
    allobj.(errind) = newerr
  endelse  ; combine multiple exposures for this filter
endfor ; unique filter loop

;stop

end
