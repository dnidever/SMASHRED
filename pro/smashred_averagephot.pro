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
;  /deeponly  Only use deep exposures.
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

pro smashred_averagephot,fstr,chstr,allsrc,allobj,usecalib=usecalib,deeponly=deeponly,error=error,silent=silent

; Not enough inputs
if n_elements(fstr) eq 0 or n_elements(chstr) eq 0 or n_elements(allsrc) eq 0 or n_elements(allobj) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_averagephot,fstr,chstr,allstr,allobj,usecalib=usecalib,deeponly=deeponly,error=error,silent=silent'
  return
endif

lallobjtags = strlowcase(tag_names(allobj))
nallobj = n_elements(allobj)

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
    error = '/USECALIB set but NO CMAG/CERR columns in ALLSRC'
    if not keyword_set(silent) then print,error
    return
  endif
endif

; Combine photometry from same filter
if not keyword_set(silent) then if keyword_set(deeponly) then print,'Combining all of the DEEP photometry' else $
    print,'Combining all of the photometry'

for i=0,n_elements(ufilter)-1 do begin

  ; Using all exposures
  if not keyword_set(deeponly) then begin
    ; Number of exposures for this filter
    filtind = where(fstr.filter eq ufilter[i],nfiltind)
    ; Chips for this filter
    chind = where(chstr.filter eq ufilter[i],nchind)

  ; Only deep exposures
  endif else begin
    ; Number of exposures for this filter
    filtind = where(fstr.filter eq ufilter[i] and fstr.exptime gt 100,nfiltind)
    ; Chips for this filter
    chind = where(chstr.filter eq ufilter[i] and chstr.exptime gt 100,nchind)
  endelse

  ; No exposures in this band, skip
  if nfiltind eq 0 then begin
    if keyword_set(deeponly) then cmt='deep ' else cmt=''
    print,'No '+cmt+'exposures in filter='+ufilter[i]
    goto,BOMB
  endif

  ; Indices for the magnitude and errors in ALLOBJ
  magind = where(lallobjtags eq ufilter[i])
  errind = where(lallobjtags eq ufilter[i]+'err')
  scatind = where(lallobjtags eq ufilter[i]+'scatter')
  nobsind = where(lallobjtags eq 'ndet'+ufilter[i])

  ; Only one exposure for this filter, copy
  if nfiltind eq 1 then begin

    ; All bad to start
    allobj.(magind) = 99.99
    allobj.(errind) = 9.99
    allobj.(scatind) = 99.99

    ; Now copy in the values, ALLSRC only had "good" detections
    for k=0,nchind-1 do begin
      ind = lindgen(chstr[chind[k]].nsrc)+chstr[chind[k]].allsrcindx
      if keyword_set(usecalib) then begin  ; calibrated phot
        allobj[allsrc[ind].cmbindx].(magind) = allsrc[ind].cmag
        allobj[allsrc[ind].cmbindx].(errind) = allsrc[ind].cerr
      endif else begin                ; instrumental phot
        allobj[allsrc[ind].cmbindx].(magind) = allsrc[ind].mag
        allobj[allsrc[ind].cmbindx].(errind) = allsrc[ind].err
      endelse
    endfor

  ; Multiple exposures for this filter to average
  endif else begin

    ; Loop through all of the chips and add up the flux, totalwt, etc.
    totalwt = dblarr(nallobj)
    totalfluxwt = dblarr(nallobj)
    for k=0,nchind-1 do begin
      ind = lindgen(chstr[chind[k]].nsrc)+chstr[chind[k]].allsrcindx
      if keyword_set(usecalib) then begin  ; calibrated phot
        totalwt[allsrc[ind].cmbindx] += 1.0d0/allsrc[ind].cerr^2
        totalfluxwt[allsrc[ind].cmbindx] += 2.5118864d^allsrc[ind].cmag * (1.0d0/allsrc[ind].cerr^2)
      endif else begin                ; instrumental phot 
        totalwt[allsrc[ind].cmbindx] += 1.0d0/allsrc[ind].err^2
        totalfluxwt[allsrc[ind].cmbindx] += 2.5118864d^allsrc[ind].mag * (1.0d0/allsrc[ind].err^2)
      endelse
    endfor
    newflux = totalfluxwt/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bdmag = where(finite(newmag) eq 0,nbdmag)
    if nbdmag gt 0 then begin
      newmag[bdmag] = 99.99
      newerr[bdmag] = 9.99
    endif

    ; Measure scatter, RMS
    ;  sqrt(mean(diff^2))
    totaldiff = dblarr(nallobj)
    numobs = lonarr(nallobj)
    for k=0,nchind-1 do begin
      ind = lindgen(chstr[chind[k]].nsrc)+chstr[chind[k]].allsrcindx
      if keyword_set(usecalib) then begin
        totaldiff[allsrc[ind].cmbindx] += (newmag[allsrc[ind].cmbindx] - allsrc[ind].cmag)^2
      endif else begin
        totaldiff[allsrc[ind].cmbindx] += (newmag[allsrc[ind].cmbindx] - allsrc[ind].mag)^2
      endelse
      numobs[allsrc[ind].cmbindx]++
    endfor
    newscatter = sqrt( totaldiff/(numobs>1) )
    if nbdmag gt 0 then newscatter[bdmag]=99.99

    ; Set scatter=99.99 for numobs=1
    oneobs = where(numobs eq 1,noneobs)
    if noneobs gt 0 then newscatter[oneobs]=99.99

    allobj.(magind) = newmag
    allobj.(errind) = newerr
    allobj.(scatind) = newscatter
    allobj.(nobsind) = numobs

  endelse  ; combine multiple exposures for this filter
  BOMB:
endfor ; unique filter loop

;stop

end
