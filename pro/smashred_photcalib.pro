;+
;
; SMASHRED_PHOTCALIB
;
; This program calibrates photometry for a field using traditional
; photometric transformation equations, ubercal and additional
; photometric "anchor" data (i.e. 0.9m data).
;
; INPUTS:
;  fstr    The structure with information for each exposure.
;  chstr   The structure with information for each chip.
;  allsrc  The structure with information for each source detection.
;  allobj  The structure with information for each unique object.
;
; OUTPUTS:
;  ??
;  ??
;
; USAGE:
;  IDL>smashred_photcalib,fstr,chstr,allsrc,allobj
;
; D.Nidever  March 2016
;-

pro smashred_photcalib,fstr,chstr,allsrc,allobj

nfstr = n_elements(fstr)
nchstr = n_elements(chstr)
nallsrc = n_elements(allsrc)
nallobj = n_elements(allobj)

; Not enough inputs
if nfstr eq 0 or nchstr eq 0 or nallsrc eq 0 or nallobj eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_photcalib,fstr,chstr,allsrc,allobj'
  return
endif

; also need the list of photometric nights, transformation equations,
; 0.9m data

; HOW DOES THE 0.9M DATA FIT IN????  Can't be done just as an
; afterburner, since this will affect the mean mags which will in turn
; affect the color terms.  So needs to be a 


; Initialize the transformation equations
; Add the photometric transformation equation columns to CHSTR
if tag_exist(chstr,'photometric') eq 0 then begin
  newtags = ['photometric','band','colband','colsign','zpterm','zptermsig','amterm',$
             'amtermsig','colterm','coltermsig','amcolterm','amcoltermsig','colsqterm','colsqtermsig']
  temp = chstr & undefine,chstr
  add_tags,temp,newtags,['0B','""','""','-1',replicate('0.0d0',10)],chstr
  undefine,temp
endif
; Stuff in the transformation equation information
;  I got these from stdred_201405/n1.trans
chstr.photometric = 1
chstr.band = chstr.filter
nchstr = n_elements(chstr)
for i=0,nchstr-1 do begin
  case chstr[i].filter of
  'u': begin
    chstr[i].colband = 'g'
    chstr[i].colsign = 1
    chstr[i].zpterm = 1.5137
    chstr[i].amterm = 0.3787
    chstr[i].colterm = 0.0094 
    chstr[i].amcolterm = 0.0
    chstr[i].colsqterm = 0.0
  end
  'g': begin
    ;chstr[i].colband = 'u'
    ;chstr[i].colsign = 2
    chstr[i].colband = 'r'
    chstr[i].colsign = 1
    chstr[i].zpterm = -0.3136
    chstr[i].amterm = 0.1582
    chstr[i].colterm = -0.0566
    chstr[i].amcolterm = 0.0
    chstr[i].colsqterm = 0.0
  end
  'r': begin
    chstr[i].colband = 'i'
    chstr[i].colsign = 1
    chstr[i].zpterm = -0.4785
    chstr[i].amterm = 0.0758
    chstr[i].colterm = -0.1484
    chstr[i].amcolterm = 0.0
    chstr[i].colsqterm = 0.0
  end
  'i': begin
    chstr[i].colband = 'z'
    chstr[i].colsign = 1
    chstr[i].zpterm = -0.3545
    chstr[i].amterm = 0.0380
    chstr[i].colterm = -0.2933
    chstr[i].amcolterm = 0.0
    chstr[i].colsqterm = 0.0
  end
  'z': begin
    chstr[i].colband = 'i'
    chstr[i].colsign = 2
    chstr[i].zpterm = -0.0423
    chstr[i].amterm = 0.0498
    chstr[i].colterm = -0.0675
    chstr[i].amcolterm = 0.0
    chstr[i].colsqterm = 0.0
  end
  else: stop
  endcase
endfor
chstr.zptermsig = 0.001
chstr.amtermsig = 0.001
chstr.coltermsig = 0.001
chstr.amcoltermsig = 0.0
chstr.colsqtermsig = 0.0

; SHOULD I ALSO COMPUTE EXPOSURE-LEVEL APERTURE CORRECTIONS??
; MAYBE FIT A WEAK SPATIAL POLYNOMIAL??



; Get unique filters
ui = uniq(fstr.filter,sort(fstr.filter))
ufilter = fstr[ui].filter
nufilter = n_elements(ufilter)

; Iterate until we converge
;--------------------------
lastcmag = allsrc.cmag
doneflag = 0
niter = 0
last_zpterm = chstr.zpterm
WHILE (doneflag eq 0) do begin

  print & print,'Global field photometric calibration.  Iteration = ',strtrim(niter+1,2) & print

  ;----------------------------------------------------------
  ; Step 1. Regular photometric calibration with trans eqns.
  ;==========================================================
  ; "regular" trans calib, use allobj for the other band
  ; needed to construct the color, but use the magnitude
  ; from that exposure for the other band to construct
  ; the color (just like photcalib.pro does), otherwise
  ; I'd have to run averagemag many times
  ; This also averages mags for allobj to get the color terms
  print,'--- Step 1. Regular photometric calibration with trans eqns. ---'

  ; Do regular calibration with transformation equations
  ;  for non-photometric data only color-terms are applied
  ; this updates CMAG/CERR in ALLSRC
  SMASHRED_APPLY_PHOTTRANSEQN,fstr,chstr,allsrc,allobj


  ;----------------------------------------------------------
  ; Step 2. Ubercal measurement and offsets to individual
  ;           exposure zeropoints
  ;==========================================================
  ; ubercal measurement and application, with filter loop
  ; add these offsets to the zpterm terms for each chip
  ; in the TRANS structure so they'll improve the photometry
  ; calculated in step 1 the next time around
  ; use trans.photometric to help anchor the offsets
  print,'--- Step 2. Ubercal measurement and application ---'

  ; Filter loop
  For f=0,nufilter-1 do begin
    ifilter = ufilter[f]
    chfiltind = where(chstr.filter eq ifilter,nchfiltind)
    chfiltstr = chstr[chfiltind]
    print,'- ',strtrim(f+1,2),' FILTER = ',ifilter,' ',strtrim(nchfiltind,2),' nchips'

    ; Step 2a. Measure photometric offsets between chip pairs
    print,'2a. Measuring relative magnitude offsets between chip pairs'
    SMASHRED_MEASURE_MAGOFFSET,chfiltstr,allsrc,overlapstr,/usecalib

    ; Step 2b. Determine relative magnitude offsets per chip using ubercal 
    print,'2b. Determine relative magnitude offsets per chip using ubercal'
    SMASHRED_SOLVE_UBERCAL,overlapstr,ubercalstr

    ; Step 2c. Set absolute zeropoint of magnitute offsets with
    ;            photometric data and/or anchor points
    print,'2c. Set absolute zeropoint with photometric data and achor points'
    SMASHRED_SET_ZEROPOINTS,chfiltstr,ubercalstr

    ; Stuff it back in
    chstr[chfiltind] = chfiltstr
  Endfor  ; filter loop

  ; BRIGHTER-FATTER CORRECTION??
  ; maybe make it a (instrumental) magnitude-dependent correction, bfcorr
  ; compare psf to aperture photometry to figure it out
  ; BUT CORRECT IT TO WHAT???  what's the fiducial value.
  ; it will likely depend on the seeing, so we might need to figure it
  ; out on an exposure-by-exposure level

  ; Check for convergence
  ;  check changes in ALLSRC CMAG values, changes in CHSTR ZPTERM and iterations
  maxiter = 50
  maxdiff_thresh = 0.0001
  maxzpterm_thresh = 0.0001
  diffcmag = abs(allsrc.cmag-lastcmag)
  diffzpterm = abs(chstr.zpterm-last_zpterm)
  if (max(diffcmag) lt maxdiff_thresh and max(diffzpterm) lt maxzpterm_thresh) or niter gt maxiter then doneflag=1
  last_zpterm = chstr.zpterm  ; save for next time
  lastcmag = allsrc.cmag

  niter++

  stop

ENDWHILE


stop

end
