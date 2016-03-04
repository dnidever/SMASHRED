;+
;
; SMASHRED_PHOTCALIB
;
; This program calibrates photometry for a field using traditional
; photometric transformation equations, ubercal and additional
; photometric "anchor" data (i.e. 0.9m data).
;
; INPUTS:
;
; OUTPUTS:
;
; USAGE:
;  IDL>smashred_photcalib,fstr,chstr,allsrc,allobj,
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

; Get unique filters
ui = uniq(fstr.filter,sort(fstr.filter))
ufilter = fstr[ui].filter
nufilter = n_elements(ufilter)

; Iterate until we converge
;--------------------------
flag = 0
count = 0
WHILE (flag eq 0) do begin

  ; Filter loop
  For f=0,nufilter-1 do begin

    ifilter = ufilter[i]
    chfiltind = where(chstr.filter eq ifilter,nchfiltind)
    chfiltstr = chstr[chfiltind]

    print,'--- FILTER = ',ifilter,' ',strtrim(nchfiltind,2),' nchips'

    ; Step 1. Measure photometric offsets between chip pairs
    ;SMASHRED_MEASURE_MAGOFFSET,chfiltstr,overlapstr
    SMASHRED_MEASURE_MAGOFFSET,chfiltstr,allsrc,overlapstr,/usecalib

    ; Step 2. determine relative magnitude offsets per chip using ubercal 

    ; Solve the UBERCAL problem
    PERFORM_UBERCAL_CALIB,overlapstr,fmagoff

    ; Apply these magnitude offsets to the data
    for k=0,nchfiltind-1 do begin
      (*chfiltstr[k].data).mag += fmagoff[k]
      chfiltstr[k].magoffset = fmagoff[k]
    endfor


    ; Step 3. Set zeropoint per band with "photometric" data and/or anchor data

    ; Step 4. Calculate average photometry per object per filter
    SMASHRED_AVERAGEPHOT,fstr,chstr,allsrc,allobj,/usecalib

    ; Step 5. calibrate instrumental photometry with trans eqns. and ubercal mag offsets

    stop

  Endfor ; filter loop

  ; Have we converged?

  count++

ENDWHILE


stop

end
