;+
;
; SMASHRED_APPLY_PHOTTRANSEQN
;
; This program applies photometric transformation equations to data.
; Very similar to photcalib.pro
;
; INPUTS:
;  fstr     The structure giving informaton on each exposure.
;  chstr    The structure giving information on each chip exposure
;             including the transformation equations.
;  allsrc   The structure with information on each detected source.
;  allobj   The structure with information on each unique object.
;
; OUTPUTS:
;  The CMAG and CERR values in ALLSRC are updated.
;
; USAGE:
;  IDL>smashred_apply_phottranseqn,fstr,chstr,allsrc,allobj
;
; By D. Nidever  April 2016
;-

function smashred_apply_phottranseqn_simplestar,inmag,clr,t,mask
;   V = mV - v1 - v2 * XV - v3 * (B-V) - v4 * XV*(B-V) - v5 * (B-V) * (B-V)
;         +(aperture correction) + (time correction) 
;   The aperture corrections need to be POSITIVE
outmag = inmag - t.zpterm - t.amterm*t.airmass - t.colterm*clr - t.amcolterm*t.airmass*clr - t.colsqterm*clr*clr - t.apcor + 2.5*alog10(t.exptime)
if n_elements(mask) gt 0 then outmag = outmag*mask + (mask-1)*99.99
return,outmag
end

;---------------------------------------------

function smashred_apply_phottranseqn_simplerr,inerr,clr,clrerr,t,mask
; sigma(V)=sqrt[ sigma(mV)^2 + + sigma(v1)^2 + (XV*sigma(v2))^2+((B-V)*sigma(v3))^2 
;                + (XV*(B-V)*sigma(v4))^2 + (B-V*sigma(v5))^2]
;                +(v3+v4*XV*+v5*2*(B-V))^2*sigma(B-V)^2  
outerr = sqrt( inerr^2 + t.zptermsig^2 + (t.airmass*t.amtermsig)^2  + $
               (clr*t.coltermsig)^2 + (t.airmass*clr*t.amcoltermsig)^2 + $
               (clr*clr*t.colsqtermsig)^2 + $
               (t.colterm+t.amcolterm*t.airmass+2*t.colsqterm*clr)^2*clrerr^2 )
if n_elements(mask) gt 0 then outerr = outerr*mask + (mask-1)*9.99
return,outerr
end

;---------------------------------------------

pro smashred_apply_phottranseqn,fstr,chstr,allsrc,allobj

; Not enough inputs
if n_elements(fstr) eq 0 or n_elements(chstr) eq 0 or n_elements(allsrc) eq 0 or n_elements(allobj) eq 0 then begin
  print,'Syntax - smashred_apply_phottranseqn,fstr,chstr,allsrc,allobj'
  return
endif

   
;=====================================================================  
;
;   The heart of the program.  This iteratively solves for each star.
;   It applies the transformation equation assuming a color of zero.  It then
;   averages the passbands that are common, solves for the color and resolves for
;   each magnitude given the new color, gradually iterating until
;   convergence.
;   
;   The solved magnitude is in the form: (V in the example)
;    V = mV - v1 - v2 * XV - v3 * (B-V) - v4 * XV*(B-V)
;           - v5 * (B-V) * (B-V) + (aperture correction)+(time correction)  
;
;  It sends back to the code outstar, an array of average values in
;  each passband and tempstar, an array of the individual solved magnitudes
;
;=====================================================================

nchstr = n_elements(chstr)
allobjtags = tag_names(allobj)

;##################################
; Iteration loop
niter = 0
converge = 0
WHILE (converge eq 0) do begin

  chstr_maxdiff = fltarr(nchstr)+99.0  ; initalize the difference for the outer while loop

  ; Chip loop
  For i=0,nchstr-1 do begin

    ; Other band for the color
    colband = chstr[i].colband
    colsign = chstr[i].colsign
    nsrc = chstr[i].nsrc

    ; The input magnitudes and errors
    allsrcind = lindgen(chstr[i].nsrc)+chstr[i].allsrcindx
    inmag = allsrc[allsrcind].mag
    inerr = allsrc[allsrcind].err

    ; Initializing some arrays
    clr = fltarr(nsrc)*0.
    clrerr = fltarr(nsrc)*0.
    laststar = fltarr(nsrc)*0.

    ; Get the color information from ALLOBJ
    colmagind = where(allobjtags eq strupcase(colband),ncolmagind)
    colerrind = where(allobjtags eq strupcase(colband+'err'),ncolerrind)
    colmag = allobj[allsrc[allsrcind].cmbindx].(colmagind[0])
    colmagerr = allobj[allsrc[allsrcind].cmbindx].(colerrind[0])

    ; MAGMASK, COLMASK, good photometry for mag and col
    magmask = (inmag lt 50)
    colmask = (inmag lt 50 and colmag lt 50)

    ;#############################
    ; Initial guess with zero color
    ; Calculate transformed mags and errors
    tempmag = SMASHRED_APPLY_PHOTTRANSEQN_SIMPLESTAR(inmag,clr,chstr[i],magmask)
    ;temperr = SMASHRED_APPLY_PHOTTRANSEQN_SIMPLERR(inerr,clr,clrerr,chstr[i],magmask)

    ; Inner while loop, converge with the color band we have
    lastmag = tempmag
    doneflag = 0
    niter2 = 0
    While (doneflag eq 0) do begin

      ; Make the color
      if colsign eq 1 then begin
        clr = tempmag - colmag      
      endif else begin
        clr = colmag - tempmag
      endelse
      clrerr = sqrt( inerr^2 + colmagerr^2 )
      ; Stars with bad color, set clr=0 and clerr=0
      clr *= colmask
      clrerr *= colmask

      ; Calculate transformed mags again
      tempmag = SMASHRED_APPLY_PHOTTRANSEQN_SIMPLESTAR(inmag,clr,chstr[i],magmask)

      ; Are we done?
      absdiffmag = abs(tempmag-lastmag)
      if max(absdiffmag) lt 0.0001 or niter2 gt 10 then doneflag=1
      lastmag = tempmag   ; save the last values

      niter2++
    Endwhile

    ; Calculate final errors
    temperr = SMASHRED_APPLY_PHOTTRANSEQN_SIMPLERR(inerr,clr,clrerr,chstr[i],magmask)

    ; Calculate the differences
    chstr_maxdiff[i] = max(abs(tempmag-allsrc[allsrcind].cmag))

    ; Stuff final values into CMAG/CERR columns of ALLSRC
    allsrc[allsrcind].cmag = tempmag
    allsrc[allsrcind].cerr = temperr

  Endfor  ; chstr loop

  ; Now calculate the average mags for ALLOBJ
  SMASHRED_AVERAGEPHOT,fstr,chstr,allsrc,allobj,/usecalib,/silent

  ; Check for convergence
  ;  iterate until all changes are <0.0001 or niter = 50
  maxiter = 50
  maxdiff_thresh = 0.0001
  if max(chstr_maxdiff) lt maxdiff_thresh or niter gt maxiter then converge=1

  if not keyword_set(silent) then $
    print,'Iter = ',strtrim(niter+1,2),'  max diff = ',strtrim(max(chstr_maxdiff),2),' mag'

  ; Increment
  niter++
ENDWHILE


;stop

end
