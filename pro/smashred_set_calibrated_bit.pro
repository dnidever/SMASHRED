;+
;
; SMASHRED_SET_CALIBRATED_BIT
;
; This determines which chips have been calibrated and which ones
; haven't.
;
; INPUTS:
;  chstr       Structure with information on each individual chip.
;  overlapstr  The structure containing the relative magnitude offsets
;                of chip pairs.
;
; OUTPUTS:
;  The CALIBRATED column is updated in the CHSTR structure for
;  chips that are calibrated.
;
; USAGE:
;  IDL>smashred_set_calibrated_bit,chstr,overlapstr
;
; By D. Nidever  August 2016
;-

pro smashred_set_calibrated_bit,chstr,overlapstr

; Not enough inputs
if n_elements(chstr) eq 0 or n_elements(overlapstr) eq 0 then begin
  print,'Syntax - smashred_set_calibrated_bit,chstr,overlapstr'
  return
endif

chstr.calibrated = 0   ; all uncalibrated to start with

data = overlapstr.data
allchipid = strtrim(chstr.expnum,2)+'-'+strtrim(chstr.chip,2)  ; for matching

; Get calibrated chips
;   ZPCALIBTYPE = 1, 2 and 3 are "real" zeropoint calibrations
;     0 - uncalibrated, 4 - GAIA calibrated
gdphot = where(chstr.zpcalibflag ge 1 and chstr.zpcalibflag le 3,ngdphot,comp=bdphot,ncomp=nbdphot)
;gdphot = where(chstr.photometric ge 1 and chstr.badsoln eq 0,ngdphot,comp=bdphot,ncomp=nbdphot)
if ngdphot gt 0 then begin

  ; Set CALIBRATED for well-calibrated chips
  chstr[gdphot].calibrated = 1

  ; Some are uncalibrated that could get calibrated by ubercal
  if nbdphot gt 0 then begin

    ; Loop until done
    doneflag = 0
    niter = 0
    last_calibrated = chstr.calibrated
    WHILE (doneflag eq 0) do begin

      ; Get the calibrated chips
      MATCH,chstr.calibrated,1,ind1,ind2,/sort,count=nmatch
      ; Loop through the calibrated chips and find overlap chips
      For i=0,nmatch-1 do begin
        ; Find overlap chips with good crossmatches
        ind = overlapstr.index[overlapstr.ind0[ind1[i]]:overlapstr.ind1[ind1[i]]]
        revind = overlapstr.revindex[overlapstr.ind0[ind1[i]]:overlapstr.ind1[ind1[i]]]
        gd = where(data[ind].nmatch gt 3 and data[ind].magofferr lt 0.05,ngd)
        ; Some good overlaps
        if ngd gt 0 then chstr[data[ind[gd]].chstrindex2].calibrated = 1
      Endfor  ; calibrated chip loop

      ; Are we done?
      maxiter = 1000
      diff_calibrated = last_calibrated - chstr.calibrated
      if total(abs(diff_calibrated)) eq 0 or niter gt maxiter then doneflag=1
      last_calibrated = chstr.calibrated

      ; Increment the counter
      niter++
    ENDWHILE

  endif  ; some uncalibrated chips

; No calibrated chips at all
endif else begin
  print,'   No calibrated chips'
endelse

; Print out the results
calind = where(chstr.calibrated eq 1,ncalind)
print,' ',strtrim(ncalind,2),' of ',strtrim(n_elements(chstr),2),' chips calibrated'

;stop

end
