pro smashred_set_calibrated_bit,chstr,overlapstr

; This determines which chips have been calibrated and which ones
; haven't

chstr.calibrated = 0   ; all uncalibrated to start with

data = overlapstr.data
allchipid = strtrim(chstr.expnum,2)+'-'+strtrim(chstr.chip,2)  ; for matching

; Get calibrated chips
gdphot = where(chstr.photometric eq 1 and chstr.badsoln eq 0,ngdphot,comp=bdphot,ncomp=nbdphot)
if ngdphot gt 0 then begin

  ; Set CALIBRATED for photometric chips
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
print,'   ',strtim(ncalind,2),' of ',strtrim(n_elements(chstr),2),' calibrated'

;stop

end
