;+
;
; SMASHRED_SET_ZEROPOINTS
;
; This determines the absolute zeropint of the magnitude offsets
; using photometric and/or anchor data and applies these to the
; chip-level transformation equations (zeropoints) in CHSTR.
;
; INPUTS:
;  chstr        The structure with information on each chip including
;                 the PHOTOMETRIC column and transformation equations.
;  ubercalstr   The structure with information on the ubercal-solved
;                  offsets for each chip.
;
; OUTPUTS:
;  The values in the ZPTERM column in CHSTR are updated.
;
; USAGE:
;  IDL>smashred_set_zeropoints,chstr,ubercalstr
;
; By D. Nidever  April 2016
;-

pro smashred_set_zeropoints,chstr,ubercalstr

; Not enough inputs
if n_elements(chstr) eq 0 or n_elements(ubercalstr) eq 0 then begin
  print,'Syntax - smashred_set_zeropoints,chstr,ubercalstr'
  return
endif

; Set absolute zeropoint using photometry and/or anchor data
gdphot = where(chstr.photometric eq 1,ngdphot)
if ngdphot eq 0 then begin
  print,'WARNING!!!  No photometric data to anchor the relative magnitude zeropoints'
  return
endif

; Calculate median mag offset for photometric data
;  Use median, don't remove outliers, basically want mean of all
;  photometric data
zpmagoff = median(ubercalstr[gdphot].magoff)
print,'Removing median mag offset for photometric data = ',strtrim(zpmagoff,2),' mag'
ubercalstr.magoff -= zpmagoff

; Apply the new zeropoints to the chip-level transformation equations
;  we want to add this offset to the magnitudes
;  zpterm is subtracted from instrumental mags, so subtract mag offset
chstr.zpterm -= ubercalstr.magoff

; DO WE UPDATED THE ZPTERM ERROR??

;stop

end
