;+
;
; SMASHRED_MEASURE_GAIA_OFFSET
;
; Measure a photometric magnitude for a single band
; for the SMASH survey using the GAIA G-band magnitudes
; and SMASH-GAIA color terms.
;
; INPUTS:
;  gaia        The structure of GAIA sources.
;  allobj      The structure of average SMASH source values.
;  filter      The filter name.
;  gaiacolstr  The structure of GAIA-SMASH color terms.
;  /silent     Don't print anything to the screen.
;
; OUTPUTS:
;  magoff      The magnitude offset.  This value should be
;                 SUBTRACTED from the data to correct it.
;  magofferr   The error in the magnitude offset.
;
; USAGE:
;  IDL>smashred_measure_gaia_offset,gaia,allobj,'g',gaiastr,magoff,magofferr
;
; D. Nidever  Sep. 2016
;-

pro smashred_measure_gaia_offset,gaia,allobj,filter,gaiacolstr,magoff,magofferr,silent=silent

undefine,magoff
undefine,magofferr

; Not enough inputs
if n_elements(gaia) eq 0 or n_elements(allobj) eq 0 or n_elements(filter) eq 0 or $
   n_elements(gaiacolstr) eq 0 then begin
  print,'Syntax - smashred_measure_gaia_offset,gaia,allobj,filter,gaiacolstr,magoff,magofferr,silent=silent'
  return
endif

; Get GAIA-SMASH color terms for this filter
colind = where(gaiacolstr.filter eq filter,ncolind)
if ncolind eq 0 then begin
  print,'No GAIA-SMASH color terms found for Filter=',filter
  return
endif
gaiacolstr1 = gaiacolstr[colind[0]]

tags = tag_names(allobj)
magind = where(tags eq strupcase(filter),nmagind)
errind = where(tags eq strupcase(filter+'err'),nerrind)
colmagind1 = where(tags eq strupcase(gaiacolstr1.colnames[0]),ncolmagind1)
colmagind2 = where(tags eq strupcase(gaiacolstr1.colnames[1]),ncolmagind2)

; Match up the stars
SRCMATCH,allobj.ra,allobj.dec,gaia.ra_icrs,gaia.de_icrs,0.5,ind1,ind2,/sph,count=nmatch
if nmatch eq 0 then begin
  print,'No matches'
  return
endif
mallobj = allobj[ind1]
mgaia = gaia[ind2]
if not keyword_set(silent) then print,'  ',strtrim(nmatch,2),' matches to GAIA'

; Get the good mags and sharp/chi
magmask = (mallobj.(magind[0]) lt 50 and abs(mallobj.sharp) lt 1 and mallobj.chi lt 4)
ngdmag = long(total(magmask))
if ngdmag eq 0 then begin
  print,'No good magnitudes in ',filter
  return
endif

; Get the color
colmask = (mallobj.(colmagind1) lt 50 and mallobj.(colmagind2) lt 50)
ngdcol = long(total(colmask))
if ngdcol gt 0 then begin
  gdphot = where(magmask eq 1 and colmask eq 1 and mallobj.(colmagind1)-mallobj.(colmagind2) ge gaiacolstr1.bestcolrange[0] and $
                 mallobj.(colmagind1)-mallobj.(colmagind2) le gaiacolstr1.bestcolrange[1],ngdphot)
  col = mallobj[gdphot].(colmagind1)-mallobj[gdphot].(colmagind2)
  mag = mallobj[gdphot].(magind[0])
  magerr = mallobj[gdphot].(errind[0])
  gmag = mgaia[gdphot]._gmag_
  gmagerr = 2.5*alog10(1.0+mgaia[gdphot].e__fg_/mgaia[gdphot]._fg_)
; no color info
endif else begin
  gdphot = where(magmask eq 1,ngdphot)
  col = mallobj[gdphot].g*0
  mag = mallobj[gdphot].(magind[0])
  magerr = mallobj[gdphot].(errind[0])
  gmag = mgaia[gdphot]._gmag_
  gmagerr = 2.5*alog10(1.0+mgaia[gdphot].e__fg_/mgaia[gdphot]._fg_)
endelse

; Get the color terms
coef = gaiacolstr1.bestcoef

; x-G vs. g-i
; x = f(g-i)+G
magmodel = poly(col,coef)+gmag
diff = mag-magmodel
magoff0 = median(diff)
sig = mad(diff)
gdind = where( abs(diff-magoff0) lt (sig>0.02),ngdind)
differr = sqrt(magerr[gdind]^2+gmagerr[gdind]^2)
magoff = dln_poly_fit(col[gdind],diff[gdind],0,measure_errors=differr,sigma=sigma,/bootstrap,status=status)
magoff = magoff[0]
magofferr = sigma[0]

if not keyword_set(silent) then $
  print,'  GAIA offset = ',strtrim(string(magoff,format='(F8.4)'),2),' +/- ',strtrim(string(magofferr,format='(F8.5)'),2),' mag'

;stop

end
