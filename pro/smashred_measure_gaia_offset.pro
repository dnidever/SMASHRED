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
;
; OUTPUTS:
;  magoff      The magnitude offset.
;  magofferr   The error in the magnitude offset.
;
; USAGE:
;  IDL>smashred_measure_gaia_offset,gaia,allobj,'g',gaiastr,magoff,magofferr
;
; D. Nidever  Sep. 2016
;-

pro smashred_measure_gaia_offset,gaia,allobj,filter,gaiacolstr,magoff,magofferr

undefine,magoff
undefine,magofferr

; Not enough inputs
if n_elements(gaia) eq 0 or n_elements(allobj) eq 0 or n_elements(filter) eq 0 or $
   n_elements(gaiacolstr) eq 0 then begin
  print,'Syntax - smashred_measure_gaia_offset,gaia,allobj,filter,gaiacolstr,magoff,magofferr'
  return
endif

; Color terms
;ucoef = [0.79317753d0,  0.72506918d0,  1.4458541d0, -0.47691866d0]
;gcoef = [0.076057942d0, 0.64317519d0,  0.089028260d0,  -0.015999274d0]
;rcoef = [0.038497136d0,  0.093925705d0,  -0.12002761d0,  0.058146257d0]
;icoef = [0.077870198d0,  -0.36229412d0,   0.093845879d0,  -0.017227423d0]
;zcoef = [0.19274132d0,  -0.72134595d0,   0.21075378d0,  -0.050022044d0]

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

; Match up the stars
SRCMATCH,allobj.ra,allobj.dec,gaia.ra_ircs,gaia.de_ircs,0.5,ind1,ind2,/sph,count=nmatch
if nmatch eq 0 then begin
  print,'No matches'
  return
endif
mallobj = allobj[ind1]
mgaia = gaia[ind2]

; Get the good mags
gdmag = where(mallobj.(magind[0]) lt 50,ngdmag)
if ngdmag eq 0 then begin
  print,'No good magnitudes in ',filter
  return
endif
mallobj2 = mallobj[gdmag]
mgaia2 = mgaia[gdmag]

; Get the color
gdcol = where(mallobj2.g lt 50 and mallobj2.i lt 50,ngdcol)
if ngdcol gt 0 then begin
  col = mallobj2[gdcol].g-mallobj2[gdcol].i
  mag = mallobj2[gdcol].(magind[0])
  magerr = mallobj2[gdcol].(errind[0])
  gmag = mgaia2[gdcol]._gmag_
  gmagerr = 2.5*alog10(1.0+mgaia2[gdcol].e__fg_/mgaia2[gdcol]._fg_)
; no color info
endif else begin
  col = mallobj2.g*0
  mag = mallobj2.(magind[0])
  magerr = mallobj2.(errind[0])
  gmag = mgaia2._gmag_
  gmagerr = 2.5*alog10(1.0+mgaia2.e__fg_/mgaia2._fg_)
endelse

; Get the color terms
coef = gaiacolstr1.bestcoef

; x-G vs. g-i
; x = f(g-i)+G
xmodel = poly(col,coef)+gmag
diff = mag-xmodel
magoff = median(diff)
sig = mad(diff)
gdind = where( abs(diff-magoff) lt (sig>0.02),ngdind)
differr = sqrt(magerr[gdind]^2+gmagerr[gdind]^2)
magoff = dln_poly_fit(col[gdind],diff[gdind],0,measure_errors=differr,sigma=sigma,/bootstrap,status=status)
magofferr = sigma[0]


; NEED TO ITERATE!!!!


stop

end
