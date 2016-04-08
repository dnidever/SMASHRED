;+
;
; SMASHRED_SOLVE_UBERCAL
;
; Solve the ubercal problem with iteration
;
; INPUTS:
;  overlapstr  The structure containing the relative magnitude offsets
;                of chip pairs.
;  /verbose    Print extra information to the screen.
;
; OUTPUTS:
;  fmagoff     The array of relative photometric offsets, one per chip.
;
; USAGE:
;  IDL>smashred_solve_ubercal,overlapstr,fmagoff,verbose=verbose
;
; By D. Nidever  April 2016
;-

pro smashred_solve_ubercal,overlapstr,fmagoff,verbose=verbose

undefine,fmagoff

; Not enough inputs
if n_elements(overlapstr) eq 0 then begin
  print,'Syntax - smashred_solve_ubercal,overlapstr,fmagoff,verbose=verbose'
  return
endif

sz = size(overlapstr)
nchips = sz[1]

; Now iterate to put them all on the same system
count = 0
flag = 0
davg0 = 100
dmax0 = 100
fmagoff = fltarr(nchips)
fmagflag = lonarr(nchips)  ; 1-good, 0-no good
magoff = overlapstr.magoff  ; initialize
bd = where(magoff gt 50,nbd,comp=gd)
if nbd gt 0 then magoff[bd]=!values.f_nan
WHILE (flag eq 0) do begin

  ; maybe calculate for each chip how (on average) it compares to
  ; its neighbors, then remove that and start again
  medoff = median(magoff,dim=1,/even)
  ; we actually use the weighted mean now, but just in case

  ; update the catalogs with these offsets
  ; the trick is only to go halfway, otherwise we facilate a lot
  ;for k=0,nchips-1 do (*chfiltstr[k].data).mag += medoff[k]*0.5

  ; Compute weighted mean offset and remove
  ;  from the magoffset values
  mnoff = fltarr(nchips)
  for k=0,nchips-1 do begin
    ; get good matches with low sigma
    gd = where(finite(reform(magoff[*,k])) eq 1 and overlapstr[*,k].nmatch gt 3 and overlapstr[*,k].magofferr lt 0.05,ngd)
    ;gd = where(finite(reform(magoff[*,k])) eq 1 and overlapstr[*,k].nmatch gt 10 and overlapstr[*,k].magofferr lt 0.05,ngd)
;if ngd eq 0 then stop,'no good matches'
;if max(magoff[gd,k]) gt 0.5 then stop,'big offset'
    if ngd gt 0 then begin
      ; calculate the weighted mean
      ROBUST_MEAN,magoff[gd,k],robmean,robsig,sig=overlapstr[gd,k].magofferr
      if finite(robmean) eq 0 then robmean=medoff[k] ; just in case
      magoff[k,gd] += robmean*0.5
      magoff[gd,k] -= robmean*0.5
      fmagoff[k] += robmean*0.5
      mnoff[k] = robmean
      fmagflag[k] = 1
    endif
  endfor

  ; How much have things changed
  davg = total(abs(mnoff))/nchips
  dmax = max(abs(mnoff))
  if keyword_set(verbose) then print,count,davg,dmax

  ; Do we need to stop
  if count ge 300 or (abs(davg-davg0)/davg0 lt 1e-2 and abs(dmax-dmax0)/dmax0 lt 1e-2) then flag=1

  davg0 = davg
  dmax0 = dmax

  count++

  ;stop

ENDWHILE

; PRINT OUT FINAL AVERAGE DIFFERENCE
print,'Final avg. diff = ',strtrim(davg,2),' mag   max diff = ',strtrim(dmag,2),' mag'

; Calculate the relative offset for each chip
;  first calculate the relative mag offset from exposure median
rel_fmagoff = fmagoff
;uiexp = uniq(chfiltstr.expnum,sort(chfiltstr.expnum))
allexpnum = overlapstr[*,0].expnum1
uiexp = uniq(allexpnum,sort(allexpnum))
uexp = allexpnum[uiexp]
nuexp = n_elements(uexp)
for k=0,nuexp-1 do begin
  expind = where(allexpnum eq uexp[k],nexpind)
  if nexpind gt 0 then rel_fmagoff[expind] -= median([fmagoff[expind]])
endfor
allchips = overlapstr[*,0].chip1
uichips = uniq(allchips,sort(allchips))
uchips = allchips[uichips]
nuchips = n_elements(uchips)
deltamagoff_chip = fltarr(nuchips)
for k=0,nuchips-1 do begin
  chind = where(allchips eq uchips[k],nchind)
  if nchind gt 0 then deltamagoff_chip[k] = median([rel_fmagoff[chind]])
endfor

; chips with no overlap
;   use the median offset for all chips of that exposure
totoverlap = total( (overlapstr.overlap eq 1 and overlapstr.magoff lt 50),1)
bdoverlap = where(totoverlap eq 0,nbdoverlap)
for k=0,nbdoverlap-1 do begin
  ;expind = where(chfiltstr.expnum eq chfiltstr[bdoverlap[k]].expnum and totoverlap gt 0 and fmagflag eq 1,nexpind)
  ;chind = where(uchips eq chfiltstr[bdoverlap[k]].chip,nchind)
  expind = where(allexpnum eq allexpnum[bdoverlap[k]] and totoverlap gt 0 and fmagflag eq 1,nexpind)
  chind = where(uchips eq allchips[bdoverlap[k]],nchind)
;stop
  fmagoff[bdoverlap[k]] = median([fmagoff[expind]]) + deltamagoff_chip[chind]
endfor


;stop

end
