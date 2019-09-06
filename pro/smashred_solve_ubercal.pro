;+
;
; SMASHRED_SOLVE_UBERCAL
;
; Solve the ubercal problem with iteration
;
; INPUTS:
;  overlapstr  The structure containing the relative magnitude offsets
;                of chip pairs.
;  =chiptrans  A structure with the chip-level relative magnitude
;                offsets.  We use this for chips that don't have any
;                good overlap magnitude offsets.
;  /verbose    Print extra information to the screen.
;
; OUTPUTS:
;  ubercalstr  The ubercal structure with information on the final
;                magnitude offset, errors and number of overlaps.
;
; USAGE:
;  IDL>smashred_solve_ubercal,overlapstr,ubercalstr,verbose=verbose
;
; By D. Nidever  April 2016
;-

pro smashred_solve_ubercal,overlapstr,ubercalstr,chiptrans=chiptrans,verbose=verbose

undefine,ubercalstr

; Not enough inputs
if n_elements(overlapstr) eq 0 then begin
  print,'Syntax - smashred_solve_ubercal,overlapstr,ubercalstr,chiptrans=chiptrans,verbose=verbose'
  return
endif

sz = size(overlapstr.expnum)
nchips = sz[1]

; Initialize ubercal structure
ubercalstr = replicate({expnum:'',chip:-1L,magoff:0.0,magofferr:0.0,noverlap:-1L,flag:-1},nchips)
ubercalstr.expnum = overlapstr.expnum
ubercalstr.chip = overlapstr.chip

; No overlaps
if n_elements(overlapstr.expnum) eq 1 and overlapstr.noverlap[0] le 0 then begin
  print,'No overlaps'
  ubercalstr = {expnum:overlapstr.expnum,chip:long(overlapstr.chip),magoff:0.0,magofferr:0.0,noverlap:0L,flag:-1}
  return
endif

; Now iterate to put them all on the same system
data = overlapstr.data
count = 0
flag = 0
davg0 = 100
dmax0 = 100
fmagflag = lonarr(nchips)  ; 1-good, 0-no good
magoff = overlapstr.data.magoff  ; initialize
;bd = where(magoff gt 50,nbd,comp=gd)
;if nbd gt 0 then magoff[bd]=!values.f_nan  ; set bad values to NAN
WHILE (flag eq 0) do begin

  ; The method is basically to calculate for each chip how it
  ; (on average) compares to its overlapping neighbors.  Then
  ; remove that and iterate.  The TRICK is to only go halfway,
  ; otherwise you faciliate a lot.

  ; Compute weighted mean offset and remove
  ;  from the magoffset values
  mnoff = fltarr(nchips)
  for k=0,nchips-1 do begin
    ngd = 0
    if overlapstr.noverlap[k] gt 0 then begin
      ind = overlapstr.index[overlapstr.ind0[k]:overlapstr.ind1[k]]
      revind = overlapstr.revindex[overlapstr.ind0[k]:overlapstr.ind1[k]]
      ; Get good matches with low sigma
      gd = where(data[ind].nmatch gt 3 and data[ind].magofferr lt 0.05,ngd)
      if ngd eq 0 then $  ; raise the error threshold slightly
        gd = where(data[ind].nmatch gt 3 and data[ind].magofferr lt 0.07,ngd)
      ; Found some good mag offsets
      if ngd gt 0 then begin
        ; calculate the weighted mean
        ROBUST_MEAN,magoff[ind[gd]],robmean,robsig,sig=data[ind[gd]].magofferr
        if ngd eq 1 or robsig lt 1e-5 then robsig=median([data[ind[gd]].magofferr])>1e-2  ; lower limit to robsig
        if finite(robmean) eq 0 then robmean=median([magoff[ind[gd]]])  ; just in case
        ; remove offset from all overlaps 
        magoff[ind] -= robmean*0.5   ; offset all of the values
        magoff[revind] += robmean*0.5
        ubercalstr[k].magoff -= robmean*0.5
        ubercalstr[k].magofferr = robsig/sqrt(ngd)
        ubercalstr[k].noverlap = ngd
        mnoff[k] = robmean
        fmagflag[k] = 1  ; we have a good mag offset for this one
        if keyword_set(verbose) then print,count,k+1,ubercalstr[k].expnum,ubercalstr[k].chip,ngd,robmean,robsig,$
                                           format='(I5,I5,A11,I5,I7,F10.6,F10.6)'
      endif
    endif ; no overlap

    ; No good mag offsets
    if ngd eq 0 then begin
      if keyword_set(verbose) then print,count,k+1,ubercalstr[k].expnum,ubercalstr[k].chip,ngd,format='(I5,I5,A11,I5,I7)'
      ubercalstr[k].magoff = 99.99
      ubercalstr[k].magofferr = 9.99
    endif

  endfor ; chip loop

  ; How much have things changed
  davg = total(abs(mnoff))/nchips
  dmax = max(abs(mnoff))
  if keyword_set(verbose) then print,count,' avg. diff = ',strtrim(davg,2),'  max diff =',strtrim(dmax,2)

  ; Do we need to stop
  if count ge 300 or (abs(davg-davg0)/davg0 lt 1e-2 and abs(dmax-dmax0)/dmax0 lt 1e-2) then flag=1

  davg0 = davg
  dmax0 = dmax

  count++

  ;stop

ENDWHILE

; PRINT OUT FINAL AVERAGE DIFFERENCE
print,'Final avg. diff = ',strtrim(davg,2),' mag   max diff = ',strtrim(dmax,2),' mag'

; Set flag
ubercalstr.flag = fmagflag   ; 0-nothing set, 1-good weighted mean offset

; Calculate the relative offset for each chip
;  first calculate the relative mag offset from exposure median
rel_magoff = ubercalstr.magoff
uiexp = uniq(ubercalstr.expnum,sort(ubercalstr.expnum))
uexp = ubercalstr[uiexp].expnum
nuexp = n_elements(uexp)
expmagoff = fltarr(nuexp)
expmagofferr = fltarr(nuexp)
for k=0,nuexp-1 do begin
  ;; Measure median magoff for this exposure
  expindgd = where(ubercalstr.expnum eq uexp[k] and ubercalstr.magoff lt 50,nexpindgd)
  if nexpindgd gt 0 then begin
    expmagoff[k] = median([ubercalstr[expindgd].magoff])
    expmagofferr[k] = mad([ubercalstr[expindgd].magoff]) / sqrt(nexpindgd)  ; error in mean
  endif
  if keyword_set(verbose) then print,k+1,' ',uexp[k],' ',expmagoff[k],' ',stddev([ubercalstr[expindgd].magoff])
  ;; Measure relative offset for chips of this exposure
  expindall = where(ubercalstr.expnum eq uexp[k],nexpindall)
  if nexpindall gt 0 then rel_magoff[expindall] -= expmagoff[k]
endfor
; For unique chips measure the median relative offset
uichips = uniq(ubercalstr.chip,sort(ubercalstr.chip))
uchips = ubercalstr[uichips].chip
nuchips = n_elements(uchips)
chipdeltamagoff = fltarr(nuchips)
chipdeltamagofferr = fltarr(nuchips)+9.99
for k=0,nuchips-1 do begin
  chind = where(ubercalstr.chip eq uchips[k] and abs(rel_magoff) lt 10,nchind)
  if nchind gt 0 then begin
    chipdeltamagoff[k] = median([rel_magoff[chind]])
    chipdeltamagofferr[k] = mad([rel_magoff[chind]]) / sqrt(nchind)  ; error in mean
  endif
endfor

; Fix chips with no overlap or no good offsets
;   use the median offset for all chips of that exposure
;totoverlap = total( (overlapstr.overlap eq 1 and overlapstr.magoff lt 50),1)  ; 1 for chips with good offsets, 0 otherwise
;bdind = where(totoverlap eq 0 or (totoverlap gt 0 and ubercalstr.magoff gt 50),nbdind)
; The relative chip offsets derived above tend to be better
; than the relative ones derived for the whole survey.
bdind = where(ubercalstr.flag eq 0,nbdind)
for k=0,nbdind-1 do begin
  expind = where(uexp eq ubercalstr[bdind[k]].expnum,nexpind)
  chind = where(uchips eq ubercalstr[bdind[k]].chip,nchind)
  ; Some values for this chip and exosure
  if nexpind gt 1 and nchind gt 1 then begin
    ind = overlapstr.index[overlapstr.ind0[bdind[k]]:overlapstr.ind1[bdind[k]]]   
    ubercalstr[bdind[k]].magoff = expmagoff[expind] + chipdeltamagoff[chind]
    ubercalstr[bdind[k]].magofferr = sqrt( expmagofferr[expind]^2 + chipdeltamagofferr[chind]^2 )  ; add errors in quadrature
    ubercalstr[bdind[k]].noverlap = total(overlapstr.data[ind].overlap eq 1)
    ubercalstr[bdind[k]].flag = 2  ; set the flag
  ; Use value from CHIPTRANS
  ;  Use chiptrans only if necessary
  endif else begin
    MATCH,chiptrans.chip,ubercalstr[bdind[k]].chip,ind1,ind2,/sort,count=nmatch
    if nexpind gt 1 then begin  ; there are other chips in this exp
      ubercalstr[bdind[k]].magoff = expmagoff[expind] + chiptrans[ind1].zpterm
      ubercalstr[bdind[k]].magofferr = sqrt( expmagofferr[expind]^2 + chiptrans[ind1].zptermerr^2 )  ; add errors in quadrature
    endif else begin
      ubercalstr[bdind[k]].magoff = chiptrans[ind1].zpterm
      ubercalstr[bdind[k]].magofferr = chiptrans[ind1].zptermerr
    endelse
    nind = overlapstr.ind1[bdind[k]]-overlapstr.ind0[bdind[k]]+1
    if nind gt 0 then begin  ; do we need this?
      ind = overlapstr.index[overlapstr.ind0[bdind[k]]:overlapstr.ind1[bdind[k]]]   
      ubercalstr[bdind[k]].noverlap = total(overlapstr.data[ind].overlap eq 1)
    endif
    ubercalstr[bdind[k]].flag = 3  ; set the flag
  endelse
endfor

;; Check for disconnected "continents" of chips
;;----------------------------------------------
;;  they might be split into contiguously connected "continents" with
;;  gaps between them
gover = where(overlapstr.noverlap gt 0,ngover)
If ngover gt 0 then begin
  group = intarr(n_elements(overlapstr.chip))
  smashred_find_connected_chips,overlapstr,index,count=nindex
  if nindex gt 0 then group[index] = 1

  ;; Loop until done
  doneflag = 0
  niter = 0
  last_group = group
  WHILE (doneflag eq 0) do begin
    ;; -1 are known single chips, 0 not sure yet
    left = where(group eq 0 and overlapstr.noverlap gt 0,nleft)
    if nleft gt 0 then begin
      smashred_find_connected_chips,overlapstr,index,start=left[0],count=nindex
      if nindex gt 1 then group[index] = max(group)+1 > 1
      if nindex eq 1 then group[index] = -1   ;; single chip
    endif
    if nleft eq 0 then doneflag=1
    diff_group = last_group - group
    if total(abs(diff_group)) eq 0 then doneflag=1
    last_group = group
    niter++
  ENDWHILE

  ;; Multiple groups/continents
  ;;----------------------------
  ngroups = max(group)>1
  if ngroups gt 1 then begin
    print,'  There are ',strtrim(ngroups,2),' groups of connected chips'
    ;; build group index that excludes chips with no overlap
    ingroups = where(group gt 0,numingroups)
    grpindex = {index:lonarr(numingroups),value:strarr(ngroups),$
                num:lonarr(ngroups),lo:lonarr(ngroups),hi:lonarr(ngroups)}
    grpindex.value = lindgen(ngroups)+1
    count = 0L
    for i=0,ngroups-1 do begin
      ind = where(group eq i+1,nind)
      grpindex.num[i] = nind
      grpindex.lo[i] = count
      grpindex.hi[i] = count+nind-1
      grpindex.index[count:count+nind-1] = ind
      count += nind
      print,'   GROUP'+strtrim(i+1,2),': ',strtrim(nind,2),' chips'
    endfor

    ;; We can connect them if the chips belong to the same exposure
    ;; Loop over exposures
    eindex = create_index(overlapstr.expnum)
    expgrpmagoff = fltarr(nuexp,ngroups)+!values.f_nan
    for i=0,nuexp-1 do begin
      ind = eindex.index[eindex.lo[i]:eindex.hi[i]]
      egrp = group[ind]
      bindata,egrp-1,ubercalstr[ind].magoff,xbin,ybin,bin=1,/med,xmnbin=xmnbin
      if n_elements(xbin) gt 1 then $
        expgrpmagoff[i,long(xmnbin)] = ybin
    endfor
  
    ;; Apply offsets until we have converged
    if max(total(finite(expgrpmagoff),2)) gt 0 then begin
      doneflag = 0
      niter = 0
      last_expgrpmagoff = expgrpmagoff
      WHILE (doneflag eq 0) do begin
        ;; Subtract median across each exposure
        medexp = median(expgrpmagoff,dim=2,/even)
        for i=0,nuexp-1 do begin
           if finite(medexp[i]) then expgrpmagoff[i,*]-=medexp[i]
        endfor
        ;; Calculate mean offset for this group across all exposures and apply
        medgrp = median(expgrpmagoff,dim=1,/even)
        ;; Subtract from each groups chips MAGOFF and from EXPGRPMAGOFF as well
        for i=0,ngroups-1 do begin
          if finite(medgrp[i]) eq 1 then begin
            ind = grpindex.index[grpindex.lo[i]:grpindex.hi[i]]
            ubercalstr[ind].magoff -= medgrp[i]
            expgrpmagoff[*,i] -= medgrp[i]
          endif
        endfor
        ;; Are we done?
        diff_expgrpmagoff = last_expgrpmagoff - expgrpmagoff
        maxdiff = max(abs(diff_expgrpmagoff),/nan)
        if maxdiff lt 1e-4 or finite(maxdiff) eq 0 then doneflag=1
        last_expgrpmagoff = expgrpmagoff
        ;print,niter+1,' ',max(abs(diff_expgrpmagoff),/nan)
        niter++
      ENDWHILE
    endif  ;; exposure-level offsets to apply
  endif ;; multiple groups/continents
Endif  ;; some overlaps, check for "continents"

; Some chips have no magnitude offsets
bdflag = where(ubercalstr.flag eq 0,nbdflag)
if nbdflag gt 0 then print,strtrim(nbdflag,2),'/',strtrim(nchips,2),' chips have NO measured offsets'

;stop

end
