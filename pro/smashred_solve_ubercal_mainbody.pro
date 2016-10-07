;+
;
; SMASHRED_SOLVE_UBERCAL_MAINBODY
;
; Solve the ubercal problem with iteration for the main body field
; overlap structure.
;
; INPUTS:
;  overlapstr  The structure containing the relative magnitude offsets
;                of chip pairs.
;  /verbose    Print extra information to the screen.
;
; OUTPUTS:
;  ubercalstr  The ubercal structure with information on the final
;                magnitude offset, errors and number of overlaps.
;
; USAGE:
;  IDL>smashred_solve_ubercal_mainbody,overlapstr,ubercalstr,verbose=verbose
;
; By D. Nidever  April 2016
;-

pro smashred_solve_ubercal_mainbody,overlapstr,ubercalstr,verbose=verbose

undefine,ubercalstr

; Not enough inputs
if n_elements(overlapstr) eq 0 then begin
  print,'Syntax - smashred_solve_ubercal_mainbody,overlapstr,ubercalstr,verbose=verbose'
  return
endif

sz = size(overlapstr.field)
nfields = sz[1]

; Initialize ubercal structure
ubercalstr = replicate({field:'',magoff:fltarr(5),magofferr:fltarr(5),noverlap:lonarr(5),flag:intarr(5)-1},nfields)
ubercalstr.field = overlapstr.field

filters = ['u','g','r','i','z']
nfilters = n_elements(filters)

; Filter loop
For f=0,nfilters-1 do begin
  ifilter = filters[f]
  print,'FILTER = ',ifilter

  ; Now iterate to put them all on the same system
  data = overlapstr.data
  count = 0
  flag = 0
  davg0 = 100
  dmax0 = 100
  fmagflag = lonarr(nfields)  ; 1-good, 0-no good
  magoff = overlapstr.data.magoff[f]  ; initialize
  ;bd = where(magoff gt 50,nbd,comp=gd)
  ;if nbd gt 0 then magoff[bd]=!values.f_nan  ; set bad values to NAN
  WHILE (flag eq 0) do begin

    ; The method is basically to calculate for each chip how it
    ; (on average) compares to its overlapping neighbors.  Then
    ; remove that and iterate.  The TRICK is to only go halfway,
    ; otherwise you faciliate a lot.

    ; Compute weighted mean offset and remove
    ;  from the magoffset values
    mnoff = fltarr(nfields)
    for k=0,nfields-1 do begin
      ngd = 0
      if overlapstr.noverlap[k] gt 0 then begin
        ind = overlapstr.index[overlapstr.ind0[k]:overlapstr.ind1[k]]
        revind = overlapstr.revindex[overlapstr.ind0[k]:overlapstr.ind1[k]]
        ; Get good matches with low sigma
        gd = where(data[ind].ngood[f] gt 3 and data[ind].magofferr[f] lt 0.05,ngd)
        if ngd eq 0 then $  ; raise the error threshold slightly
          gd = where(data[ind].ngood[f] gt 3 and data[ind].magofferr[f] lt 0.07,ngd)
        ; Found some good mag offsets
        if ngd gt 0 then begin
          ; calculate the weighted mean
          ROBUST_MEAN,magoff[ind[gd]],robmean,robsig,sig=data[ind[gd]].magofferr[f]
          if ngd eq 1 or robsig lt 1e-5 then robsig=median([data[ind[gd]].magofferr[f]])>1e-2  ; lower limit to robsig
          if finite(robmean) eq 0 then robmean=median([magoff[ind[gd]]])  ; just in case
          ; remove offset from all overlaps 
          magoff[ind] -= robmean*0.5   ; offset all of the values
          magoff[revind] += robmean*0.5
          ubercalstr[k].magoff[f] -= robmean*0.5
          ubercalstr[k].magofferr[f] = robsig/sqrt(ngd)
          ubercalstr[k].noverlap[f] = ngd
          mnoff[k] = robmean
          fmagflag[k] = 1  ; we have a good mag offset for this one
          if keyword_set(verbose) then print,count,k+1,ubercalstr[k].field,ngd,robmean,robsig,$
                                             format='(I5,I5,A11,I7,F10.6,F10.6)'
        endif
      endif ; no overlap

      ; No good mag offsets
      if ngd eq 0 then begin
        if keyword_set(verbose) then print,count,k+1,ubercalstr[k].field,ngd,format='(I5,I5,A11,I7)'
        ubercalstr[k].magoff[f] = 99.99
        ubercalstr[k].magofferr[f] = 9.99
      endif

    endfor ; chip loop

    ; How much have things changed
    davg = total(abs(mnoff))/nfields
    dmax = max(abs(mnoff))
    if keyword_set(verbose) then print,count,' avg. diff = ',strtrim(davg,2),'  max diff =',strtrim(dmax,2)

    ; Do we need to stop
    if count ge 300 or (abs(davg-davg0)/davg0 lt 1e-2 and abs(dmax-dmax0)/dmax0 lt 1e-2) then flag=1

    davg0 = davg
    dmax0 = dmax

    count++

    ;stop

  ENDWHILE

  ; Print out the final offsets
  if not keyword_set(silent) then begin
    print,'' & print,'Final offsets for filter=',ifilter
    for k=0,nfields-1 do begin
      if overlapstr.noverlap[k] gt 0 then begin
        print,k+1,ubercalstr[k].field,ubercalstr[k].noverlap[f],ubercalstr[k].magoff[f],ubercalstr[k].magofferr[f],$
                                             format='(I5,A11,I7,F10.6,F10.6)'
      endif else begin
        print,k+1,ubercalstr[k].field,ubercalstr[k].noverlap[f],format='(I5,A11,I7)'
      endelse
    endfor
  endif

  ; PRINT OUT FINAL AVERAGE DIFFERENCE
  print,' Final avg. diff = ',strtrim(davg,2),' mag   max diff = ',strtrim(dmax,2),' mag'

  ; Set flag
  ubercalstr.flag[f] = fmagflag   ; 0-nothing set, 1-good weighted mean offset

  ; Some chips have no magnitude offsets
  bdflag = where(ubercalstr.flag[f] eq 0,nbdflag)
  if nbdflag gt 0 then print,strtrim(nbdflag,2),'/',strtrim(nfields,2),' fields have NO measured offsets'

Endfor ; filter loop

;stop

end
