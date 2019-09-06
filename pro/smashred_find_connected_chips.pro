pro smashred_find_connected_chips,overlapstr,cindex,start=start,count=count

;; Check which chips are connected/overlap in one contiguous region
  
cindex = -1
count = 0
  
if n_elements(start) gt 0 then $
   if overlapstr.noverlap[start] eq 0 then return

;; Pick START value
if n_elements(start) eq 0 then begin
  gover = where(overlapstr.noverlap gt 0,ngover)
  if ngover eq 0 then return
  start = gover[0]
endif

expnum = overlapstr.expnum
chip = overlapstr.chip
data = overlapstr.data
  
;; Start CONNECT array
connect = intarr(n_elements(expnum))
connect[start] = 1

;; Loop until done
doneflag = 0
niter = 0
last_connect = connect
WHILE (doneflag eq 0) do begin

  ;; Get connected chips
  MATCH,connect,1,ind1,ind2,/sort,count=nmatch
  ;; Loop through the calibrated chips and find overlap chips
  For i=0,nmatch-1 do begin
    ;; Find overlap chips with good crossmatches
    ind = overlapstr.index[overlapstr.ind0[ind1[i]]:overlapstr.ind1[ind1[i]]]
    gd = where(data[ind].nmatch gt 3 and data[ind].magofferr lt 0.05,ngd)
    ;; Some good overlaps
    if ngd gt 0 then begin
       MATCH,expnum+'-'+strtrim(chip,2),$
             data[ind[gd]].expnum2+'-'+strtrim(data[ind[gd]].chip2,2),indx1,indx2,/sort
       connect[indx1] = 1
    endif
  Endfor

  ;; Are we done?
  maxiter = 1000
  diff_connect = last_connect - connect
  if total(abs(diff_connect)) eq 0 or niter gt maxiter then doneflag=1
  last_connect = connect

  ;; Increment the counter
  niter++
ENDWHILE

;; This index INCLUDES start
cindex = where(connect eq 1,count)

;stop

end
  
