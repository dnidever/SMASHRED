pro run_model_lmc_field

; Run model_lmc_field on multiple LMC fields

  ;; fields where we see the LMC main-sequence
  smash_detect = 'Field'+strtrim([56, 61, 60, 64, 33, 31, 25, 26, 24, 69, 66, 67, 63, 59, 157, 58, 57, 246, 156],2)
  ;match,smash.field,smash_detect,ind1,ind2,/sort
  ;smash[ind1].lmcdetect = 1
  ;; not clear
  smash_notclear = 'Field'+strtrim([65, 23, 17, 62, 68, 155, 154, 153, 159, 87, 115, 106, 247],2)
  ;match,smash.field,smash_notclear,ind1,ind2,/sort
  ;smash[ind1].lmcdetect = 2
  ;; no 
  ;smash_nondetect = 'Field'+strtrim([83, 80, 130],1)
  ;match,smash.field,smash_nondetect,ind1,ind2,/sort
  ;smash[ind1].lmcdetect = 0

  fields = [smash_detect,smash_notclear]
  nfields = n_elements(fields)
  for i=0,nfields-1 do begin
    model_lmc_field,fields[i]
  endfor  

stop

end
