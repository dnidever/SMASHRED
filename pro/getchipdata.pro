function getchipdata,chstr,allsrc,expnum,chip
  if n_elements(chip) eq 0 then begin
    dum = strtsplit(expnum,'-',/extract)
    expnum = dum[0]
    chip = dum[1]
  endif
  gd = where(chstr.expnum eq expnum and chstr.chip eq chip,ngd)
  if ngd eq 0 then begin
    print,expnum,' ',chip,' NOT FOUND'
    return,-1
  endif
  return,allsrc[chstr[gd[0]].allsrcindx:chstr[gd[0]].allsrcindx+chstr[gd[0]].nsrc-1]

end
