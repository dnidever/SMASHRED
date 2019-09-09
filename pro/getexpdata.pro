function getexpdata,chstr,allsrc,expnum
  gd = where(chstr.expnum eq expnum,ngd)
  if ngd eq 0 then begin
    print,expnum,' ',' NOT FOUND'
    return,-1
  endif
  chstr1 = chstr[gd]
  ntot = total(chstr1.nsrc,/integer)
  ind = lonarr(ntot)
  cnt = 0LL
  for i=0,ngd-1 do begin
    n = chstr1[i].nsrc
    lo = chstr1[i].allsrcindx
    ind[cnt:cnt+n-1] = lindgen(n)+lo
    cnt += n
  endfor
  return,allsrc[ind]

end
