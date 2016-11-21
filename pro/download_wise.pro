pro download_wise,field,redo=redo,compress=compress

outdir = '/data/smash/cp/red/photred/wise/'
outfile = outdir+field+'_wise.fits'
if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  print,outfile,' EXISTS and /redo not set'
  return
endif

smashred_getredinfo,info,/silent
ind = where(info.field eq field,nind)
info1 = info[ind]

; Load the summary chip files
undefine,chips
for i=0,n_elements(info1)-1 do begin
  chips1 = mrdfits(info1[i].file,2,/silent)
  push,chips,chips1
endfor

cendec = mean(minmax(chips.dec))
decr = range(chips.dec)*1.1 > 2.3
if range(chips.ra) gt 180 then begin
  ra = chips.ra
  over = where(ra gt 180,nover,comp=under,ncomp=nunder)
  if nover gt 0 then ra[over]-=360
  cenra = mean(minmax(ra))
  rar = range(ra)*cos(cendec/!radeg)*1.1 > 2.3
endif else begin
  cenra = mean(minmax(chips.ra))
  rar = range(chips.ra)*cos(cendec/!radeg)*1.1 > 2.3
endelse
print,field,' ',cenra,cendec,rar,decr
wise = queryvizier('ALLWISE',[cenra,cendec],[rar*60,decr*60],/canada,/all)
print,strtrim(n_elements(wise),2),' ALLWISE sources found'

MWRFITS,wise,outfile,/create
if keyword_set(compress) then spawn,['gzip',outfile],/noshell

end
