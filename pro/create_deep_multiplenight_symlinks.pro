pro create_deep_multiplenight_symlinks

; Create symlinks in deep/ for Fields
; that were observed on multiple nights
; and were reprocessed.  This is to save
; space.

rootdir = smashred_rootdir()+'cp/red/photred/'

check = mrdfits(rootdir+'catalogs/pro/check_calibrated_v5.fits',1)
check.field = strtrim(check.field,2)
ncheck = n_elements(check)
smashred_getredinfo,info,/silent
for i=0,ncheck-1 do begin
  gd = where(info.field eq check[i].field and info.sh eq 0,ngd)
  if ngd gt 1 then push,fields,check[i].field
endfor

;; Remove the short LMC fields
;fnum = long(strmid(check.field,5))
;;shlmc = where(check.nuchips eq 0 and check.ngchips lt 150 and check.nrchips lt 150 and check.nichips lt 150 and check.nzchips lt 150,nshlmc)
;gshlmc = where(fnum ge 184 and fnum lt 246 and check.calflag gt 0,ngshlmc)
;shlmcfields = check[gshlmc].field
;MATCH,shlmcfields,fields,ind1,ind2,/sort,count=nmatch
;if nmatch gt 0 then remove,ind2,fields

;fieldnum = [100,101,104,109,113,116,76,77,84,92,98]
;fields = 'Field'+strtrim(fieldnum,2)
nfields = n_elements(fields)

;smashred_getredinfo,info,/silent

for i=0,nfields-1 do begin
  ifield = fields[i]
  print,strtrim(i+1,2),' ',ifield
  ind = where(info.field eq ifield and info.sh eq 0,nind)
  if nind eq 1 then begin
    stop,'only one night'
  endif
  object = info[ind[0]].object   ; original object name
  info1 = info[ind[0]]
  indir = file_dirname(info1.file)+'/'
  outdir = rootdir+'deep/'+ifield+'/'
  origfield = info1.object

  ;deep/Field110/F1/F1-00519030_01
  ;20160218/F2/F2-00519030_01
  ;The header of the FITS file is slightly modified because the AIRMASS line was updated
  ;(just the date in the comment changed), which changed the CHECKSUM line and the
  ;DATASUM comment line.
  ;All other files are identical except for the '.alf' and '.log' files.

  stop

  ; Link to summary.fits and final catalog
  file_delete,outdir+ifield+'_summary.fits'
  file_link,info1.file,outdir+ifield+'_summary.fits'
  file_delete,outdir+ifield+'.fits.gz'
  file_link,indir+origfield+'.fits.gz',outdir+ifield+'.fits.gz'

  ; Link to apcor.lst and photred.setup file
  if file_test(outdir+'apcor.lst') eq 0 then $
    file_link,indir+'apcor.lst',outdir+'apcor.lst'
  if file_test(outdir+'photred.setup') eq 0 then $
    file_link,indir+'photred.setup',outdir+'photred.setup'

  ; Create new "fields" file
  readcol,indir+'fields',shnames,lnames,format='A,A',/silent
  ;gdname = where(lnames eq ifield,ngdname)
  ;if ngdname eq 0 then stop,ifield,' not found'
  gdname = where(lnames eq object,gdname)
  if ngdname eq 0 then stop,object,' not found'
  shnames1 = shnames[gdname[0]]

  ; Link to subdirectory with the data
  if file_test(outdir+shnames1+'/') eq 0 then $
    file_link,indir+shnames1+'/',outdir+shnames1+'/'

  ;stop

  BOMB:

endfor

stop

end
