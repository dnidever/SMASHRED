pro create_deep_singlenight_symlinks

; Create symlinks in deep/ for Fields
; that were observed on only ONE night
; and didn't need PHOTRED reprocessing

;       0  Field100  20150317 20150317
;       1  Field101  20150318 20150318
;       2  Field104  20140602 20140602
;       3  Field109  20140602 20140602
;       5  Field113  20160218 20160218
;       6  Field116  20160216 20160216
;      33  Field76  20150316 20150316
;      34  Field77  20160216 20160216
;      36  Field84  20160217 20160217
;      40  Field92  20150317 20150317
;      43  Field98  20140527 20140527

rootdir = smashred_rootdir()+'cp/red/photred/'

check = mrdfits(rootdir+'catalogs/pro/check_calibrated_v4.fits',1)
check.field = strtrim(check.field,2)
ncheck = n_elements(check)
smashred_getredinfo,info,/silent
for i=0,ncheck-1 do begin
  gd = where(info.field eq check[i].field and info.sh eq 0,ngd)
  if ngd eq 1 then push,fields,check[i].field
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
  if nind gt 1 then begin
    stop,'more than one night'
  endif
  info1 = info[ind[0]]
  indir = file_dirname(info1.file)+'/'
  outdir = rootdir+'deep/'+ifield+'/'
  origfield = info1.object

  ;if ifield eq 'Field157' then stop,ifield
  ;if ifield eq 'Field157' then goto,BOMB
  ;if ifield eq 'Field176' then goto,BOMB

  ;if file_test(outdir,/directory) eq 1 then begin
  ;  print,outdir,' already exists.  Skipping'
  ;  goto,BOMB
  ;endif

  ; Create directory
  if file_test(outdir,/directory) eq 0 then file_mkdir,outdir

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
  gdname = where(lnames eq ifield,ngdname)
  shnames1 = shnames[gdname[0]]

  ; Link to subdirectory with the data
  if file_test(outdir+shnames1+'/') eq 0 then $
    file_link,indir+shnames1+'/',outdir+shnames1+'/'

  ;stop

  BOMB:

endfor

stop

end
