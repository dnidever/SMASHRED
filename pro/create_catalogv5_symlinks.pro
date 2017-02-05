pro create_catalogv5_symlinks

; Create symlinks in catalogs/final/v5/ to v4/ for those fields
; that didn't need to be reprocessed.

rootdir = smashred_rootdir()+'cp/red/photred/'

check = mrdfits(rootdir+'catalogs/pro/check_calibrated_v4.fits',1)
check.field = strtrim(check.field,2)
ncheck = n_elements(check)
smashred_getredinfo,info,/silent

; Get the 73 fields with deep data observed on duplicate nights
for i=0,ncheck-1 do begin
  gd = where(info.field eq check[i].field and info.sh eq 0,ngd)
  if ngd gt 1 then push,dupfields,check[i].field
endfor
fields = check.field
MATCH,fields,dupfields,ind1,ind2,/sort,count=nmatch
remove,ind1,fields
nfields = n_elements(fields)


for i=0,nfields-1 do begin
  ifield = fields[i]
  print,strtrim(i+1,2),' ',ifield

  v4dir = rootdir+'catalogs/final/v4/'
  v5dir = rootdir+'catalogs/final/v5/'
  srcfiles = v4dir+ifield+'_combined_'+['exposures.fits.gz','chips.fits.gz','allsrc.fits.gz','allobj.fits.gz',$
                                        'expmap.fits.gz','allobj_xmatch.fits.gz','allobj_bright.fits']
  relsrcfiles = '../v4/'+ifield+'_combined_'+['exposures.fits.gz','chips.fits.gz','allsrc.fits.gz','allobj.fits.gz',$
                                        'expmap.fits.gz','allobj_xmatch.fits.gz','allobj_bright.fits']
  destfiles = v5dir+ifield+'_combined_'+['exposures.fits.gz','chips.fits.gz','allsrc.fits.gz','allobj.fits.gz',$
                                        'expmap.fits.gz','allobj_xmatch.fits.gz','allobj_bright.fits']

  bd1 = where(file_test(srcfiles) eq 0,nbd1)
  if nbd1 gt 0 then stop,strtrim(nbd1,2)+' source files are missing'

  bd2 = where(file_test(destfiles) eq 1,nbd2)
  ;if nbd2 gt 0 then stop,strtrim(nbd2,2)+' destination files already exist'
  print,strtrim(nbd2,2)+' destination files already exist'
  if nbd2 gt 0 then file_delete,destfiles[bd2]

  file_link,relsrcfiles,destfiles,/allow

  ;stop

endfor

stop

end
