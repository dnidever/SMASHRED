pro create_catalogv6_symlinks

; Create symlinks in catalogs/final/v6/ to v5/ for those fields
; that didn't need to be recalibrated.

rootdir = smashred_rootdir()+'cp/red/photred/'

check = mrdfits(rootdir+'catalogs/final/v5/check_calibrated_v5.fits',1)
check.field = strtrim(check.field,2)
ncheck = n_elements(check)
smashred_getredinfo,info,/silent

; 26 fields were calibrated and/or got new data
newfieldnum = [115, 117, 118, 121, 123, 130, 156, 158, 161, 164, 166, 167, 168,$
               172, 183, 62, 63, 67, 69, 70, 74, 76, 85, 87, 92, 118]
newfields = 'Field'+strtrim(newfieldnum,2)
MATCH,check.field,newfields,ind1,ind2,/sort,count=nmatch
fields = check.field
remove,ind1,fields
nfields = n_elements(fields)

v5dir = rootdir+'catalogs/final/v5/'
v6dir = rootdir+'catalogs/final/v6/'
if file_test(v6dir,/directory) eq 0 then file_mkdir,v6dir
if file_test(v6dir+'ast',/directory) eq 0 then file_mkdir,v6dir+'ast'
if file_test(v6dir+'stars1',/directory) eq 0 then file_mkdir,v6dir+'stars1'

for i=0,nfields-1 do begin
  ifield = fields[i]
  print,strtrim(i+1,2),' ',ifield

  srcfiles = v5dir+ifield+'_combined_'+['exposures.fits.gz','chips.fits.gz','allsrc.fits.gz','allobj.fits.gz',$
                                        'expmap.fits.gz','allobj_xmatch.fits.gz','allobj_bright.fits']
  srcfiles = [srcfiles,v5dir+'ast/'+ifield+'_complete.fits.gz',v5dir+'stars1/'+ifield+'_allobj_stars.fits.gz',$
              v5dir+'stars1/'+ifield+'_allobj_deep_stars.fits.gz']
  relsrcfiles = '../v5/'+ifield+'_combined_'+['exposures.fits.gz','chips.fits.gz','allsrc.fits.gz','allobj.fits.gz',$
                                        'expmap.fits.gz','allobj_xmatch.fits.gz','allobj_bright.fits']
  relsrcfiles = [relsrcfiles,'../v5/ast/'+ifield+'_complete.fits.gz','../v5/stars1/'+ifield+'_allobj_stars.fits.gz',$
                 '../v5/stars1/'+ifield+'_allobj_deep_stars.fits.gz']
  destfiles = v6dir+ifield+'_combined_'+['exposures.fits.gz','chips.fits.gz','allsrc.fits.gz','allobj.fits.gz',$
                                        'expmap.fits.gz','allobj_xmatch.fits.gz','allobj_bright.fits']
  destfiles = [destfiles,v6dir+'ast/'+ifield+'_complete.fits.gz',v6dir+'stars1/'+ifield+'_allobj_stars.fits.gz',$
               v6dir+'stars1/'+ifield+'_allobj_deep_stars.fits.gz']

  bd1 = where(file_test(srcfiles) eq 0,nbd1,comp=gd1,ncomp=ngd1)
  if nbd1 gt 0 then print,strtrim(nbd1,2)+' source files are missing: ',srcfiles[bd1]

  bd2 = where(file_test(destfiles) eq 1,nbd2)
  ;if nbd2 gt 0 then stop,strtrim(nbd2,2)+' destination files already exist'
  print,strtrim(nbd2,2)+' destination files already exist'
  if nbd2 gt 0 then file_delete,destfiles[bd2]

  ; copy the source files that exist
  file_link,relsrcfiles[gd1],destfiles[gd1],/allow

  ;stop

endfor

stop

end
