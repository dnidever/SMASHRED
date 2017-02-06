pro make_allobj_deeponly,version,outversion,redo=redo

;; Remake the allobj SMASH catalogs with only the
;; the deep exposures.

rootdir = smashred_rootdir()+'cp/red/photred/'
if n_elements(version) eq 0 then version='v5'
if n_elements(outversion) eq 0 then outversion=version

; Do this for ALL the fields
files = file_search(rootdir+'catalogs/final/'+version+'/Field*_combined_allobj.fits.gz',count=nfiles)
print,strtrim(nfiles,2),' allobj files to fix'

for i=0,nfiles-1 do begin

  ifield = file_basename(files[i],'_combined_allobj.fits.gz')
  print,strtrim(i+1,2),' ',ifield

  outfile = rootdir+'catalogs/final/'+outversion+'/'+ifield+'_combined_allobj_deep.fits'
  if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
    print,outfile,' exists and /redo not set'
    goto,BOMB
  endif

  fstr = mrdfits(rootdir+'catalogs/final/'+version+'/'+ifield+'_combined_exposures.fits.gz',1)
  chstr = mrdfits(rootdir+'catalogs/final/'+version+'/'+ifield+'_combined_chips.fits.gz',1)
  allobj = mrdfits(rootdir+'catalogs/final/'+version+'/'+ifield+'_combined_allobj.fits.gz',1)
  allsrc = mrdfits(rootdir+'catalogs/final/'+version+'/'+ifield+'_combined_allsrc.fits.gz',1)

  ; Redo the average photometry, morphology parameters and coordinates
  SMASHRED_AVERAGEPHOT,fstr,chstr,allsrc,allobj,/deeponly
  SMASHRED_AVERAGEMORPHCOORD,fstr,chstr,allsrc,allobj,/deeponly

  ; Write out new allobj file
  print,'Writing to ',outfile
  file_delete,[outfile,outfile+'.gz'],/allow
  MWRFITS,allobj,outfile,/create
  spawn,['gzip',outfile],out,errout,/noshell

  ;stop
  BOMB:

endfor

stop

end
