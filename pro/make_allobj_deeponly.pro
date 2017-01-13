pro make_allobj_deeponly

;; Remake the v4 allobj SMASH catalogs with only the
;; the deep exposures.
;; Put the final catalogs inf v4c/

rootdir = smashred_rootdir()

; Do the LMC periphery fields FIRST
fields = 'Field'+strtrim([17,23,24,25,26,33,56,57,58,59,60,61,62,63,64,65,66,67,68,69,80,106,130,157,160,181],2)
files = rootdir+'cp/red/photred/catalogs/final/v4/'+fields+'_combined_allobj.fits.gz'
nfiles = n_elements(files)

;files = file_search(rootdir+'cp/red/photred/catalogs/final/v4/Field*_combined_allobj.fits.gz',count=nfiles)
print,strtrim(nfiles,2),' allobj files to fix'

for i=0,nfiles-1 do begin

  ifield = file_basename(files[i],'_combined_allobj.fits.gz')
  print,strtrim(i+1,2),' ',ifield
  fstr = mrdfits(rootdir+'cp/red/photred/catalogs/final/v4/'+ifield+'_combined_exposures.fits.gz',1)
  chstr = mrdfits(rootdir+'cp/red/photred/catalogs/final/v4/'+ifield+'_combined_chips.fits.gz',1)
  allobj = mrdfits(rootdir+'cp/red/photred/catalogs/final/v4/'+ifield+'_combined_allobj.fits.gz',1)
  allsrc = mrdfits(rootdir+'cp/red/photred/catalogs/final/v4/'+ifield+'_combined_allsrc.fits.gz',1)

  smashred_averagephot,fstr,chstr,allsrc,allobj,/deeponly
  smashred_averagemorphcoord,fstr,chstr,allsrc,allobj,/deeponly

  ; Write out new allobj file
  outfile = rootdir+'cp/red/photred/catalogs/final/v4c/'+ifield+'_combined_allobj.fits'
  print,'Writing to ',outfile
  file_delete,[outfile,outfile+'.gz'],/allow
  MWRFITS,allobj,outfile,/create
  spawn,['gzip',outfile],out,errout,/noshell

  ;stop

endfor

stop

end
