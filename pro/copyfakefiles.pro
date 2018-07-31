pro copyfakefiles,deepdirs,version=version

; Copy the files needed to run FAKERED on the deep data
rootdir = smashred_rootdir()+'cp/red/photred/'
fakedir = rootdir+'addfakes/'
if n_elements(version) eq 0 then version='v6'

ndeepdirs = n_elements(deepdirs)
if ndeepdirs eq 0 then $
  deepdirs = file_search(rootdir+'deep/Field*',/test_directory,count=ndeepdirs)
deepfields = file_basename(deepdirs)
ndeep = ndeepdirs
print,strtrim(ndeep,2),' deep fields to copy FAKERED files for'

fakefield = 'T1'

; Loop through the deep fields
for i=0,ndeep-1 do begin
  ifield = deepfields[i]
  print,strtrim(i+1,2),' ',ifield

  ; Check if the addfakes directory already exists
  fakedir1 = fakedir+ifield+'/'
  if file_test(fakedir1,/directory) eq 1 then begin
    print,fakedir1,' already exists. Skipping it.'
    goto, BOMB
  endif

  chipfile = rootdir+'catalogs/final/'+version+'/'+ifield+'_combined_chips.fits.gz'
  if file_test(chipfile) eq 0 then begin
    print,chipfile,' NOT FOUND'
    goto,BOMB
  endif

  ; Get short field name
  deepdir1 = deepdirs[i]+'/'
  shdirs = file_search(deepdir1+'F*',/test_directory,count=nshdirs)
  if nshdirs gt 1 then stop,'More than ONE directory in '+deepdir1
  shfield = file_basename(shdirs[0])

  ; 1) Make directory
  ;   mkdir ../../addfakes/Field71/F1/
  print,' 1) Make directories'
  file_mkdir,fakedir1
  file_mkdir,fakedir1+shfield

  ; 2) Copy over the files
  ;   cp -a F1/*_01* ../../addfakes/Field71/F1/
  print,' 2) Copy over the 01 files'
  cpfiles = file_search(deepdir1+shfield+'/'+shfield+'-*_01*',count=ncpfiles)
  file_copy,cpfiles,fakedir1+shfield+'/'

  ; Run fix_fitsheaders.pro, otherwise ADDSTAR won't run
  print,' Running fix_fitsheaders'
  outfitsfiles = file_search(fakedir1+shfield+'/'+shfield+'-*_01.fits',count=noutfitsfiles)
  if noutfitsfiles eq 0 then $
    outfitsfiles = file_search(fakedir1+shfield+'/'+shfield+'-*_01.fits.fz',count=noutfitsfiles)
  FILE_CHMOD,outfitsfiles,'755'o
  FIX_FITSHEADERS,outfitsfiles

  ; 3) Copy over apcor.lst and photred.setup files
  ;   cp apcor.lst ../../addfakes/Field71/apcor.lst.orig
  ;   cp photred.setup ../../addfakes/Field71/
  print,' 3) copy over apcor.lst and photred.setup files'
  file_copy,deepdir1+'apcor.lst',fakedir1+'apcor.lst.orig'
  file_copy,deepdir1+'photred.setup',fakedir1+'photred.setup'

  ; 4) Copy some fakered files
  ;   cd ../../addfakes/Field71
  ;   cp ../Field100/fakered.setup .
  ;   cp ../Field100/fakered.batch .
  ;   cp ../Field100/nocalib.trans .
  print,' 4) Copy some fakered files'
  file_copy,fakedir+'Field100/'+['fakered.setup','fakered.batch','nocalib.trans'],fakedir1

  ; Update fakered.setup with important photred.setup settings
  ;  alftiletype, daopsfva, daofitradfwhm
  READLINE,fakedir1+'photred.setup',photlines
  READLINE,fakedir1+'fakered.setup',fakelines
  photopts = ['alftiletype','daopsfva','daofitradfwhm']
  nphotopts = n_elements(photopts)
  for j=0,nphotopts-1 do begin
    ind1 = where(stregex(photlines,photopts[j],/boolean,/fold_case) eq 1,nind1)
    if nind1 gt 0 then begin
      brkind = where(stregex(fakelines,'##### STAGES #####',/boolean,/fold_case) eq 1,nbrkind)
      ; Put the line before the stages
      if nbrkind gt 0 then begin
        fakelines1 = fakelines[0:brkind[0]-1]
        fakelines2 = fakelines[brkind[0]:*]
        fakelines = [fakelines1, photlines[ind1[0]], fakelines2]
      ; Put the line at the end
      endif else begin
        fakelines = [fakelines, photlines[ind1[0]]]
      endelse
    endif
  endfor  
  WRITELINE,fakedir1+'fakered.setup',fakelines

  ; 5) Make fields file
  ;   F1T1  Field71
  print,' 5) Make "fields" file'
  writeline,fakedir1+'fields',shfield+fakefield+'  '+ifield

  ; 6) Run addfakes.pro
  ;   addfakes,'F1'
  print,' 6) Run addfakes.pro.'
  cd,fakedir1
  ADDFAKES,shfield

  ; 7) Add files to inlist
  ;   mkdir logs
  ;   ls F1/F1T1-*_01.fits >> logs/DAOPHOT.inlist
  print,' 7) Add files to inlist'
  file_mkdir,'logs'
  fitsfiles = file_search(shfield+'/'+shfield+fakefield+'-*_01.fits',count=nfitsfiles)
  writeline,'logs/DAOPHOT.inlist',fitsfiles

  ; 8) Make sure fakered.setup is good
  ; 9) run fakered.batch
  ;   idlbatch fakered.batch
  ; 10) run getcomplete
  ;   getcomplete,'F1'

  ;stop

  BOMB:
endfor

stop

end
