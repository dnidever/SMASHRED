pro run_sex,field,redo=redo

; Run sextractor on combined images

rootdir = smashred_rootdir()+'cp/red/photred/'

if n_elements(field) eq 0 then begin
  print,'Syntax - run_sex,field,redo=redo'
  return
endif

print,'Running Source Extractor on field ',field

fdir = rootdir+'deep/'+field+'/'
; Make the main directory read/write-able
file_chmod,fdir,'755'o

; Load the "fields" file
READCOL,fdir+'fields',shnames,lnames,format='A,A',/silent
ind = where(strtrim(lnames,2) eq field,nind)
if nind eq 0 then begin
  print,field,' not found in "fields" file'
  return
endif
dir = fdir+shnames[ind[0]]+'/'
file_chmod,dir,'755'o  ; make writeable

; Make the "sex/" directory
sdir = dir+'sex/'
if file_test(sdir,/directory) eq 0 then file_mkdir,sdir
file_chmod,sdir,'755'o

; Copy the default files
file_copy,'/datalab/users/dnidever/nsc/instcal/default.'+['config','param','conv','nnw'],sdir,/allow,/over
cd,current=curdir
cd,sdir

; Fix default.config
READLINE,'default.config',clines
badlines = where(stregex(clines,'WEIGHT_IMAGE',/boolean) eq 1 or stregex(clines,'WEIGHT_TYPE',/boolean) eq 1 or $
                 stregex(clines,'FLAG_IMAGE',/boolean) eq 1 or stregex(clines,'FLAG_TYPE',/boolean) eq 1,nbadlines)
if nbadlines gt 0 then remove,badlines,clines
WRITELINE,'default.config',clines
; Fix default.param
READLINE,'default.param',plines
badlines = where(stregex(plines,'^IMAFLAGS_',/boolean) eq 1 or stregex(plines,'^NIMAFLAGS_',/boolean) eq 1,nbadlines)
if nbadlines gt 0 then plines[badlines]='#'+plines[badlines]
WRITELINE,'default.param',plines

; Get reference image from .dered filename
dfile = file_search(dir+'F*-*.dered',count=ndfile)
if ndfile gt 1 then begin
  print,'More than 1 dered file found'
  stop
endif
refname = file_basename(dfile,'.dered')

files = file_search(dir+refname+'_??_comb.fits',count=nfiles)
if nfiles gt 62 then begin
  print,strtrim(nfiles,2),' files.  Should be 60-62.'
  stop
endif

; Loop through the chips
for i=0,nfiles-1 do begin
  dir1 = file_dirname(files[i])+'/'
  base = file_basename(files[i],'.fits')
  head = headfits(files[i])
  print,strtrim(i+1,2),' ',base

  if file_test(base+'_cat.fits') eq 1 and not keyword_set(redo) then begin
    print,base+'_cat.fits exists already and /redo not set'
    goto,BOMB
  endif

  ; Read the .opt file
  READLINE,dir1+base+'.opt',optlines
  optarr = strsplitter(optlines,' ',/extract)
  g = where(stregex(optlines,'HI =',/boolean) eq 1,ng)
  satlevel = optarr[2,g[0]]
  g = where(stregex(optlines,'GA =',/boolean) eq 1,ng)
  gain = optarr[2,g[0]]
  g = where(stregex(optlines,'FW =',/boolean) eq 1,ng)
  fwhm = optarr[2,g[0]]
  scale = 0.27

  ; Load the the default.config file and modify
  readline,'default.config',sexlines
  ; Modify for this file
  sexlines2 = sexlines
  ; CATALOG_NAME
  g = where(stregex(sexlines2,'^CATALOG_NAME',/boolean) eq 1,ng)
  ;catfile = base+'.cat'
  catfile = base+'_cat.fits'
  sexlines2[g[0]] = 'CATALOG_NAME    '+catfile+' # name of the output catalog'
  ; SATUR_LEVEL
  g = where(stregex(sexlines2,'^SATUR_LEVEL',/boolean) eq 1,ng)
  sexlines2[g[0]] = 'SATUR_LEVEL     '+satlevel+'         # level (in ADUs) at which arises saturation'
  ; GAIN
  g = where(stregex(sexlines2,'^GAIN',/boolean) eq 1,ng)
  sexlines2[g[0]] = 'GAIN            '+gain+'             # detector gain in e-/ADU.'
  ; PIXEL_SCALE
  g = where(stregex(sexlines2,'^PIXEL_SCALE',/boolean) eq 1,ng)
  sexlines2[g[0]] = 'PIXEL_SCALE     '+strtrim(scale,2)+'             # size of pixel in arcsec (0=use FITS WCS info).'
  ; SEEING_FWHM
  g = where(stregex(sexlines2,'^SEEING_FWHM',/boolean) eq 1,ng)
  fwhmas = float(fwhm)*float(scale)
  sexlines2[g[0]] = 'SEEING_FWHM     '+strtrim(fwhmas,2)+'            # stellar FWHM in arcsec'
  ; MASK/WEIGHT file
  maskfile = base+'.mask.fits'
  PUSH,sexlines2,'WEIGHT_IMAGE  '+maskfile
  PUSH,sexlines2,'WEIGHT_TYPE   MAP_WEIGHT'
  ; Write the file
  sexconfigfile = base+'.config'
  WRITELINE,sexconfigfile,sexlines2

  ; Make symlinks
  file_link,'../'+file_basename(files[i]),sdir,/allow
  file_link,'../'+base+'.mask.fits',sdir,/allow

  ; Run sextractor
  print,'  Running SExtractor on ',base,'.fits'
  spawn,'sex -c '+sexconfigfile+' '+base+'.fits',out,errout

  hd = headfits(catfile,exten=1)
  ncat = sxpar(hd,'naxis1')
  print,'  ',strtrim(ncat,2),' sources found'

  ;stop
  BOMB:

endfor

; Load all of the files
base = file_basename(files,'.fits')
catfiles = base+'_cat.fits'
for i=0,nfiles-1 do begin
  base1 = file_basename(files[i],'_comb.fits')
  dum = strsplit(base1,'_',/extract)
  ccdnum = long(first_el(dum,/last))
  cat = mrdfits(catfiles[i],2,/silent)
  add_tag,cat,'ccdnum',0,cat
  cat.ccdnum = ccdnum
  push,all,cat
endfor
nall = n_elements(all)
print,strtrim(nall,2),' final sources'

ref = first_el(strsplit(file_basename(files[0],'_comb.fits'),'_',/extract))
;print,'Writing combined catalog to ',ref,'_cat.fits'
MWRFITS,all,ref+'_cat.fits',/create
print,'Writing combined catalog to ',field,'_sex.fits'
MWRFITS,all,fdir+field+'_sex.fits',/create

cd,curdir

stop

end
