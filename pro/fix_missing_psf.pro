pro fix_missing_psf,nmulti=nmulti

; Rerun daophot on the chips with missing PSF files

; 1) Find the missing PSF files
; 2) Make the "missingpsf/" directory
; 3) Copy the necessary files to missingpsf/
; 4) Run daophot.sh with pbs_daemon
; 5) Remove temporary files
; 6) Copy psf files to individual directories
; 7) Fix summary.fits files, one in individual directory and main directory

cd,current=curdir
print,'Running fix_missing_psf in ',curdir

if n_elements(nmulti) eq 0 then nmulti=15

; 1) Find the missing PSF files
print,'1) Find the missing PSF files'
alsfiles = file_search('F*/F*-*_??.als',count=nalsfiles)
allfield = file_dirname(alsfiles)
allbase = file_basename(alsfiles,'.als')
allpsffiles = allfield+'/'+allbase+'.psf'
missind = where(file_test(allpsffiles) eq 0,nmissing)
print,strtrim(nmissing,2),' psf files are missing'
if nmissing eq 0 then return
field = allfield[missind]
base = allbase[missind]

;goto,step5

; 2) Make the "missingpsf/" directory
print,'2) Make the "missingpsf/" directory'
if file_test('missingpsf',/directory) eq 0 then file_mkdir,'missingpsf'

; 3) Copy the necessary files to missingpsf/
print,'3) Copy the necessary files to missingpsf/'
for i=0,nmissing-1 do begin
  file_copy,field[i]+'/'+base[i]+'*','missingpsf/',/over
  ;file_copy,field[i]+'/'+base[i]+['.fits','.opt','.als.opt','.cmn.lst'],'missingpsf/'
endfor
; copy the scripts
file_copy,field[0]+'/'+['daophot.sh','apcor.opt','photo.opt','goodpsf.pro','lstfilter'],'missingpsf/',/verbose,/over

stop

; 4) Run daophot.sh with pbs_daemon
print,'4) Run daophot.sh with pbs_daemon'
cmd = './daophot.sh '+base
dirs = strarr(nmissing)+curdir+'/missingpsf/'
pbs_daemon,cmd,dirs,/hyper,nmulti=nmulti,waittime=1,/cdtodir

; 5) Check that the als files are the same as the originals
step5:
print,'5) Check that the als files are the same as the originals'
for i=0,nmissing-1 do begin
  origals = field[i]+'/'+base[i]+'.als'
  newals = 'missingpsf/'+base[i]+'.als'
  spawn,['diff','origals','newals'],out,errout,/noshell
  if out[0] ne '' then print,origals,' ',newals,' NOT THE SAME'
endfor

; 6) Remove temporary files
print,'6) Remove temporary files'
files = file_search('missingpsf/F*-*_??*')
b = where(stregex(files,'.psf$',/boolean) eq 1,nb)  ; the files we want to keep
if nb gt 0 then remove,b,files
file_delete,files

; 7) Copy psf files to individual directories
print,'7) Copy psf files to individual directories'
for i=0,nmissing-1 do file_copy,'missingpsf/'+base[i]+'.psf',field[i],/verbose

; 8) Fix summary.fits files, one in individual directory and main directory
print,'8) Fix summary.fits files'
ui = uniq(field,sort(field))
ufield = field[ui]
nufield = n_elements(ufield)
for i=0,nufield-1 do begin
  print,strtrim(i+1,2),' ',ufield[i]
  sumfile = file_search(ufield[i]+'/*_summary.fits')
  sumfile = sumfile[0]
  file_copy,sumfile,sumfile+'.bak',/over
  print,'Fixing ',sumfile
  head = headfits(sumfile,exten=0)
  expstr = mrdfits(sumfile,1)
  chipstr = mrdfits(sumfile,2)

  ; Loop through chips for this field
  chipind = where(field eq ufield,nchipind)
  for j=0,nchipind-1 do begin
    ; Load DAOPHOT PSF file
    psffile = 'missingpsf/'+base[chipind[j]]+'.psf'
    READLINE,psffile,psflines
    ; PENNY1     69    4    6    0   14.048       1091.621   1022.5   2046.5
    psfarr = strsplit(psflines[0],' ',/extract)
    chipstr_ind = where(chipstr.base eq base[chipind[j]],nchipstr_ind)
    chipstr[chipstr_ind].dao_psftype = psfarr[0]
    chipstr[chipstr_ind].dao_psfboxsize = long(psfarr[1])
  endfor

  file_delete,sumfile
  print,'Rewriting ',sumfile
  mwrfits,0,sumfile,head,/create
  mwrfits,expstr,sumfile,/silent
  mwrfits,chipstr,sumfile,/silent

  ; copy to main directory
  file_copy,sumfile,'.',/over
  ;stop
endfor

print,'Finished'

stop

end
