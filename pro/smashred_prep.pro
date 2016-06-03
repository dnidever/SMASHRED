pro smashred_prep,dir,nmulti=nmulti,keepstandards=keepstandards,refname=refname

; Get "calibrated" images ready for PHOTRED

; STEPS:
; -rename tu to c4d
; -move standards
; -make PHOTRED-ready split files (using c4d image files)
; -stack the deep exposures
; -group files into fields (short/long) and rename them (using split files)
; -create PHOTRED WCS input file list (using renamed split files)
; -WCS (using renamed split files)

; Set /keepstandards to specify that you DON'T want to move the standards.
  
; SETTINGS
;----------
dodeepstack = 0   ; Are we stacking the deep frames or not
sepfielddir = 1   ; Fields in separate directories
if n_elements(nmulti) eq 0 then nmulti=10

;goto,movesepfielddir


; Start the logfile
;------------------
; format is photred.DATETIME.log
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
smonth = strtrim(month,2)
if month lt 10 then smonth = '0'+smonth
sday = strtrim(day,2)
if day lt 10 then sday = '0'+sday
syear = strmid(strtrim(year,2),2,2)
shour = strtrim(hour,2)
if hour lt 10 then shour='0'+shour
sminute = strtrim(minute,2)
if minute lt 10 then sminute='0'+sminute
ssecond = strtrim(round(second),2)
if second lt 10 then ssecond='0'+ssecond
logfile = 'smashred.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile


t0 = systime(1)

if n_elements(dir) eq 0 then cd,current=dir

; CD to the desired directory
cd,current=origdir
cd,dir[0]


; Step 1. Rename "tu" files to "c4d" files if necessary
;------------------------------------------------------
print,'Step 1. Rename "tu" files to "c4d" files if necessary'
files = file_search('tu*.fits*',count=nfiles)
if nfiles gt 0 then SMASHRED_PREP_FILEINFO,fstr,files=files

; Rename tuXXXXXX.fits.fz to the new DSID names
; DTNSANAM= 'c4d_140108_024739_ori.fits'  /  file name in NOAO Science Archive    
undefine,renamelines
for i=0,nfiles-1 do begin
  if strmid(fstr[i].file,0,2) eq 'tu' then begin
    outfile = file_basename(fstr[i].dsidname)
    if strmid(fstr[i].file,1,2,/reverse) eq 'fz' and strmid(outfile,3,4,/reverse) eq 'fits' then outfile+='.fz'
    ; image, ooi
    ;if fstr[i].desext eq '0' then outfile = repstr(outfile,'ori','ooi')
    ;if (fstr[i].desext eq '0') and (fstr[i].bitpix eq -32) then outfile = repstr(outfile,'ori','ooi')
    if (fstr[i].desext eq '0') and ((fstr[i].bitpix eq -32) or (fstr[i].bitpix eq -64)) then outfile = repstr(outfile,'ori','ooi')  ; added -64 for 08/2013 data
    ; dqmap, ood
    ;if fstr[i].desext2 eq 'MASK' then outfile = repstr(outfile,'ori','ood')
    if (fstr[i].desext2 eq 'MASK') or (fstr[i].bitpix eq 32) then outfile = repstr(outfile,'ori','ood')
    ; weightmap, oow
    ;if fstr[i].desext2 eq 'WEIGHT' then outfile = repstr(outfile,'ori','oow')
    if (fstr[i].desext2 eq 'WEIGHT' and fstr[i].bitpix eq -32) then outfile = repstr(outfile,'ori','oow')
    ; add filter name, c4d_140108_024739_ooi_g.fits
    outfile = repstr(outfile,'.fits','_'+fstr[i].filter+'.fits')
    if file_test(outfile) eq 0 then begin
      file_move,fstr[i].file,outfile,/verbose
      ;print,fstr[i].file,' -> ',outfile
      fstr[i].file = outfile
      push,renamelines,fstr[i].archivename+'  '+fstr[i].file
    endif else print,'CANNOT move ',fstr[i].file,' because output filename already exists'
  endif
endfor

if n_elements(renamelines) gt 0 then begin
  print,'Writing renaming to >>file_renaming.txt<<'
  ;writecol,'file_renaming.txt',fstr.archivename,'  '+fstr.file
  writeline,'file_renaming.txt',renamelines
endif else print,'No "tu" files to rename'


; Get info on all the c4d files
;-------------------------------
files = file_search('c4d_*_*.fits*',count=nfiles)
if nfiles eq 0 then begin
  print,'No ARCHIVE files.  Skipping Steps 2+3 and going to STACKING.'
  goto,dostack
endif

SMASHRED_PREP_FILEINFO,fstr,files=files
nfiles = n_elements(files)

;stop


; Step 2. Move standards into their own directory
;-------------------------------------------------
print,'Step 2. Move standards into their own directory'
if not keyword_set(keepstandards) then begin
  stdind = where( strmid(strupcase(fstr.object),0,4) eq 'SDSS' or $
                  strmid(strupcase(fstr.object),0,2) eq 'SA' or $
                  strmid(strupcase(fstr.object),0,6) eq 'MAXVIS' or $
                  stregex(fstr.object,'[0-9]{6}-[0-9]{6}',/boolean) eq 1 or $
                  strmid(strupcase(fstr.object),0,1) eq 'E',nstdind)
  if nstdind gt 0 then begin
    print,'Moving ',strtrim(nstdind,2),' standard star field images to standards/'
    writecol,-1,lindgen(nstdind)+1,'  '+fstr[stdind].file,'  '+fstr[stdind].object
    if file_test('standards/') eq 0 then file_mkdir,'standards'
    file_move,fstr[stdind].file,'standards/',/verbose
    stdfstr = fstr[stdind]
    remove,stdind,fstr
  endif else print,'No standard star field images to move.'
endif else print,'Keeping standards'
print,''
nfstr = n_elements(fstr)

; Move Junk
; we should move the 30s r cal/focus exposures and junk someplace else
;jnkind = where( strtrim(fstr.object,2) eq '' or fstr.exptime lt 10 or $
;                (fstr.filter eq 'r' and fstr.exptime eq 30),njnkind)
;if njnkind gt 0 then begin
;  print,'Moving ',strtrim(njnkind,2),' junk to junk/'
;  writecol,-1,lindgen(njnkind)+1,'  '+fstr[jnkind].file,'  '+fstr[jnkind].object
;  if file_test('junk/') eq 0 then file_mkdir,'junk'
;  file_move,fstr[jnkind].file,'junk/',/verbose
;  jnkfstr = fstr[jnkind]
;  remove,jnkind,fstr
;endif
;print,''
;nfstr = n_elements(fstr)

print,'Check that everything is okay (i.e. remove test exposures).  Then type ".c" to continue.'

stop

; Step 3. Now make the PHOTRED-ready FITS files for the science exposures
;-------------------------------------------------------------------------
print,'Step 3. Make the PHOTRED-ready split chip files'

; get info for c4d image files
files = file_search('c4d_*_ooi*.fits*',count=nfiles)
if nfiles eq 0 then begin
  print,'No ARCHIVE image files to split.  Skipping to STACKING.'
  goto,dostack
endif
SMASHRED_PREP_FILEINFO,fstr,files=files
si = sort(fstr.expnum)  ; chronological order
fstr = fstr[si]
nfstr = n_elements(fstr)

undefine,cmd,cmddir
for i=0,nfstr-1 do begin
  if file_test(fstr[i].expnum+'_01.fits') eq 0 or keyword_set(redo) then begin
    push,cmd,"smashred_imprep_single,'"+file_basename(fstr[i].file)+"'"
  endif else push,cmd,''
endfor
todo = where(cmd ne '',ntodo)
ncmd = n_elements(cmd)
if ntodo gt 0 then begin
  cd,current=curdir
  cmddir = strarr(ncmd)+curdir
  print,'Running SMASHRED_IMPREP_SINGLE via PBS_DEAMON on ',strtrim(ncmd,2),' exposures'
  PBS_DAEMON,cmd[todo],cmddir[todo],nmulti=nmulti,prefix='imprep',/hyperthread,/idle,waittime=5 ; 30
endif

; Check if it worked and then move them to orig/ directory
if file_test('orig',/directory) eq 0 then file_mkdir,'orig'
for i=0,ntodo-1 do begin
  fstr1 = fstr[todo[i]]
  chipfiles = file_search(fstr1.expnum+'_??.fits',count=nchipfiles)
  ;if nchipfiles ge 60 then begin
  if nchipfiles ge 58 then begin  ; some exposures only have 58 for some reason
    print,strtrim(nchipfiles,2),' chip files for ',fstr1.expnum,'.  Moving ',fstr1.base,' to orig/'
    origfiles = file_search(fstr1.base+'*',count=norigfiles)
    if norigfiles gt 0 then FILE_MOVE,origfiles,'orig/',/over,/verbose
  endif else begin
    if nchipfiles eq 0 then print,'NO chip files for ',fstr1.expnum else $
      print,'ONLY ',strtrim(nchipfiles,2),' chip files for ',fstr1.expnum,'.  Need 60.'
  endelse
endfor

;stop

; Step 4. Stack the deep exposures
;----------------------------------
dostack:
if keyword_set(dodeepstack) then begin
  print,'Step 4. Stack the deep exposures'

  ; get split files, chip 1 only
  files = file_search('????????_01.fits',count=nfiles)
  if nfiles eq 0 then begin
    print,'No SPLIT images to stack/group/rename.  Skipping Step 4-6 and going to WCS.'
    goto,dowcs
  endif
  SMASHRED_PREP_FILEINFO,fstr,files=files

  ; stack the deep images
  SMASHRED_IMSTACK,fstr

endif else print,'NOT stacking the deep exposures'

;stop

; Step 5. Now figure out the fields and make shallow/deep groups
;----------------------------------------------------------------
dofields:
print,"Step 5. Group the exposures into 'fields' (shallow/deep) and rename"

; ***use the FSTR structure from smashred_imstack***

ui = uniq(fstr.object,sort(fstr.object))
fields = fstr[ui].object  ; unique fields
nfields = n_elements(fields)

fcount = 0
undefine,fieldstr  ; the names for each "field", shallow/deep are separate
for i=0,nfields-1 do begin
  ; Shallow exposures
  gshallow = where(fstr.object eq fields[i] and fstr.exptime lt 80,nshallow)
  if nshallow gt 0 then begin
    shname = 'F'+strtrim(fcount+1,2)
    fname = fields[i]+'sh'
    fstr[gshallow].shname = shname
    fstr[gshallow].fname = fname
    fstr[gshallow].newfile = shname+'-'+fstr[gshallow].expnum
    push,fieldstr,{shname:shname,name:fname}
    fcount++
  endif

  ; Deep stacks
  if keyword_set(dodeepstack) then begin
    ; Deep exposures, only rename COMBINED files
    ;gdeep = where(fstr.object eq fields[i] and fstr.exptime ge 80,ndeep)
    gdeep = where(fstr.object eq fields[i] and $
                ((fstr.deep eq 1 and fstr.combfile ne '') or (fstr.ncombine gt 1 and stregex(fstr.file,'c_[0-9][0-9].fits',/boolean) eq 1)) ,ndeep)
    if ndeep gt 0 then begin
      shname = 'F'+strtrim(fcount+1,2)
      fname = fields[i]
      fstr[gdeep].shname = shname
      fstr[gdeep].fname = fname
      fstr[gdeep].newfile = shname+'-'+fstr[gdeep].expnum
      push,fieldstr,{shname:shname,name:fname}
      fcount++
    endif
  ; NO deep stacks
  endif else begin
    ; Deep exposures, only rename COMBINED files
    gdeep = where(fstr.object eq fields[i] and fstr.exptime ge 80,ndeep)
    if ndeep gt 0 then begin
      shname = 'F'+strtrim(fcount+1,2)
      fname = fields[i]
      fstr[gdeep].shname = shname
      fstr[gdeep].fname = fname
      fstr[gdeep].newfile = shname+'-'+fstr[gdeep].expnum
      push,fieldstr,{shname:shname,name:fname}
      fcount++
    endif
  endelse
endfor
print,'There are ',strtrim(n_elements(fieldstr),2),' field groupings'
writecol,-1,fieldstr.shname+'  ',fieldstr.name
print,'Writing to "fields" file'
WRITELINE,'fields',fieldstr.shname+'   '+fieldstr.name

;stop

; Now move/rename the files
if not keyword_set(keepstandards) then begin
  print,''
  print,'Renaming the files'
  for i=0,nfstr-1 do begin
    ; This should get renamed
    if fstr[i].newfile ne '' then begin
      files1 = file_search(fstr[i].expnum+'_*.fits',count=nfiles1)
      if nfiles1 gt 0 then begin
        print,strtrim(i+1,2),'/',strtrim(nfstr,2),'  ',fstr[i].expnum,' -> ',fstr[i].newfile
        newfiles1 = fstr[i].shname+'-'+files1
        ;writecol,-1,files1,'  '+newfiles1
        file_move,files1,newfiles1
      endif else print,'NO files for ',fstr[i].expnum
    endif
  endfor

; don't rename standards
endif else begin
  print,'Standards. Not renaming'
  fstr.newfile = fstr.expnum  ; we are using the original filenames
endelse

save,fstr,fields,fieldstr,file='smashred_prep.dat'

;stop

; Move short dithers into separate directory
;shortdithers:
;restore,'smashred_prep.dat'
;add_tag,fstr,'radeg',0.0d0,fstr
;add_tag,fstr,'decdeg',0.0d0,fstr
;add_tag,fstr,'dist',0.0,fstr
;add_tag,fstr,'toremove',0,fstr
;fstr.radeg = 15*sexig2ten(fstr.ra)
;fstr.decdeg = sexig2ten(fstr.dec)
;ui = uniq(fstr.object,sort(fstr.object))
;fields = fstr[ui].object  ; unique fields
;nfields = n_elements(fields)
;bands = ['u','g','r','i','z']
;nbands = n_elements(bands)
;undefine,toremove
;for i=0,nfields-1 do begin
;  print,fields[i]
;  ind = where(fstr.object eq fields[i] and fstr.exptime gt 100,nind)
;
;  ; Get the deep images for all bands and find the mean RA and DEC
;  mndeepra = mean(fstr[ind].radeg)
;  mndeepdec = mean(fstr[ind].decdeg)
;
;  for j=0,nbands-1 do begin
;
;    ; For each band find the short exposure that is closest do 
;    ; the mean deep RA/DEC
;    shind = where(fstr.object eq fields[i] and fstr.filter eq bands[j] and fstr.exptime lt 100 and fstr.exptime gt 31,nshind)
;    dist = sphdist(fstr[shind].radeg,fstr[shind].decdeg,mndeepra,mndeepdec,/deg)
;    fstr[shind].dist = dist
;    si = reverse(sort(dist))
;    if nshind gt 1 then fstr[shind[si[0:nshind-2]]].toremove=1
;    ;if nshind gt 1 then push,toremove,shind[si[0:nshind-2]]
;  endfor
;  allsh = where(fstr.object eq fields[i] and fstr.exptime lt 100 and fstr.exptime gt 31,nallsh)
;  writecol,-1,fstr[allsh].expnum+' ',fstr[allsh].filter,fstr[allsh].exptime,fstr[allsh].dist,fstr[allsh].toremove
;
;  stop
;endfor


; Step 6: Get astrometic catalogs
;--------------------------------
dowcs:
print,'Step 6: Get astrometic catalogs'
;restore,'smashred_prep.dat'
fieldstr = importascii('fields',fieldnames=['SHNAME','NAME'],fieldtypes=[7,7])
nfields = n_elements(fieldstr)

; get renamed split files, chip 1 only
if not keyword_set(keepstandards) then begin
  files = file_search('F*-*_01.fits',count=nfiles)
  if nfiles eq 0 then begin
   print,'No RENAMED SPLIT/STACKED images to run WCS on.  Quitting.'
     return
  endif
endif else files=file_search('*_01.fits',count=nfiles)
SMASHRED_PREP_FILEINFO,fstr,files=files
if keyword_set(keepstandards) then fstr.newfile=fstr.expnum

;ui = uniq(fstr.object,sort(fstr.object))
;fieldstr = fstr[ui]  ; unique fields
;nfields = n_elements(fieldstr)


; Reference catalog
if n_elements(refname) eq 0 then begin
  if not keyword_set(keepstandards) then refcatname='USNO-B1' else refcatname='2MASS-PSC'
endif else refcatname=refname


;for i=0,nfields-1 do begin
; This does all of them
  SMASHRED_WCSPREP,fstr,refcatname=refcatname
  ; maybe use pbs_daemon as well.
;endfor


; Step 7: Move to separate field directories
;-------------------------------------------
movesepfielddir:
if keyword_set(sepfielddir) and not keyword_set(keepstandards) then begin
  print,'Step 7.  Moving files to separate field directories'
  readline,'fields',fieldlines
  dum = strsplitter(fieldlines,' ',/extract)
  shnames = reform(dum[0,*])
  fieldnames = reform(dum[1,*])
  nfields = n_elements(shnames)
  print,strtrim(nfields,2),' Fields'
  for i=0,nfields-1 do begin
    print,strtrim(i+1,2),' ',shnames[i],' ',fieldnames[i]
    if file_test(shnames[i],/directory) eq 0 then file_mkdir,shnames[i]
    fieldfiles = file_search(shnames[i]+'-*_*',count=nfieldfiles)
    print,'Moving ',strtrim(nfieldfiles,2),' files'
    if nfieldfiles gt 0 then file_move,fieldfiles,shnames[i]
  endfor
endif


; Step 8: Create PHOTRED file list
;----------------------------------
print,'Step 8. Creating PHOTRED WCS.inlist file'
if file_test('logs',/directory) eq 0 then file_mkdir,'logs'
if keyword_set(sepfielddir) and not keyword_set(keepstandards) then begin
  files = file_search('F*/F*-*_??.fits',/fully)
endif else begin
  files = file_search('F*-*_??.fits',/fully)
endelse
writeline,'logs/WCS.inlist',files



print,'SMASHRED_PREP FINISHED.'
dt = systime(1)-t0
print,'dt = ',strtrim(dt,2),' sec.'


; End logfile
;------------
JOURNAL

stop

cd,origdir

end
