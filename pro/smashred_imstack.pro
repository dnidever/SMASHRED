pro smashred_imstack,fstr,redo=redo

; Stack the deep SMASH images

if n_elements(fstr) eq 0 then begin
  print,'Syntax - smashred_imstack,fstr'
  return
endif

cd,current=curdir

; Copy photo.opt and daofindphot.sh
if file_test('photo.opt') eq 0 then file_copy,'~/daophot/photo.opt','.'
if file_test('daofindphot.sh') eq 0 then file_copy,'~/daophot/daofindphot.sh','.'

ui = uniq(fstr.object,sort(fstr.object))
fieldstr = fstr[ui]  ; unique fields 

bands = ['u','g','r','i','z']
nbands = n_elements(bands)

undefine,cmd,cmddir,cmdfile,cmdcombfile,cmdfcombfile,cmdindivfiles

; Fields loop
nfields = n_elements(fieldstr)
for i=0,nfields-1 do begin
  print,strtrim(i+1,2),'/',strtrim(nfields,2),'  ',fieldstr[i].object

  ; Bands loop
  for j=0,nbands-1 do begin
    ind = where(fstr.object eq fieldstr[i].object and fstr.filter eq bands[j] and fstr.exptime ge 100 and $
                stregex(fstr.file,'c_[0-9][0-9].fits',/boolean) eq 0,nind)  ; can't already be a combined image
    stackind = where(fstr.object eq fieldstr[i].object and fstr.filter eq bands[j] and fstr.exptime ge 100 and $
                     stregex(fstr.file,'c_[0-9][0-9].fits',/boolean) eq 1,nstackind)
    if nind eq 0 then begin
      txt = 'No deep '+bands[j]+' exposures for '+fieldstr[i].object
      if nstackind gt 0 then txt+='.  Stack already exists '+strjoin(fstr[stackind].expnum,' ')
      print,txt
      goto,BOMB1
    endif
    fstr1 = fstr[ind]
    fstr[ind].deep = 1
    fstr[ind[0]].combfile = fstr[ind[0]].expnum+'c'  ; that's the final combined name (renamed below)

    chipfiles = file_search(fstr1[0].expnum+'_??.fits',count=nchipfiles)
    ext = strmid(file_basename(chipfiles,'.fits'),1,2,/reverse_offset)

    ; Chip loop
    for k=0,nchipfiles-1 do begin
      outfile1 = fstr1[0].expnum+'_'+ext[k]+'_comb.fits'
      outfile2 = fstr1[0].expnum+'c_'+ext[k]+'.fits'
      if (file_test(outfile1) eq 0 and file_test(outfile2) eq 0) or keyword_set(redo) then begin
        push,cmd,"imstack,['"+strjoin(fstr1.expnum+"_"+ext[k]+".fits","','")+"']"
      endif else begin
        push,cmd,''
        if file_test(outfile1) then print,outfile1,' already exists and /redo not set'
        if file_test(outfile2) then print,outfile2,' already exists and /redo not set'
      endelse
      push,cmddir,curdir
      push,cmdfile,fstr1[0].expnum+'_'+ext[k]
      push,cmdcombfile,outfile1
      push,cmdfcombfile,outfile2
      push,cmdindivfiles,strjoin(fstr1.expnum+'_'+ext[k]+'.fits',',')
    endfor

    BOMB1:

  endfor ; band loop
endfor ; field loop

;stop

; Run IMSTACK with PBS_DEAMON
;----------------------------
ncmd = n_elements(cmd)
if ncmd gt 0 then todo=where(cmd ne '',ntodo) else ntodo=0
if ntodo gt 0 then print,'Running IMSTACK on images'
print,strtrim(ntodo,2),' files to run IMSTACK on'

if ntodo gt 0 then begin
  ; Submit the jobs to the daemon
  nmulti = 10
  PBS_DAEMON,cmd[todo],cmddir[todo],nmulti=nmulti,prefix='imstk',/hyperthread,/idle,waittime=1  ;2, 5
endif


; Create directory for individual exposures if necessary
if file_test('deepexp',/directory) eq 0 then file_mkdir,'deepexp'

; Check for success/failures, and move individual exposures
for i=0,ncmd-1 do begin
  test1 = file_test(cmdcombfile[i])
  test2 = file_test(cmdfcombfile[i])

  ; Move individual images and rename combined file
  if test1 eq 1 then begin
    indivfiles = strtrim(strsplit(cmdindivfiles[i],',',/extract),2)
    tomv = where(file_test(indivfiles),ntomv)
    if ntomv gt 0 then begin
      FILE_MOVE,indivfiles[tomv],'deepexp'
      print,'Moving ',strjoin(indivfiles[tomv],' '),' to deepexp/'
    endif
    ; rename the combined file
    ; 00272326_01_comb.fits -> 00272326c_01.fits
    dum = strsplit(combfile[i],'_',/extract)
    newfile = dum[0]+'c_'+dum[1]+'.fits'
    FILE_MOVE,combfile[i],newfile,/over,/verbose
    ; Rename BPM file
    bpmfile = file_basename(combfile[i],'.fits')+'.bpm.fits'
    newbpmfile = dum[0]+'c_'+dum[1]+'.bpm.fits'
    if file_test(bpmfile) then FILE_MOVE,bpmfile,newbpmfile,/over,/verbose
  endif else begin
    if test1 eq 0 and test2 eq 0 then print,cmdcombfile[i],' and ',cmdfcombfile[i],' NOT FOUND'
    ;print,combfile[i],' NOT FOUND'
  endelse
endfor

;stop

end
